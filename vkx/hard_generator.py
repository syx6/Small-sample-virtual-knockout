from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd


def train_hard_constrained_generator(
    state_csvs: list[str | Path],
    ko_col: str,
    target_kos: list[str],
    prior_dir: str | Path,
    out_dir: str | Path,
    features: list[str] | None = None,
    samples_per_ko: int = 300,
    max_residual_fraction: float = 0.35,
    anchor_method: str = "vkx",
    epochs: int = 80,
    seed: int = 11,
) -> dict[str, pd.DataFrame | str]:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    frames = [pd.read_csv(path) for path in state_csvs]
    features = _common_features(frames, ko_col, features)
    residual_bank = _collect_residual_bank(frames, ko_col, features)
    residual_bank.to_csv(out / "generator_residual_bank.csv", index=False)

    baseline = _baseline_delta(frames[0], ko_col, target_kos, prior_dir, features, seed, anchor_method=anchor_method)
    model_status, sampler = _fit_residual_sampler(residual_bank[features].to_numpy(dtype=float), epochs=epochs, seed=seed)
    samples = _sample_hard_constrained(
        frame=frames[0],
        ko_col=ko_col,
        features=features,
        baseline=baseline,
        sampler=sampler,
        samples_per_ko=samples_per_ko,
        max_residual_fraction=max_residual_fraction,
        seed=seed,
    )
    intervals = _intervals(samples, features)
    metrics = _evaluate_samples(samples, frames[0], ko_col, features)
    samples.to_csv(out / "hard_generator_virtual_cells.csv", index=False)
    intervals.to_csv(out / "hard_generator_intervals.csv", index=False)
    metrics.to_csv(out / "hard_generator_metrics.csv", index=False)
    _plot_generator_outputs(samples, intervals, metrics, out)
    _write_report(out, state_csvs, model_status, features, target_kos, max_residual_fraction, anchor_method)
    return {"samples": samples, "intervals": intervals, "metrics": metrics, "model_status": model_status}


def _common_features(frames: list[pd.DataFrame], ko_col: str, explicit: list[str] | None) -> list[str]:
    if explicit:
        return explicit
    ignored = {ko_col, "cell_id", "dataset", "state", "batch", "cell_type"}
    feature_sets = []
    for frame in frames:
        feature_sets.append({col for col in frame.columns if col not in ignored and pd.api.types.is_numeric_dtype(frame[col])})
    common = set.intersection(*feature_sets)
    if not common:
        raise ValueError("No common numeric state features were found across state CSV files.")
    return sorted(common)


def _collect_residual_bank(frames: list[pd.DataFrame], ko_col: str, features: list[str]) -> pd.DataFrame:
    from .core import control_mask

    rows = []
    for dataset_id, frame in enumerate(frames):
        for ko, group in frame.loc[~control_mask(frame[ko_col])].groupby(ko_col, observed=True):
            values = group[features].to_numpy(dtype=float)
            if len(values) < 3:
                continue
            residual = values - np.nanmean(values, axis=0, keepdims=True)
            tmp = pd.DataFrame(residual, columns=features)
            tmp["source_dataset_id"] = dataset_id
            tmp["ko_target"] = str(ko)
            rows.append(tmp)
    if not rows:
        raise ValueError("No KO residuals could be collected. Need labelled perturbation states with at least a few cells per KO.")
    return pd.concat(rows, ignore_index=True)


def _baseline_delta(
    frame: pd.DataFrame,
    ko_col: str,
    target_kos: list[str],
    prior_dir: str | Path,
    features: list[str],
    seed: int,
    anchor_method: str,
) -> pd.DataFrame:
    method = anchor_method.lower().replace("-", "").replace("_", "")
    if method == "vkx":
        from .core import run_virtual_ko

        result = run_virtual_ko(
            frame=frame,
            ko_col=ko_col,
            holdouts=target_kos,
            prior_dir=prior_dir,
            features=features,
            dataset_name="hard constrained generator baseline",
            modality="state score table",
            representation="state scores",
            calibrate="auto",
            shape_calibrate="none",
            seed=seed,
        )
        return result.delta_table
    from .formal_benchmark import _predict_constrained_ensemble, _predict_prior_model

    if method in {"pls", "ridge"}:
        pred, status = _predict_prior_model(frame, ko_col, target_kos, prior_dir, features, model_type=method)
    elif method in {"ensemble", "constrainedensemble", "vkxensemble"}:
        pred, status = _predict_constrained_ensemble(frame, ko_col, target_kos, prior_dir, features, seed=seed)
    else:
        raise ValueError("anchor_method must be one of vkx, pls, ridge, or ensemble.")
    if pred.empty:
        raise ValueError(f"Could not build generator anchor using {anchor_method}: {status.get('reason')}")
    rows = []
    for _, row in pred.iterrows():
        out = {"ko_target": row["ko_target"], "prediction_source": f"{anchor_method}_anchor"}
        for feature in features:
            out[f"pred_delta_{feature}"] = row.get(f"pred_delta_{feature}", np.nan)
        rows.append(out)
    return pd.DataFrame(rows)


def _fit_residual_sampler(residuals: np.ndarray, epochs: int, seed: int):
    residuals = np.nan_to_num(residuals, nan=0.0)
    centered = residuals - residuals.mean(axis=0, keepdims=True)
    try:
        import torch
        from torch import nn
        from torch.utils.data import DataLoader, TensorDataset
    except Exception:
        return "pca_residual_fallback_no_torch", _fit_pca_sampler(centered, seed)

    class ResidualVAE(nn.Module):
        def __init__(self, dim: int, latent: int, hidden: int):
            super().__init__()
            self.encoder = nn.Sequential(nn.Linear(dim, hidden), nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU())
            self.mu = nn.Linear(hidden, latent)
            self.logvar = nn.Linear(hidden, latent)
            self.decoder = nn.Sequential(nn.Linear(latent, hidden), nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU(), nn.Linear(hidden, dim))

        def forward(self, x):
            h = self.encoder(x)
            mu = self.mu(h)
            logvar = self.logvar(h).clamp(-6, 4)
            z = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
            return self.decoder(z), mu, logvar

        def sample(self, n: int):
            z = torch.randn(n, self.mu.out_features)
            return self.decoder(z).detach().cpu().numpy()

    torch.manual_seed(seed)
    x = torch.tensor(centered.astype("float32"))
    dim = x.shape[1]
    model = ResidualVAE(dim=dim, latent=max(2, min(12, dim // 3)), hidden=max(32, min(128, dim * 2)))
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
    loader = DataLoader(TensorDataset(x), batch_size=min(128, max(16, len(x) // 4)), shuffle=True)
    for _ in range(max(1, epochs)):
        for (xb,) in loader:
            recon, mu, logvar = model(xb)
            recon_loss = ((recon - xb) ** 2).mean()
            kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
            loss = recon_loss + 0.015 * kl
            opt.zero_grad()
            loss.backward()
            opt.step()

    def sampler(n: int, rng: np.random.Generator) -> np.ndarray:
        del rng
        return model.sample(n)

    return "torch_residual_vae_hard_constrained", sampler


def _fit_pca_sampler(centered: np.ndarray, seed: int):
    rng = np.random.default_rng(seed)
    u, s, vt = np.linalg.svd(centered, full_matrices=False)
    k = max(1, min(12, vt.shape[0]))
    components = vt[:k]
    scales = s[:k] / max(1, np.sqrt(centered.shape[0] - 1))

    def sampler(n: int, local_rng: np.random.Generator | None = None) -> np.ndarray:
        use_rng = local_rng or rng
        z = use_rng.normal(size=(n, k)) * scales.reshape(1, -1)
        return z @ components

    return sampler


def _sample_hard_constrained(
    frame: pd.DataFrame,
    ko_col: str,
    features: list[str],
    baseline: pd.DataFrame,
    sampler,
    samples_per_ko: int,
    max_residual_fraction: float,
    seed: int,
) -> pd.DataFrame:
    from .core import control_mask

    rng = np.random.default_rng(seed)
    control = frame.loc[control_mask(frame[ko_col]), features].to_numpy(dtype=float)
    rows = []
    for _, row in baseline.iterrows():
        ko = str(row["ko_target"])
        delta = np.asarray([row.get(f"pred_delta_{feature}", np.nan) for feature in features], dtype=float)
        delta = np.nan_to_num(delta, nan=0.0)
        ctrl = control[rng.integers(0, len(control), size=samples_per_ko)]
        residual = sampler(samples_per_ko, rng)
        residual = _hard_bound_residual(residual, delta, max_residual_fraction)
        values = ctrl + delta.reshape(1, -1) + residual
        tmp = pd.DataFrame(values, columns=features)
        tmp["ko_target"] = ko
        tmp["state"] = "hard-constrained generator cells"
        tmp["prediction_source"] = "baseline_delta_plus_bounded_residual"
        rows.append(tmp)
    return pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()


def _hard_bound_residual(residual: np.ndarray, delta: np.ndarray, max_fraction: float) -> np.ndarray:
    residual = residual - residual.mean(axis=0, keepdims=True)
    delta_norm = np.linalg.norm(delta) + 1e-9
    norms = np.linalg.norm(residual, axis=1, keepdims=True) + 1e-9
    max_norm = max_fraction * delta_norm
    residual = residual * np.minimum(1.0, max_norm / norms)
    direction = delta.reshape(1, -1) / delta_norm
    projection = residual @ direction.T
    residual = residual - np.minimum(projection, 0.0) * direction
    residual = residual - residual.mean(axis=0, keepdims=True)
    return residual


def _intervals(samples: pd.DataFrame, features: list[str]) -> pd.DataFrame:
    rows = []
    for ko, group in samples.groupby("ko_target", observed=True):
        for feature in features:
            values = group[feature].to_numpy(dtype=float)
            rows.append(
                {
                    "ko_target": ko,
                    "feature": feature,
                    "q05": float(np.nanquantile(values, 0.05)),
                    "q50": float(np.nanquantile(values, 0.50)),
                    "q95": float(np.nanquantile(values, 0.95)),
                }
            )
    return pd.DataFrame(rows)


def _evaluate_samples(samples: pd.DataFrame, frame: pd.DataFrame, ko_col: str, features: list[str]) -> pd.DataFrame:
    from .core import control_mask

    rows = []
    control = frame.loc[control_mask(frame[ko_col]), features].mean().to_numpy(dtype=float)
    for ko, generated in samples.groupby("ko_target", observed=True):
        true = frame.loc[frame[ko_col].astype(str) == str(ko), features]
        if true.empty:
            continue
        true_delta = true.mean().to_numpy(dtype=float) - control
        pred_delta = generated[features].mean().to_numpy(dtype=float) - control
        rows.append(
            {
                "ko_target": ko,
                "mae": float(np.mean(np.abs(pred_delta - true_delta))),
                "direction_cosine": _cosine(true_delta, pred_delta),
                "r2": _r2(true_delta, pred_delta),
            }
        )
    return pd.DataFrame(rows)


def _plot_generator_outputs(samples: pd.DataFrame, intervals: pd.DataFrame, metrics: pd.DataFrame, out: Path) -> None:
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
    except Exception:
        _plot_generator_fallback(intervals, metrics, out)
        return
    if not metrics.empty:
        fig, axes = plt.subplots(1, 3, figsize=(13, 4.5), constrained_layout=True)
        for ax, metric, title in zip(axes, ["direction_cosine", "r2", "mae"], ["Direction", "R2", "MAE"]):
            sns.barplot(data=metrics, x=metric, y="ko_target", color="#4C78A8", ax=ax)
            ax.set_title(title)
            ax.set_ylabel("")
        fig.suptitle("Hard-constrained Generator Evaluation")
        fig.savefig(out / "01_hard_generator_metric_panel.png", bbox_inches="tight", dpi=300)
        plt.close(fig)
    if not intervals.empty:
        top = intervals.assign(width=intervals["q95"] - intervals["q05"]).sort_values("width", ascending=False).head(30)
        fig, ax = plt.subplots(figsize=(10, max(4, 0.25 * len(top) + 1.8)), constrained_layout=True)
        ax.hlines(y=np.arange(len(top)), xmin=top["q05"], xmax=top["q95"], color="#4C78A8", linewidth=2)
        ax.scatter(top["q50"], np.arange(len(top)), color="#E76F51", s=18, zorder=3)
        ax.set_yticks(np.arange(len(top)))
        ax.set_yticklabels([f"{row.ko_target} | {_short(row.feature)}" for row in top.itertuples()], fontsize=8)
        ax.axvline(0, color="0.3", linewidth=1)
        ax.set_title("Generator uncertainty intervals")
        ax.set_xlabel("generated state value")
        fig.savefig(out / "02_hard_generator_intervals.png", bbox_inches="tight", dpi=300)
        plt.close(fig)


def _plot_generator_fallback(intervals: pd.DataFrame, metrics: pd.DataFrame, out: Path) -> None:
    try:
        from PIL import Image, ImageDraw, ImageFont
    except Exception:
        return
    img = Image.new("RGB", (1100, 520), "white")
    draw = ImageDraw.Draw(img)
    font = ImageFont.load_default()
    draw.text((30, 25), "Hard-constrained Generator Evaluation", fill=(20, 20, 20), font=font)
    y = 75
    if metrics.empty:
        draw.text((30, y), "No true KO labels available for generator accuracy metrics.", fill=(60, 60, 60), font=font)
    else:
        for _, row in metrics.iterrows():
            draw.text(
                (30, y),
                f"{row['ko_target']}: direction {row['direction_cosine']:.2f}, R2 {row['r2']:.2f}, MAE {row['mae']:.3f}",
                fill=(30, 30, 30),
                font=font,
            )
            y += 30
    img.save(out / "01_hard_generator_metric_panel.png")

    if intervals.empty:
        return
    top = intervals.assign(width=intervals["q95"] - intervals["q05"]).sort_values("width", ascending=False).head(30)
    row_h = 22
    img2 = Image.new("RGB", (1200, 80 + row_h * len(top)), "white")
    draw2 = ImageDraw.Draw(img2)
    draw2.text((30, 25), "Generator uncertainty intervals", fill=(20, 20, 20), font=font)
    min_v = float(np.nanmin(top["q05"]))
    max_v = float(np.nanmax(top["q95"]))
    span = max(max_v - min_v, 1e-9)
    for i, row in enumerate(top.itertuples()):
        y = 65 + i * row_h
        label = f"{row.ko_target} | {_short(row.feature, 42)}"
        draw2.text((30, y), label, fill=(30, 30, 30), font=font)
        x1 = 470 + int((row.q05 - min_v) / span * 560)
        x2 = 470 + int((row.q95 - min_v) / span * 560)
        xm = 470 + int((row.q50 - min_v) / span * 560)
        draw2.line((x1, y + 8, x2, y + 8), fill=(76, 120, 168), width=3)
        draw2.ellipse((xm - 3, y + 5, xm + 3, y + 11), fill=(231, 111, 81))
    img2.save(out / "02_hard_generator_intervals.png")


def _write_report(
    out: Path,
    state_csvs: list[str | Path],
    status: str,
    features: list[str],
    target_kos: list[str],
    max_fraction: float,
    anchor_method: str,
) -> None:
    text = f"""# Hard-constrained Residual Generator

This generator keeps the selected baseline KO direction fixed and only learns single-cell residual variation around that direction.

## Training data

- State CSV files: {len(state_csvs)}
- Common state features: {len(features)}
- Target KOs: {', '.join(target_kos)}
- Anchor method: `{anchor_method}`
- Residual backend: `{status}`
- Hard residual bound: residual norm <= {max_fraction:.2f} x baseline delta norm

## Outputs

- `generator_residual_bank.csv`
- `hard_generator_virtual_cells.csv`
- `hard_generator_intervals.csv`
- `hard_generator_metrics.csv`
- `01_hard_generator_metric_panel.png`
- `02_hard_generator_intervals.png`

## Interpretation

This is not an unconstrained free generator. A generated cell is:

```text
control cell + selected baseline KO delta + bounded learned residual
```

The residual is centered and norm-bounded so it cannot overwrite the predicted KO direction.
"""
    (out / "hard_generator_report.md").write_text(text, encoding="utf-8")


def _r2(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    denom = float(np.sum((y_true - np.mean(y_true)) ** 2))
    if denom <= 1e-12:
        return np.nan
    return float(1.0 - np.sum((y_true - y_pred) ** 2) / denom)


def _cosine(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    denom = float(np.linalg.norm(y_true) * np.linalg.norm(y_pred))
    return float(np.dot(y_true, y_pred) / denom) if denom > 1e-12 else np.nan


def _short(value: str, max_len: int = 34) -> str:
    value = str(value)
    return value if len(value) <= max_len else value[: max_len - 3] + "..."
