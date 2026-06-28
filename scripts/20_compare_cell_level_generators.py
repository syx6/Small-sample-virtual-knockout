from __future__ import annotations

import math
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import wasserstein_distance
from sklearn.cross_decomposition import PLSRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
import torch
from torch import nn
from torch.utils.data import DataLoader, TensorDataset


GENE_RE = re.compile(r"^[A-Z0-9.-]+$")
FIG_DIR = Path("results/figures")
HOLDOUT_KOS = ["STAT1", "JAK2", "IFNGR2", "IRF1"]
MODEL_LABELS = {
    "Residual baseline": "Residual",
    "Conditional VAE": "VAE",
    "Flow matching": "Flow",
    "Diffusion": "Diffusion",
    "Guided Conditional VAE": "Guided VAE",
    "Guided Flow matching": "Guided Flow",
    "Guided Diffusion": "Guided Diffusion",
}
STATE_LABELS = {
    "control cells": "control",
    "true KO cells": "true KO",
    "Residual baseline cells": "Residual",
    "Conditional VAE cells": "VAE",
    "Flow matching cells": "Flow",
    "Diffusion cells": "Diffusion",
    "Guided Conditional VAE cells": "Guided VAE",
    "Guided Flow matching cells": "Guided Flow",
    "Guided Diffusion cells": "Guided Diffusion",
}


def setup_plot() -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="talk")
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    for font in ["Microsoft YaHei", "SimHei", "Noto Sans CJK SC", "Arial Unicode MS"]:
        if font in available_fonts:
            plt.rcParams["font.family"] = font
            break
    plt.rcParams["axes.unicode_minus"] = False
    plt.rcParams["figure.dpi"] = 140
    plt.rcParams["savefig.dpi"] = 300


def control_mask(values: pd.Series) -> pd.Series:
    text = values.astype(str).str.lower()
    return text.eq("ctrl") | text.str.contains("nt|control|non|safe|neg")


def split_ko(label: str) -> list[str]:
    text = str(label).replace("+", "_").replace("|", "_").replace(",", "_")
    genes = [part.strip().upper() for part in text.split("_") if part.strip()]
    return [gene for gene in genes if not gene.lower().startswith("nt") and gene.lower() != "ctrl"]


def parse_gmt(path: Path, include_term_gene: bool = False) -> list[tuple[str, set[str]]]:
    terms = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            term = parts[0]
            genes = {gene.upper() for gene in parts[2:] if GENE_RE.match(gene.upper())}
            if include_term_gene:
                first = term.split()[0].upper() if term.split() else ""
                if GENE_RE.match(first):
                    genes.add(first)
            if genes:
                terms.append((term, genes))
    return terms


def select_prior_terms(priors_dir: Path, perturb_genes: set[str], max_terms_per_library: int = 160) -> list[tuple[str, set[str]]]:
    selected = []
    for path in sorted(priors_dir.glob("*.gmt")):
        scored = []
        for term, genes in parse_gmt(path, include_term_gene=path.stem == "ppi_hub"):
            overlap = len(genes & perturb_genes)
            if overlap == 0 or len(genes) < 5 or len(genes) > 800:
                continue
            scored.append(((overlap, -len(genes)), f"{path.stem}:{term}", genes))
        scored.sort(reverse=True, key=lambda item: item[0])
        selected.extend((name, genes) for _, name, genes in scored[:max_terms_per_library])
    return selected


def ko_prior_vector(label: str, terms: list[tuple[str, set[str]]]) -> np.ndarray:
    genes = set(split_ko(label))
    denom = max(1, len(genes))
    values = [len(genes & members) / denom for _, members in terms]
    values.append(float(len(genes)))
    return np.asarray(values, dtype=np.float32)


def load_papalexi_state() -> tuple[pd.DataFrame, list[str]]:
    adata = ad.read_h5ad("data/papalexi_small_pathway.h5ad")
    obs = adata.obs.copy()
    frame = pd.DataFrame(index=adata.obs_names)
    frame["ko_target"] = obs["ko_target"].astype(str).values
    for col in [c for c in obs.columns if c.startswith("pathway_")]:
        frame[col] = obs[col].astype(float).values
    if "protein" in adata.obsm:
        protein = np.asarray(adata.obsm["protein"])
        for i, name in enumerate(adata.uns["protein_names"]):
            frame[f"protein_{name}"] = protein[:, i]
    state_cols = [c for c in frame.columns if c.startswith("pathway_") or c.startswith("protein_")]
    return frame.dropna(subset=state_cols), state_cols


def fit_delta_baseline(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    holdout_kos: set[str],
) -> object:
    control_mean = frame.loc[control_mask(frame["ko_target"]), state_cols].mean().to_numpy(dtype=float)
    x_rows, y_rows = [], []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdout_kos:
            continue
        x_rows.append(ko_prior_vector(ko, terms))
        y_rows.append(group[state_cols].mean().to_numpy(dtype=float) - control_mean)
    x_train = np.vstack(x_rows)
    y_train = np.vstack(y_rows)
    n_components = min(6, x_train.shape[0] - 1, y_train.shape[1], x_train.shape[1])
    model = make_pipeline(StandardScaler(), PLSRegression(n_components=max(1, n_components), scale=True))
    model.fit(x_train, y_train)
    return model


class CVAE(nn.Module):
    def __init__(self, state_dim: int, cond_dim: int, latent_dim: int = 8, hidden: int = 64):
        super().__init__()
        self.encoder = nn.Sequential(nn.Linear(state_dim + cond_dim, hidden), nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU())
        self.mu = nn.Linear(hidden, latent_dim)
        self.logvar = nn.Linear(hidden, latent_dim)
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim + cond_dim, hidden),
            nn.SiLU(),
            nn.Linear(hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, state_dim),
        )

    def forward(self, x: torch.Tensor, c: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        h = self.encoder(torch.cat([x, c], dim=1))
        mu = self.mu(h)
        logvar = self.logvar(h).clamp(-6, 4)
        z = mu + torch.randn_like(mu) * torch.exp(0.5 * logvar)
        return self.decoder(torch.cat([z, c], dim=1)), mu, logvar

    def sample(self, c: torch.Tensor, n: int) -> torch.Tensor:
        z = torch.randn(n, self.mu.out_features)
        return self.decoder(torch.cat([z, c], dim=1))


class FlowNet(nn.Module):
    def __init__(self, state_dim: int, cond_dim: int, hidden: int = 96):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(state_dim + cond_dim + 1, hidden),
            nn.SiLU(),
            nn.Linear(hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, state_dim),
        )

    def forward(self, x: torch.Tensor, t: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        return self.net(torch.cat([x, t, c], dim=1))


class DiffusionNet(nn.Module):
    def __init__(self, state_dim: int, cond_dim: int, hidden: int = 96):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(state_dim + cond_dim + 1, hidden),
            nn.SiLU(),
            nn.Linear(hidden, hidden),
            nn.SiLU(),
            nn.Linear(hidden, state_dim),
        )

    def forward(self, x: torch.Tensor, t: torch.Tensor, c: torch.Tensor) -> torch.Tensor:
        return self.net(torch.cat([x, t, c], dim=1))


def make_training_arrays(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    holdout_kos: set[str],
) -> tuple[np.ndarray, np.ndarray]:
    rows_x, rows_c = [], []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdout_kos:
            continue
        states = group[state_cols].to_numpy(dtype=np.float32)
        cond = np.tile(ko_prior_vector(ko, terms), (len(states), 1))
        rows_x.append(states)
        rows_c.append(cond)
    return np.vstack(rows_x), np.vstack(rows_c)


def make_guided_training_arrays(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    holdout_kos: set[str],
    seed: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    rng = np.random.default_rng(seed)
    control = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=np.float32)
    control_mean = control.mean(axis=0)
    rows_target, rows_anchor, rows_residual, rows_c = [], [], [], []
    for ko, group in frame.loc[~control_mask(frame["ko_target"])].groupby("ko_target", observed=True):
        if ko in holdout_kos:
            continue
        target = group[state_cols].to_numpy(dtype=np.float32)
        delta = target.mean(axis=0) - control_mean
        ctrl = control[rng.integers(0, len(control), size=len(target))]
        anchor = ctrl + delta.reshape(1, -1)
        rows_target.append(target)
        rows_anchor.append(anchor)
        rows_residual.append(target - anchor)
        rows_c.append(np.tile(ko_prior_vector(ko, terms), (len(target), 1)))
    return np.vstack(rows_target), np.vstack(rows_anchor), np.vstack(rows_residual), np.vstack(rows_c)


def train_cvae(x: np.ndarray, c: np.ndarray, epochs: int, seed: int) -> CVAE:
    torch.manual_seed(seed)
    model = CVAE(x.shape[1], c.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
    ds = TensorDataset(torch.tensor(x), torch.tensor(c))
    loader = DataLoader(ds, batch_size=128, shuffle=True)
    for _ in range(epochs):
        for xb, cb in loader:
            recon, mu, logvar = model(xb, cb)
            recon_loss = ((recon - xb) ** 2).mean()
            kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
            loss = recon_loss + 0.02 * kl
            opt.zero_grad()
            loss.backward()
            opt.step()
    return model.eval()


def train_flow(
    x_ko: np.ndarray,
    c_ko: np.ndarray,
    control: np.ndarray,
    epochs: int,
    seed: int,
) -> FlowNet:
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = FlowNet(x_ko.shape[1], c_ko.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
    x_tensor = torch.tensor(x_ko)
    c_tensor = torch.tensor(c_ko)
    for _ in range(epochs):
        order = rng.permutation(len(x_ko))
        for start in range(0, len(order), 128):
            idx = order[start : start + 128]
            x1 = x_tensor[idx]
            cb = c_tensor[idx]
            x0 = torch.tensor(control[rng.integers(0, len(control), size=len(idx))], dtype=torch.float32)
            t = torch.rand(len(idx), 1)
            xt = (1 - t) * x0 + t * x1
            target_v = x1 - x0
            pred_v = model(xt, t, cb)
            loss = ((pred_v - target_v) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
    return model.eval()


def train_flow_pairs(
    x0_train: np.ndarray,
    x1_train: np.ndarray,
    c_train: np.ndarray,
    epochs: int,
    seed: int,
) -> FlowNet:
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = FlowNet(x1_train.shape[1], c_train.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
    x0_tensor = torch.tensor(x0_train)
    x1_tensor = torch.tensor(x1_train)
    c_tensor = torch.tensor(c_train)
    for _ in range(epochs):
        order = rng.permutation(len(x1_train))
        for start in range(0, len(order), 128):
            idx = order[start : start + 128]
            x0 = x0_tensor[idx]
            x1 = x1_tensor[idx]
            cb = c_tensor[idx]
            t = torch.rand(len(idx), 1)
            xt = (1 - t) * x0 + t * x1
            target_v = x1 - x0
            pred_v = model(xt, t, cb)
            loss = ((pred_v - target_v) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
    return model.eval()


def train_diffusion(x: np.ndarray, c: np.ndarray, epochs: int, seed: int) -> DiffusionNet:
    torch.manual_seed(seed)
    model = DiffusionNet(x.shape[1], c.shape[1])
    opt = torch.optim.AdamW(model.parameters(), lr=2e-3, weight_decay=1e-4)
    ds = TensorDataset(torch.tensor(x), torch.tensor(c))
    loader = DataLoader(ds, batch_size=128, shuffle=True)
    for _ in range(epochs):
        for xb, cb in loader:
            t = torch.rand(len(xb), 1)
            noise = torch.randn_like(xb)
            alpha = 1.0 - 0.85 * t
            xt = torch.sqrt(alpha) * xb + torch.sqrt(1 - alpha) * noise
            pred = model(xt, t, cb)
            loss = ((pred - noise) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
    return model.eval()


def flow_sample(model: FlowNet, x0: np.ndarray, cond: np.ndarray, steps: int = 20) -> np.ndarray:
    x = torch.tensor(x0, dtype=torch.float32)
    c = torch.tensor(np.tile(cond, (len(x0), 1)), dtype=torch.float32)
    with torch.no_grad():
        for i in range(steps):
            t_value = (i + 0.5) / steps
            t = torch.full((len(x0), 1), t_value)
            x = x + model(x, t, c) / steps
    return x.numpy()


def diffusion_sample(model: DiffusionNet, cond: np.ndarray, n: int, state_dim: int, steps: int = 40) -> np.ndarray:
    x = torch.randn(n, state_dim)
    c = torch.tensor(np.tile(cond, (n, 1)), dtype=torch.float32)
    with torch.no_grad():
        for i in reversed(range(steps)):
            t_value = (i + 1) / steps
            t = torch.full((n, 1), t_value)
            pred_noise = model(x, t, c)
            alpha = 1.0 - 0.85 * t
            x0_est = (x - torch.sqrt(1 - alpha) * pred_noise) / torch.sqrt(alpha)
            if i > 0:
                next_t = i / steps
                next_alpha = 1.0 - 0.85 * next_t
                x = torch.sqrt(torch.tensor(next_alpha)) * x0_est + torch.sqrt(torch.tensor(1 - next_alpha)) * pred_noise
            else:
                x = x0_est
    return x.numpy()


def evaluate_models(
    frame: pd.DataFrame,
    state_cols: list[str],
    terms: list[tuple[str, set[str]]],
    models: dict[str, object],
    scalers: tuple[StandardScaler, StandardScaler, StandardScaler],
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    rng = np.random.default_rng(seed)
    state_scaler, cond_scaler, residual_scaler = scalers
    control_raw = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=np.float32)
    control_scaled = state_scaler.transform(control_raw).astype(np.float32)
    metric_rows, cell_rows = [], []

    for offset, ko in enumerate(HOLDOUT_KOS):
        true_raw = frame.loc[frame["ko_target"].astype(str) == ko, state_cols].to_numpy(dtype=np.float32)
        if len(true_raw) == 0:
            continue
        ctrl_raw = control_raw[rng.integers(0, len(control_raw), size=len(true_raw))]
        ctrl_scaled = state_scaler.transform(ctrl_raw).astype(np.float32)
        cond_raw = ko_prior_vector(ko, terms).reshape(1, -1)
        cond_scaled = cond_scaler.transform(cond_raw).astype(np.float32).reshape(-1)

        generated = {"control": ctrl_raw}
        delta = models["Residual baseline"].predict(cond_raw).reshape(-1)
        generated["Residual baseline"] = ctrl_raw + delta.reshape(1, -1)

        c_tensor = torch.tensor(np.tile(cond_scaled, (len(true_raw), 1)), dtype=torch.float32)
        with torch.no_grad():
            cvae_scaled = models["Conditional VAE"].sample(c_tensor, len(true_raw)).numpy()
        generated["Conditional VAE"] = state_scaler.inverse_transform(cvae_scaled)
        generated["Flow matching"] = state_scaler.inverse_transform(
            flow_sample(models["Flow matching"], ctrl_scaled, cond_scaled)
        )
        generated["Diffusion"] = state_scaler.inverse_transform(
            diffusion_sample(models["Diffusion"], cond_scaled, len(true_raw), len(state_cols))
        )
        anchor_raw = generated["Residual baseline"]
        anchor_scaled = state_scaler.transform(anchor_raw).astype(np.float32)
        with torch.no_grad():
            guided_cvae_scaled = models["Guided Conditional VAE"].sample(c_tensor, len(true_raw)).numpy()
        generated["Guided Conditional VAE"] = anchor_raw + residual_scaler.inverse_transform(guided_cvae_scaled)
        generated["Guided Flow matching"] = state_scaler.inverse_transform(
            flow_sample(models["Guided Flow matching"], anchor_scaled, cond_scaled)
        )
        guided_diffusion_scaled = diffusion_sample(models["Guided Diffusion"], cond_scaled, len(true_raw), len(state_cols))
        generated["Guided Diffusion"] = anchor_raw + residual_scaler.inverse_transform(guided_diffusion_scaled)

        for model_name, matrix in generated.items():
            if model_name == "control":
                continue
            for j, feature in enumerate(state_cols):
                w_control = wasserstein_distance(true_raw[:, j], ctrl_raw[:, j])
                w_model = wasserstein_distance(true_raw[:, j], matrix[:, j])
                metric_rows.append(
                    {
                        "model": model_name,
                        "ko_target": ko,
                        "feature": feature,
                        "true_mean": true_raw[:, j].mean(),
                        "virtual_mean": matrix[:, j].mean(),
                        "control_mean": ctrl_raw[:, j].mean(),
                        "wasserstein_true_vs_virtual": w_model,
                        "wasserstein_true_vs_control": w_control,
                        "distribution_improvement": 1.0 - w_model / w_control if w_control > 1e-9 else np.nan,
                    }
                )

        for state, matrix in [("control cells", ctrl_raw), ("true KO cells", true_raw)] + [
            (f"{name} cells", mat) for name, mat in generated.items() if name != "control"
        ]:
            take = min(120, len(matrix))
            idx = rng.choice(len(matrix), size=take, replace=False)
            tmp = pd.DataFrame(matrix[idx], columns=state_cols)
            tmp["ko_target"] = ko
            tmp["state"] = state
            cell_rows.append(tmp)

    return pd.DataFrame(metric_rows), pd.concat(cell_rows, ignore_index=True)


def plot_model_summary(metrics: pd.DataFrame) -> None:
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("model", observed=True)
        .agg(
            mean_distribution_improvement=("distribution_improvement", "mean"),
            improved_fraction=("improved", "mean"),
        )
        .reset_index()
    )
    order = summary.sort_values("mean_distribution_improvement", ascending=False)["model"].tolist()
    plot = summary.melt(id_vars="model", var_name="metric", value_name="value")
    plot["model_label"] = plot["model"].map(MODEL_LABELS)
    order_labels = [MODEL_LABELS[model] for model in order]
    plot["metric"] = plot["metric"].map(
        {
            "mean_distribution_improvement": "Mean distribution improvement",
            "improved_fraction": "Fraction improved",
        }
    )
    g = sns.catplot(
        data=plot,
        x="model_label",
        y="value",
        col="metric",
        kind="bar",
        order=order_labels,
        sharey=False,
        height=4.2,
        aspect=1.15,
    )
    for ax in g.axes.flat:
        ax.axhline(0, color="0.25", linewidth=1)
        ax.tick_params(axis="x", rotation=22)
        ax.set_xlabel("")
    g.set_titles("{col_name}")
    g.fig.suptitle("Cell-level conditional generators on held-out single-gene KOs", y=1.04)
    g.savefig(FIG_DIR / "papalexi_cell_level_generator_deep_model_comparison.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def plot_model_ko_heatmap(metrics: pd.DataFrame) -> None:
    table = metrics.groupby(["model", "ko_target"], observed=True)["distribution_improvement"].mean().unstack("ko_target")
    order = metrics.groupby("model", observed=True)["distribution_improvement"].mean().sort_values(ascending=False).index
    table = table.loc[order]
    table.index = [MODEL_LABELS.get(model, model) for model in table.index]
    vmax = np.nanmax(np.abs(table.to_numpy()))
    plt.figure(figsize=(7.5, 4.7))
    ax = sns.heatmap(table, cmap="vlag", center=0, vmin=-vmax, vmax=vmax, annot=True, fmt=".2f", cbar_kws={"label": "mean distribution improvement"})
    ax.set_xlabel("held-out KO")
    ax.set_ylabel("")
    ax.set_title("Which generator works for which held-out KO?")
    plt.tight_layout()
    plt.savefig(FIG_DIR / "papalexi_cell_level_generator_model_by_ko_heatmap.png", bbox_inches="tight", dpi=300)
    plt.close()


def plot_key_distributions(cells: pd.DataFrame) -> None:
    selected = ["pathway_IFNG_JAK_STAT", "protein_PDL1", "protein_CD86"]
    selected = [feature for feature in selected if feature in cells.columns]
    keep_states = [
        "control cells",
        "true KO cells",
        "Residual baseline cells",
        "Conditional VAE cells",
        "Flow matching cells",
        "Diffusion cells",
        "Guided Conditional VAE cells",
        "Guided Flow matching cells",
        "Guided Diffusion cells",
    ]
    plot = cells.loc[cells["state"].isin(keep_states)].melt(
        id_vars=["ko_target", "state"],
        value_vars=selected,
        var_name="feature",
        value_name="score",
    )
    plot["state"] = plot["state"].map(STATE_LABELS)
    plot["panel"] = plot["ko_target"] + "\n" + plot["feature"].str.replace("pathway_", "", regex=False).str.replace("protein_", "", regex=False)
    g = sns.catplot(
        data=plot,
        x="state",
        y="score",
        col="panel",
        col_wrap=3,
        kind="box",
        showfliers=False,
        height=3.3,
        aspect=1.05,
        sharey=False,
        color="#8DA0CB",
    )
    for ax in g.axes.flat:
        ax.tick_params(axis="x", rotation=55, labelsize=8)
        ax.set_xlabel("")
        ax.set_ylabel("single-cell score")
    g.set_titles("{col_name}")
    g.fig.suptitle("Generated cell distributions from different conditional generators", y=1.02)
    g.savefig(FIG_DIR / "papalexi_cell_level_generator_deep_model_distributions.png", bbox_inches="tight", dpi=300)
    plt.close(g.fig)


def write_doc(metrics: pd.DataFrame) -> None:
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("model", observed=True)
        .agg(
            mean_distribution_improvement=("distribution_improvement", "mean"),
            improved_fraction=("improved", "mean"),
        )
        .reset_index()
        .sort_values("mean_distribution_improvement", ascending=False)
    )
    lines = [
        f"- {row.model}: 平均分布改进 {row.mean_distribution_improvement:.3f}, 改进比例 {row.improved_fraction:.1%}"
        for row in summary.itertuples(index=False)
    ]
    text = f"""# Cell-level 条件生成模型比较：Residual / VAE / Flow / Diffusion

这一步把三个更强的生成模型接到了同一个小样本多模态虚拟敲除框架里，并和稳定 residual baseline 比较。

## 比较对象

- Residual baseline：control cell + 系统先验预测的 KO 平均移动。
- Conditional VAE：学习 KO 条件下的低维潜在分布，再生成 KO-like cells。
- Flow matching：学习从 control cell 连续移动到 KO cell 的速度场。
- Diffusion：从噪声逐步去噪，生成 KO-like cells。

## 公平测试方式

完全留出这些 KO，不参与训练：

```text
STAT1, JAK2, IFNGR2, IRF1
```

所有模型都使用同一组 pathway/protein state 和同一套 Reactome/MSigDB/TF-target/PPI 条件先验。

## 当前结果

{chr(10).join(lines)}

解释：平均分布改进大于 0 表示生成细胞比 control 更接近真实 KO；改进比例表示多少个 KO-特征组合方向是有帮助的。

## 图

- `results/figures/papalexi_cell_level_generator_deep_model_comparison.png`
- `results/figures/papalexi_cell_level_generator_deep_model_distributions.png`

## 当前结论

小样本下，复杂模型不一定天然更好。这个实验的价值在于把 VAE / flow matching / diffusion 放到了同一评价框架里，能清楚看到哪种模型真的比稳定 baseline 更接近真实 KO 分布。

如果 flow matching 表现最好，下一步应把它作为主线，因为它最符合“敲除前 control cell -> 敲除后 KO cell”的状态移动假设。若 diffusion 或 VAE 表现不佳，说明当前样本量还不够支撑从噪声自由生成完整 KO 分布，需要加入 residual baseline 的方向约束或预训练。
"""
    Path("docs/cell_level_deep_generator_comparison.md").write_text(text, encoding="utf-8")


def write_doc(metrics: pd.DataFrame) -> None:
    summary = (
        metrics.assign(improved=metrics["distribution_improvement"] > 0)
        .groupby("model", observed=True)
        .agg(
            mean_distribution_improvement=("distribution_improvement", "mean"),
            improved_fraction=("improved", "mean"),
        )
        .reset_index()
        .sort_values("mean_distribution_improvement", ascending=False)
    )
    lines = [
        f"- {row.model}: 平均分布改进 {row.mean_distribution_improvement:.3f}, 改进比例 {row.improved_fraction:.1%}"
        for row in summary.itertuples(index=False)
    ]
    text = f"""# Cell-level 条件生成模型比较：Residual / VAE / Flow / Diffusion

这一步把三类生成模型放到同一个小样本多模态虚拟敲除框架里比较，并额外测试了“接在稳定 residual baseline 上”的 guided 版本。

## 比较对象

- Residual baseline：control cell + 系统先验预测的 KO 平均移动。
- Conditional VAE：直接学习 KO 条件下的潜在分布。
- Flow matching：直接学习从 control cell 到 KO cell 的移动速度场。
- Diffusion：直接从噪声去噪生成 KO-like cells。
- Guided VAE / Guided Flow / Guided Diffusion：先用 residual baseline 得到 KO anchor，再让深度模型学习剩余细胞级修正。

## 公平测试方式

完全留出这些 KO，不参与训练：

```text
STAT1, JAK2, IFNGR2, IRF1
```

所有模型都使用同一组 pathway/protein state 和同一套 Reactome/MSigDB/TF-target/PPI 条件先验。

## 当前结果

{chr(10).join(lines)}

解释：平均分布改进大于 0 表示生成细胞比 control 更接近真实 KO；改进比例表示多少个 KO-特征组合是有帮助的。

## 图

- `results/figures/papalexi_cell_level_generator_deep_model_comparison.png`
- `results/figures/papalexi_cell_level_generator_model_by_ko_heatmap.png`
- `results/figures/papalexi_cell_level_generator_deep_model_distributions.png`

## 当前结论

当前小样本条件下，复杂生成模型没有超过稳定 residual baseline。直接版 VAE、flow matching、diffusion 多数为负；guided 版本也没有改善，说明现在的细胞级残差还没有足够稳定的规律可供小模型学习。

这个结果支持一个重要开发判断：下一步不应该盲目加深模型，而应该先增强约束，包括更强的 KO 方向先验、更多训练 KO、跨数据预训练，或者把 residual baseline 作为 hard constraint，只允许生成模型学习低幅度不确定性。
"""
    Path("docs/cell_level_deep_generator_comparison.md").write_text(text, encoding="utf-8")


def main() -> None:
    setup_plot()
    seed = 53
    np.random.seed(seed)
    torch.manual_seed(seed)
    frame, state_cols = load_papalexi_state()
    perturb_genes = {gene for ko in frame["ko_target"].unique() for gene in split_ko(ko)}
    terms = select_prior_terms(Path("data/priors"), perturb_genes)

    x_raw, c_raw = make_training_arrays(frame, state_cols, terms, set(HOLDOUT_KOS))
    target_raw, anchor_raw, residual_raw, guided_c_raw = make_guided_training_arrays(
        frame, state_cols, terms, set(HOLDOUT_KOS), seed=seed
    )
    state_scaler = StandardScaler().fit(x_raw)
    cond_scaler = StandardScaler().fit(c_raw)
    residual_scaler = StandardScaler().fit(residual_raw)
    x_train = state_scaler.transform(x_raw).astype(np.float32)
    c_train = cond_scaler.transform(c_raw).astype(np.float32)
    target_train = state_scaler.transform(target_raw).astype(np.float32)
    anchor_train = state_scaler.transform(anchor_raw).astype(np.float32)
    residual_train = residual_scaler.transform(residual_raw).astype(np.float32)
    guided_c_train = cond_scaler.transform(guided_c_raw).astype(np.float32)
    control_raw = frame.loc[control_mask(frame["ko_target"]), state_cols].to_numpy(dtype=np.float32)
    control_scaled = state_scaler.transform(control_raw).astype(np.float32)

    print(f"Training cells={len(x_train)}, state_dim={x_train.shape[1]}, cond_dim={c_train.shape[1]}")
    residual = fit_delta_baseline(frame, state_cols, terms, set(HOLDOUT_KOS))
    cvae = train_cvae(x_train, c_train, epochs=450, seed=seed)
    flow = train_flow(x_train, c_train, control_scaled, epochs=500, seed=seed + 1)
    diffusion = train_diffusion(x_train, c_train, epochs=550, seed=seed + 2)
    guided_cvae = train_cvae(residual_train, guided_c_train, epochs=350, seed=seed + 3)
    guided_flow = train_flow_pairs(anchor_train, target_train, guided_c_train, epochs=400, seed=seed + 4)
    guided_diffusion = train_diffusion(residual_train, guided_c_train, epochs=450, seed=seed + 5)

    models = {
        "Residual baseline": residual,
        "Conditional VAE": cvae,
        "Flow matching": flow,
        "Diffusion": diffusion,
        "Guided Conditional VAE": guided_cvae,
        "Guided Flow matching": guided_flow,
        "Guided Diffusion": guided_diffusion,
    }
    metrics, cells = evaluate_models(frame, state_cols, terms, models, (state_scaler, cond_scaler, residual_scaler), seed=seed)
    metrics.to_csv("results/papalexi_cell_level_deep_generator_metrics.csv", index=False)
    cells.to_csv("results/papalexi_cell_level_deep_generator_cells.csv", index=False)

    plot_model_summary(metrics)
    plot_model_ko_heatmap(metrics)
    plot_key_distributions(cells)
    write_doc(metrics)

    summary = metrics.groupby("model", observed=True)["distribution_improvement"].mean().sort_values(ascending=False)
    print(summary.round(3).to_string())
    print("Saved deep generator comparison metrics, cells, figures, and docs.")


if __name__ == "__main__":
    main()
