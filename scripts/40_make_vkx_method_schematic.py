from pathlib import Path
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = ROOT / "docs" / "assets" / "draft_figures"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def font(size: int, bold: bool = False):
    candidates = [
        "C:/Windows/Fonts/msyhbd.ttc" if bold else "C:/Windows/Fonts/msyh.ttc",
        "C:/Windows/Fonts/simhei.ttf",
        "C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf",
    ]
    for path in candidates:
        if Path(path).exists():
            return ImageFont.truetype(path, size=size)
    return ImageFont.load_default()


F_TITLE = font(58, True)
F_H = font(34, True)
F_M = font(25, False)
F_S = font(21, False)
F_TINY = font(18, False)


def text_size(draw, text, fnt):
    box = draw.textbbox((0, 0), text, font=fnt)
    return box[2] - box[0], box[3] - box[1]


def wrap(draw, text, fnt, max_width):
    lines = []
    for raw in text.split("\n"):
        line = ""
        for ch in raw:
            test = line + ch
            if text_size(draw, test, fnt)[0] <= max_width or not line:
                line = test
            else:
                lines.append(line)
                line = ch
        if line:
            lines.append(line)
    return lines


def rounded_box(draw, xy, fill, outline, title, body, title_color="#103047"):
    x0, y0, x1, y1 = xy
    draw.rounded_rectangle(xy, radius=28, fill=fill, outline=outline, width=3)
    draw.text((x0 + 28, y0 + 22), title, font=F_H, fill=title_color)
    y = y0 + 78
    for line in wrap(draw, body, F_M, x1 - x0 - 56):
        draw.text((x0 + 28, y), line, font=F_M, fill="#263238")
        y += 34


def arrow(draw, start, end, color="#455A64", width=5):
    draw.line([start, end], fill=color, width=width)
    x0, y0 = start
    x1, y1 = end
    if abs(x1 - x0) >= abs(y1 - y0):
        sign = 1 if x1 > x0 else -1
        pts = [(x1, y1), (x1 - sign * 22, y1 - 12), (x1 - sign * 22, y1 + 12)]
    else:
        sign = 1 if y1 > y0 else -1
        pts = [(x1, y1), (x1 - 12, y1 - sign * 22), (x1 + 12, y1 - sign * 22)]
    draw.polygon(pts, fill=color)


def pill(draw, xy, text, fill, outline="#B0BEC5"):
    x0, y0, x1, y1 = xy
    draw.rounded_rectangle(xy, radius=18, fill=fill, outline=outline, width=2)
    w, h = text_size(draw, text, F_S)
    draw.text((x0 + (x1 - x0 - w) / 2, y0 + (y1 - y0 - h) / 2 - 2), text, font=F_S, fill="#1C2833")


def main():
    img = Image.new("RGB", (2600, 1700), "#FAFBFC")
    draw = ImageDraw.Draw(img)

    draw.text((90, 48), "VKX 小样本多模态虚拟敲除模型：算法流程图", font=F_TITLE, fill="#102A43")
    draw.text(
        (92, 125),
        "核心思想：先把单细胞/多组学矩阵转成可解释状态，再用系统先验和少量真实 perturbation 学 KO 方向，最后在 hard constraint 附近生成虚拟细胞。",
        font=F_M,
        fill="#52616B",
    )

    rounded_box(
        draw,
        (90, 220, 590, 560),
        "#E3F2FD",
        "#90CAF9",
        "1. 用户输入",
        "RNA 单细胞矩阵\n可选：ADT 蛋白\n可选：ATAC peak / gene activity\n可选：chromVAR motif\n有 KO 标签：训练/评估\n无 KO 标签：reference application",
    )
    rounded_box(
        draw,
        (720, 220, 1260, 560),
        "#E8F5E9",
        "#A5D6A7",
        "2. 可解释状态表示",
        "RNA → pathway / program score\n蛋白 → marker protein state\nATAC → gene activity\nATAC → motif + peak feature\n得到每个细胞的状态向量 S",
    )
    rounded_box(
        draw,
        (1390, 220, 1930, 560),
        "#FFF8E1",
        "#FFE082",
        "3. 系统先验",
        "Reactome / MSigDB pathway\nTF-target / PPI network\nmotif-to-peak annotation\npeak-gene linkage 和 locus 权重",
    )
    rounded_box(
        draw,
        (2060, 220, 2510, 560),
        "#F3E5F5",
        "#CE93D8",
        "4. KO 标签和真实效应",
        "control cells 与 KO cells 对比\n学习 ΔKO = KO均值 - control均值\n支持单敲、双敲、批量 KO",
    )

    arrow(draw, (590, 390), (720, 390))
    arrow(draw, (1260, 390), (1390, 390))
    arrow(draw, (1930, 390), (2060, 390))

    rounded_box(
        draw,
        (170, 720, 2430, 1090),
        "#FFFFFF",
        "#78909C",
        "5. VKX 模型核心：hard-constrained residual / PLS baseline",
        "A. 先验到 KO 方向：Δhat_z = fθ(q_z)，用 PLS / ridge / residual anchor 学少样本稳定方向\n"
        "B. 自适应选择：在训练 KO 上选择最稳 anchor，并做 feature-scale 幅度校准\n"
        "C. 双敲交互：Δhat_{a+b} = Δhat_a + Δhat_b + rhat_{a,b}，rhat 来自 interaction residual\n"
        "D. 多模态融合：RNA pathway + protein + ATAC motif/peak 共同约束同一个 KO 状态变化\n"
        "E. 生成虚拟细胞：S_virtual = S_control + Δhat_z + ε，其中 ε 只在 hard constraint 附近表示不确定性",
    )

    for i, label in enumerate(["PLS/Ridge", "Adaptive anchor", "Response boost", "Feature calibration", "Interaction residual", "Quantile ATAC shape"]):
        pill(draw, (270 + i * 360, 1125, 560 + i * 360, 1180), label, "#ECEFF1")

    arrow(draw, (1290, 560), (1290, 720))
    arrow(draw, (1290, 1180), (1290, 1300))

    rounded_box(
        draw,
        (90, 1300, 740, 1580),
        "#E0F7FA",
        "#80DEEA",
        "6. 输出结果",
        "KO delta 表格、虚拟 KO 细胞\npathway/protein/ATAC/peak 变化\nprediction-only 或 labelled benchmark 报告",
    )
    rounded_box(
        draw,
        (980, 1300, 1620, 1580),
        "#F1F8E9",
        "#C5E1A5",
        "7. 可视化",
        "ROC/AUC 曲线\nreal vs virtual KO heatmap\nbefore/after UMAP\nsingle vs double KO response map\npeak locus track 与 radar leaderboard",
    )
    rounded_box(
        draw,
        (1860, 1300, 2510, 1580),
        "#FFF3E0",
        "#FFCC80",
        "8. 评估指标",
        "方向一致性、MAE、R²、AUC\n细胞状态移动是否合理\nATAC 开放比例/分位数形状是否贴近真实 KO",
    )

    arrow(draw, (740, 1440), (980, 1440))
    arrow(draw, (1620, 1440), (1860, 1440))

    draw.text(
        (92, 1620),
        "注意：普通 10X 无 KO 标签时可做虚拟应用和状态变化展示，但不能在该数据内部报告真实准确率；有 perturbation 标签的数据才能做 labelled benchmark。",
        font=F_TINY,
        fill="#607D8B",
    )

    out = OUT_DIR / "01_vkx_model_algorithm_schematic.png"
    img.save(out, quality=95)
    print(out)


if __name__ == "__main__":
    main()
