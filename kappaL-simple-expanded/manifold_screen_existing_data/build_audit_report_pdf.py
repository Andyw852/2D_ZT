"""Build a Chinese PDF report for the structure-electronic manifold audit."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
from PIL import Image as PILImage
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4, landscape
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = ROOT / "output" / "pdf"
OUTPUT_PDF = OUTPUT_DIR / "structure_electronic_manifold_audit_report_cn.pdf"

FIG_OLD = ROOT / "figures" / "joint_structure_electronic_manifold_screen.png"
FIG_STRICT = ROOT / "figures" / "strict_and_joint_manifold.png"
FIG_AUDIT = ROOT / "figures" / "consensus_structure_electronic_audit.png"
SUMMARY_PATH = ROOT / "outputs" / "consensus_audit_summary.json"
CANDIDATES_PATH = ROOT / "outputs" / "consensus_candidates.csv"

FONT_REGULAR = "C:/Windows/Fonts/msyh.ttc"
FONT_BOLD = "C:/Windows/Fonts/msyhbd.ttc"
pdfmetrics.registerFont(TTFont("CN", FONT_REGULAR, subfontIndex=0))
pdfmetrics.registerFont(TTFont("CN-Bold", FONT_BOLD, subfontIndex=0))

PAGE = landscape(A4)
NAVY = colors.HexColor("#18324A")
BLUE = colors.HexColor("#2F6FED")
ORANGE = colors.HexColor("#F28E2B")
PURPLE = colors.HexColor("#7D3C98")
CYAN = colors.HexColor("#00A9B7")
PALE = colors.HexColor("#F4F7FA")
MID = colors.HexColor("#D9E1E8")
TEXT = colors.HexColor("#26333D")
MUTED = colors.HexColor("#65727E")
RED = colors.HexColor("#B23A48")


styles = getSampleStyleSheet()
TITLE = ParagraphStyle(
    "CNTitle", fontName="CN-Bold", fontSize=24, leading=31,
    textColor=NAVY, alignment=TA_LEFT, spaceAfter=7 * mm,
)
SUBTITLE = ParagraphStyle(
    "CNSubtitle", fontName="CN", fontSize=11.5, leading=18,
    textColor=MUTED, alignment=TA_LEFT,
)
H1 = ParagraphStyle(
    "CNH1", fontName="CN-Bold", fontSize=17, leading=22,
    textColor=NAVY, spaceAfter=4 * mm,
)
H2 = ParagraphStyle(
    "CNH2", fontName="CN-Bold", fontSize=11.5, leading=16,
    textColor=NAVY, spaceAfter=1.8 * mm,
)
BODY = ParagraphStyle(
    "CNBody", fontName="CN", fontSize=9.2, leading=14.2,
    textColor=TEXT, alignment=TA_LEFT, spaceAfter=2.2 * mm,
)
SMALL = ParagraphStyle(
    "CNSmall", fontName="CN", fontSize=7.5, leading=10.5,
    textColor=TEXT,
)
CAPTION = ParagraphStyle(
    "CNCaption", fontName="CN", fontSize=7.6, leading=11,
    textColor=MUTED, spaceBefore=1.5 * mm,
)
METRIC = ParagraphStyle(
    "CNMetric", fontName="CN-Bold", fontSize=19, leading=22,
    textColor=NAVY, alignment=TA_CENTER,
)
METRIC_LABEL = ParagraphStyle(
    "CNMetricLabel", fontName="CN", fontSize=8, leading=11,
    textColor=MUTED, alignment=TA_CENTER,
)
CALLOUT = ParagraphStyle(
    "CNCallout", fontName="CN-Bold", fontSize=11, leading=17,
    textColor=PURPLE,
)


def p(text: str, style=BODY) -> Paragraph:
    return Paragraph(text, style)


def bullets(items: list[str], style=BODY) -> list[Paragraph]:
    return [Paragraph(f"• {item}", style) for item in items]


def header_footer(canvas, doc) -> None:
    canvas.saveState()
    width, height = PAGE
    canvas.setStrokeColor(MID)
    canvas.setLineWidth(0.5)
    canvas.line(12 * mm, height - 10 * mm, width - 12 * mm, height - 10 * mm)
    canvas.setFont("CN", 7.5)
    canvas.setFillColor(MUTED)
    canvas.drawString(13 * mm, height - 7.2 * mm, "结构-电子流形热电筛选审计")
    canvas.drawRightString(width - 13 * mm, 7 * mm, f"第 {doc.page} 页")
    canvas.restoreState()


def metric_card(value: str, label: str, accent) -> Table:
    data = [[p(value, METRIC)], [p(label, METRIC_LABEL)]]
    table = Table(data, colWidths=[56 * mm], rowHeights=[13 * mm, 11 * mm])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), PALE),
                ("BOX", (0, 0), (-1, -1), 0.8, accent),
                ("LINEABOVE", (0, 0), (-1, 0), 3.0, accent),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ]
        )
    )
    return table


def image_flowable(path: Path, max_width: float, max_height: float) -> Image:
    with PILImage.open(path) as img:
        width, height = img.size
    scale = min(max_width / width, max_height / height)
    return Image(str(path), width=width * scale, height=height * scale)


def note_panel(title: str, sections: list[tuple[str, list[str]]]) -> Table:
    flow = [p(title, H2)]
    for heading, items in sections:
        flow.append(p(heading, ParagraphStyle(
            f"note-{heading}", parent=H2, fontSize=9.2, leading=13,
            textColor=PURPLE, spaceBefore=1.7 * mm, spaceAfter=1 * mm,
        )))
        flow.extend(bullets(items, SMALL))
        flow.append(Spacer(1, 1.2 * mm))
    table = Table([[flow]], colWidths=[70 * mm])
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), PALE),
                ("BOX", (0, 0), (-1, -1), 0.6, MID),
                ("LINEABOVE", (0, 0), (-1, 0), 3.0, PURPLE),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 8),
                ("RIGHTPADDING", (0, 0), (-1, -1), 8),
                ("TOPPADDING", (0, 0), (-1, -1), 8),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
            ]
        )
    )
    return table


def figure_page(
    title: str,
    figure: Path,
    sections: list[tuple[str, list[str]]],
    caption: str,
) -> list:
    img = image_flowable(figure, 190 * mm, 145 * mm)
    notes = note_panel("与图对应的说明", sections)
    content = Table(
        [[img, notes]], colWidths=[194 * mm, 72 * mm], hAlign="LEFT"
    )
    content.setStyle(
        TableStyle(
            [
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 0),
                ("RIGHTPADDING", (0, 0), (0, 0), 4 * mm),
                ("RIGHTPADDING", (1, 0), (1, 0), 0),
                ("TOPPADDING", (0, 0), (-1, -1), 0),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 0),
            ]
        )
    )
    return [p(title, H1), content, p(caption, CAPTION), PageBreak()]


def comparison_table() -> Table:
    rows = [
        ["层次", "输入/定义", "可以回答", "不能回答"],
        [
            "结构空间 S1",
            "组成统计、密度、原子体积、晶胞大小与几何畸变",
            "结构与组成是否相似",
            "不等于完整声子空间；没有声子谱、非谐性与散射率",
        ],
        [
            "电子空间 E2",
            "带隙、介电张量、电子/空穴有效质量及各向异性代理",
            "基础电子结构是否相似",
            "不等于完整输运空间；没有弛豫时间、迁移率和形变势",
        ],
        [
            "双通道得分",
            "实验 kL 近 300 K + JARVIS PF 600 K、固定载流子浓度",
            "在当前代理标签下是否同时有利",
            "不能直接代表任意温度、任意掺杂下的 zT",
        ],
        [
            "高 zT 星标",
            "StarryData2 最大 zT 按约化化学式连接 JARVIS",
            "提供家族级弱种子",
            "不能确定具体晶相、掺杂、维度、微结构和测量条件",
        ],
    ]
    data = [[p(str(cell), SMALL) for cell in row] for row in rows]
    table = Table(data, colWidths=[31 * mm, 77 * mm, 63 * mm, 91 * mm], repeatRows=1)
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), NAVY),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "CN-Bold"),
                ("GRID", (0, 0), (-1, -1), 0.45, MID),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("LEFTPADDING", (0, 0), (-1, -1), 5),
                ("RIGHTPADDING", (0, 0), (-1, -1), 5),
                ("TOPPADDING", (0, 0), (-1, -1), 6),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
            ]
        )
    )
    return table


def candidate_table(candidates: pd.DataFrame) -> Table:
    rows = [["材料", "同种子", "dual", "结构", "电子"]]
    for row in candidates.sort_values(
        ["same_seed_and_similarity", "dual_score"], ascending=False
    ).itertuples(index=False):
        rows.append(
            [
                str(row.formula),
                str(row.same_seed_formula),
                f"{row.dual_score:.3f}",
                f"{row.same_seed_structure_similarity:.3f}",
                f"{row.same_seed_electronic_similarity:.3f}",
            ]
        )
    data = [[p(cell, SMALL) for cell in row] for row in rows]
    table = Table(
        data,
        colWidths=[25 * mm, 31 * mm, 21 * mm, 23 * mm, 23 * mm],
        repeatRows=1,
    )
    table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), PURPLE),
                ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
                ("FONTNAME", (0, 0), (-1, 0), "CN-Bold"),
                ("GRID", (0, 0), (-1, -1), 0.4, MID),
                ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, PALE]),
                ("ALIGN", (2, 1), (-1, -1), "CENTER"),
                ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
                ("LEFTPADDING", (0, 0), (-1, -1), 4),
                ("RIGHTPADDING", (0, 0), (-1, -1), 4),
                ("TOPPADDING", (0, 0), (-1, -1), 3.5),
                ("BOTTOMPADDING", (0, 0), (-1, -1), 3.5),
            ]
        )
    )
    return table


def build() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    summary = json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))
    candidates = pd.read_csv(CANDIDATES_PATH)

    doc = SimpleDocTemplate(
        str(OUTPUT_PDF), pagesize=PAGE,
        leftMargin=13 * mm, rightMargin=13 * mm,
        topMargin=14 * mm, bottomMargin=12 * mm,
        title="结构-电子流形热电筛选审计报告",
        author="Codex / existing-data audit",
        subject="结构空间、电子空间、strict-AND流形与热电候选筛选",
    )
    story = []

    # Cover and executive summary
    story.append(Spacer(1, 10 * mm))
    story.append(p("结构-电子流形热电筛选审计报告", TITLE))
    story.append(p(
        "对旧联合流形、无权重 strict-AND 流形和同种子共识筛选的统一解释",
        SUBTITLE,
    ))
    story.append(Spacer(1, 8 * mm))
    cards = Table(
        [[
            metric_card(f"{summary['n_materials']:,}", "完整描述符材料数", BLUE),
            metric_card(f"{summary['n_seed_formulas']} / {summary['n_seed_rows']}", "高 zT 公式 / JARVIS结构", CYAN),
            metric_card(f"{summary['old_independent_overlap']} / 30", "旧图与独立Top-30重合", ORANGE),
            metric_card(str(summary["n_consensus_candidates"]), "严格共识候选", PURPLE),
        ]],
        colWidths=[63 * mm] * 4,
    )
    cards.setStyle(TableStyle([("VALIGN", (0, 0), (-1, -1), "TOP")]))
    story.append(cards)
    story.append(Spacer(1, 8 * mm))
    story.append(p("执行摘要", H1))
    story.append(p(
        "严格 AND 构图本身没有出现断图或明显数值失败。紫色候选与青色星标在新版二维图中分开，主要暴露了三类问题：旧图的候选得分预先包含流形邻近性；UMAP 不能可靠保持所有全维近邻；高 zT 标签仅按约化化学式映射，不能确定具体晶相和实验条件。",
        BODY,
    ))
    story.append(p(
        "因此，旧图中“紫色贴近青色”的视觉效果不能作为独立证据；新版分离也不能单独证明 strict-AND 失败。当前最可靠的做法，是在全维空间中直接要求候选相对于同一个高 zT 种子同时满足结构近邻和电子近邻，再把二维图仅作为沟通工具。",
        CALLOUT,
    ))
    story.append(Spacer(1, 4 * mm))
    story.append(p(
        "范围：仅使用现有本地数据；没有新增 DFT、BTE、声子或输运计算。报告日期：2026-08-30。",
        CAPTION,
    ))
    story.append(PageBreak())

    # Definitions and audit boundary
    story.append(p("1. 四个层次必须分开解释", H1))
    story.append(p(
        "当前工作同时出现“描述符空间、监督得分、高 zT 种子、二维投影”。它们承担的证据角色不同。如果混在同一张图里，很容易把筛选规则误读成验证结果。",
        BODY,
    ))
    story.append(comparison_table())
    story.append(Spacer(1, 6 * mm))
    mismatch = Table(
        [[
            [p("标签条件不一致", H2)] + bullets([
                "晶格通道：实验 kL 近 300 K。",
                "电子通道：JARVIS PF 在 600 K、固定载流子浓度、CRTA。",
                "星标：200-1500 K 范围内最大实验/参考 zT，掺杂与微结构未统一。",
            ], SMALL),
            [p("公式映射不等于相匹配", H2)] + bullets([
                "10 个高 zT 约化公式被投到 33 个 JARVIS 结构。",
                "C 对应 g=0.45；SnS2 对应 3L；Si 对应 nano-bulk Si(model)。",
                "留一公式检索取同公式多个结构的最高分，会产生多次尝试的乐观偏差。",
            ], SMALL),
            [p("正确证据链", H2)] + bullets([
                "候选筛选：可使用全维结构+电子同种子 AND。",
                "方法验证：必须相分辨、条件一致，并在冻结阈值后测试。",
                "UMAP：仅显示，不负责定义候选或证明性能。",
            ], SMALL),
        ]],
        colWidths=[85 * mm, 85 * mm, 85 * mm],
    )
    mismatch.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), PALE),
        ("BOX", (0, 0), (-1, -1), 0.6, MID),
        ("INNERGRID", (0, 0), (-1, -1), 0.4, MID),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
    ]))
    story.append(mismatch)
    story.append(PageBreak())

    story.extend(figure_page(
        "2. 旧联合流形：图形好看，但紫色定义包含邻近性",
        FIG_OLD,
        [
            ("图中分别是什么", [
                "左上：结构 S1；右上：电子 E2；左下：等权联合扩散流形；右下：留一公式检索。",
                "青色星为最大 zT≥1 的约化公式匹配；紫色为旧版 top-30 候选。",
            ]),
            ("为什么会显得聚集", [
                "旧紫色得分 = 流形种子邻近百分位与双通道得分的几何平均。",
                "因此“靠近星标”已经进入候选定义，不能再把图上接近当作独立验证。",
            ]),
            ("量化结果", [
                "结构 S1 留一公式中位检索约 0.95；联合 S1+E2 约 0.90。",
                "联合并未明显优于结构单视图，不能据此断言电子视图带来新增信息。",
            ]),
        ],
        "图 1。旧联合流形原图。坐标不使用 PF、kL 或 zT 构图，但紫色候选排序使用了到高 zT 种子的联合流形距离。",
    ))

    story.extend(figure_page(
        "3. 无权重 strict-AND：分离揭示独立分数与种子家族不一致",
        FIG_STRICT,
        [
            ("严格规则", [
                "每对材料取结构与电子近邻秩中较差者：r_AND=max(r_S,r_E)。",
                "一个视图非常相似不能补偿另一个视图非常不相似。",
            ]),
            ("紫色为什么不再贴星标", [
                "此图紫色仅由旧有 dual_score 选择，不使用流形距离。",
                "旧图紫色与独立紫色 top-30 只重合 1 个；二者不是同一批候选。",
            ]),
            ("构图是否失败", [
                "k=15、30、50 均为一个连通分量，没有出现图断裂。",
                "k=30 留一公式中位约 0.98、top-10% 为 7/10；但公式映射和多结构取最大值使该数值偏乐观。",
            ]),
            ("科学含义", [
                "分离不是坏图，而是说明简单 transport proxy 与当前高 zT 公式种子并不一致。",
            ]),
        ],
        "图 2。无权重 strict-AND 审计图。右上二维 UMAP 仅为全维扩散坐标的可视化，不用于候选评分。",
    ))

    story.extend(figure_page(
        "4. 同种子共识审计：把筛选、相似性和可视化分开",
        FIG_AUDIT,
        [
            ("左上：可解释的主判据", [
                "对每个材料选择使 max(r_S,r_E) 最小的同一个高 zT 种子。",
                "横轴和纵轴分别是相对于该同一种子的结构与电子相似度。",
                "右上区域要求两个视图都进入前 10%，不存在跨视图补偿。",
            ]),
            ("右上：UMAP为什么仍可能分开", [
                "紫线连接全维判据选出的对应种子；长线显示二维 UMAP 扭曲了部分全维近邻。",
                "因此以后应报告全维秩和邻域富集，不能用二维距离定量。",
            ]),
            ("左下：输运得分", [
                "紫色还必须处于现有 dual_score 前 10%。",
                "青色星在该平面分布很散，直接说明固定条件 PF、300 K kL 与最大实验 zT 不完全对齐。",
            ]),
            ("输出", [
                f"严格共识得到 {summary['n_consensus_candidates']} 个候选；独立 top-30 中仅 {summary['independent_top30_in_direct_and_top10']}/30 同时进入同种子 AND 前 10%。",
            ]),
        ],
        "图 3。同种子共识审计图。空心青色星强调其仅为约化公式匹配；紫色是筛选规则的结果，而非实验确认。",
    ))

    # Candidate table
    story.append(p("5. 严格共识候选清单", H1))
    story.append(p(
        "下表仅保留同时满足“dual_score 前 10%”以及“相对于同一个种子，结构和电子相似度均≥0.90”的材料。它们是优先核查对象，不是高 zT 预测值。",
        BODY,
    ))
    left = candidate_table(candidates.iloc[:11])
    right = candidate_table(candidates.iloc[11:])
    paired = Table([[left, right]], colWidths=[126 * mm, 126 * mm])
    paired.setStyle(TableStyle([
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (0, 0), 0),
        ("RIGHTPADDING", (0, 0), (0, 0), 4 * mm),
        ("LEFTPADDING", (1, 0), (1, 0), 4 * mm),
        ("RIGHTPADDING", (1, 0), (1, 0), 0),
    ]))
    story.append(paired)
    story.append(Spacer(1, 5 * mm))
    candidate_notes = Table(
        [[
            [p("相对更值得优先核查", H2)] + bullets([
                "HgTe、Sb2Se3、TlBiSe2、Bi4Te7Pb、GeSb4Te7 与已知窄带隙重元素热电家族更接近。",
                "这些名称仍需核对晶相、稳定性、载流子类型和实际掺杂窗口。",
            ], SMALL),
            [p("明显的精确率警告", H2)] + bullets([
                "PtI3、Cu2HgI4、TeI2、HgI2 等卤化物进入清单，可能由低 kL、重元素和简单有效质量代理共同推高。",
                "这表明当前描述符尚不能可靠排除分子晶体、绝缘体或难掺杂材料。",
            ], SMALL),
        ]],
        colWidths=[126 * mm, 126 * mm],
    )
    candidate_notes.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, -1), PALE),
        ("BOX", (0, 0), (-1, -1), 0.6, MID),
        ("INNERGRID", (0, 0), (-1, -1), 0.4, MID),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 8),
        ("RIGHTPADDING", (0, 0), (-1, -1), 8),
        ("TOPPADDING", (0, 0), (-1, -1), 8),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 8),
    ]))
    story.append(candidate_notes)
    story.append(PageBreak())

    # Final conclusion and next steps
    story.append(p("6. 正确结论与下一步", H1))
    conclusion = Table(
        [[
            [p("当前可以支持", H2)] + bullets([
                "结构与电子视图可以用非补偿的同种子近邻规则联合。",
                "现有数据中存在 21 个同时满足 transport proxy 和双视图邻近的类比候选。",
                "高 zT 材料并不必然形成单一全局二维簇，更合理的是多个家族局部邻域。",
                "二维 UMAP 适合展示，核心结论必须来自全维秩、富集和留出测试。",
            ]),
            [p("当前不能支持", H2)] + bullets([
                "不能把紫色点称为已预测高 zT 材料。",
                "不能把约化公式星标视为具体 JARVIS 晶相的实验标签。",
                "不能用旧图中紫色贴近星标证明方法有效，因为邻近性进入了候选排序。",
                "不能把 S1 和 E2 分别称为完整声子空间与完整输运空间。",
            ]),
        ]],
        colWidths=[128 * mm, 128 * mm],
    )
    conclusion.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (0, 0), colors.HexColor("#EEF7F4")),
        ("BACKGROUND", (1, 0), (1, 0), colors.HexColor("#FFF3F3")),
        ("BOX", (0, 0), (-1, -1), 0.7, MID),
        ("INNERGRID", (0, 0), (-1, -1), 0.5, MID),
        ("VALIGN", (0, 0), (-1, -1), "TOP"),
        ("LEFTPADDING", (0, 0), (-1, -1), 10),
        ("RIGHTPADDING", (0, 0), (-1, -1), 10),
        ("TOPPADDING", (0, 0), (-1, -1), 10),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 10),
    ]))
    story.append(conclusion)
    story.append(Spacer(1, 6 * mm))
    story.append(p("不增加第一性原理计算时，建议按以下顺序继续", H1))
    steps = [
        ("01", "先清理种子", "从 StarryData2 样品元数据中区分 Experiment/Reference、bulk/film/2D/model，并去除无法映射到具体相的星标。"),
        ("02", "统一条件", "至少按 300/600/900 K 分层；n 型与 p 型分开；不再用跨温度最大 zT 作为单一标签。"),
        ("03", "冻结主判据", "使用同一种子的 max(r_S,r_E) 作为主相似度；UMAP 不参与候选选择。"),
        ("04", "建立正负对照", "除高 zT 家族外加入低 zT、低 PF 的反例，报告精确率、召回率和候选富集，而不是只看星标是否聚集。"),
        ("05", "再决定是否扩展描述符", "若相分辨验证仍显示电子视图增益，再引入 DOS/谷简并/迁移率代理和更完整的结构键合描述。"),
    ]
    step_rows = []
    for number, heading, text in steps:
        step_rows.append([
            p(number, ParagraphStyle("stepn", parent=METRIC, fontSize=15, textColor=PURPLE)),
            p(heading, H2),
            p(text, BODY),
        ])
    step_table = Table(step_rows, colWidths=[18 * mm, 40 * mm, 195 * mm])
    step_table.setStyle(TableStyle([
        ("ROWBACKGROUNDS", (0, 0), (-1, -1), [PALE, colors.white]),
        ("BOX", (0, 0), (-1, -1), 0.5, MID),
        ("INNERGRID", (0, 0), (-1, -1), 0.35, MID),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 7),
        ("RIGHTPADDING", (0, 0), (-1, -1), 7),
        ("TOPPADDING", (0, 0), (-1, -1), 6),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 6),
    ]))
    story.append(step_table)
    story.append(Spacer(1, 6 * mm))
    story.append(p(
        "最终判断：当前项目最有价值的结果不是“找到一个统一高 zT 簇”，而是建立了一个可审计的双视图局部类比框架，并明确发现公式级标签和二维投影会制造过强结论。",
        CALLOUT,
    ))
    story.append(p("代码与输出测试：manifold_screen_existing_data/tests，6 passed。", CAPTION))

    doc.build(story, onFirstPage=header_footer, onLaterPages=header_footer)
    print(OUTPUT_PDF)


if __name__ == "__main__":
    build()
