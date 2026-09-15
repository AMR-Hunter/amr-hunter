from __future__ import annotations
import csv
from collections import Counter, defaultdict
from pathlib import Path
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, PathPatch
from matplotlib.path import Path as MplPath

BASE = Path("c:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\内部数据")
PER_ISOLATE = BASE / "amr_hunter_local_analysis" / "per_isolate_variants.csv"
SERVER = BASE / "final_calls_server.csv"
JUDGE = BASE / "bdq_variants_to_judge.csv"
OUT_DIR = Path(__file__).resolve().parents[1] / "figures" / "local_cohort"
OUT_DIR.mkdir(parents=True, exist_ok=True)
STRUCTURAL_GENES = {"mtrA", "mtrB", "mmpS5", "mmpL5", "lpqB"}
EXPLORATORY_GENES = {"glpK", "Rv1979c"}
TARGET_GENES = {
    "atpE",
    "Rv0678",
    "mmpS5",
    "mmpL5",
    "pepQ",
    "lpqB",
    "mtrA",
    "mtrB",
    "Rv1979c",
    "glpK",
}
LABEL_MAP = {"RESISTANT": "R", "SENSITIVE": "S", "INDETERMINATE": "INDET"}
C_R = "#d62728"
C_S = "#2ca02c"
C_I = "#9aa5b1"
C_JUDGED = "#4c72b0"
C_NJUDGE = "#c2c8d0"
C_NOTUP = "#8a919c"
plt.rcParams.update(
    {"font.size": 10, "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight"}
)


def norm(call: str) -> str:
    return LABEL_MAP.get(call.strip().upper(), call.strip().upper())


def re_resolve(gene: str, call: str) -> str:
    call = norm(call)
    if gene in STRUCTURAL_GENES and call == "R":
        return "S"
    if gene in EXPLORATORY_GENES and call == "R":
        return "INDET"
    return call


def load_flows() -> dict:
    server_rows = {}
    with SERVER.open(encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            server_rows[row["gene"], row["pos"], row["ref"], row["alt"]] = row
    in_server = set(server_rows)
    judge_keys = set()
    with JUDGE.open(encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            judge_keys.add((row["gene"], row["pos"], row["ref"], row["alt"]))
    judge_keys = {k for k in judge_keys if k[0] in TARGET_GENES}
    key_info: dict[tuple, dict] = {}
    isolate_vars: dict[str, list[tuple]] = defaultdict(list)
    with PER_ISOLATE.open(encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            gene = row["gene"]
            if gene not in TARGET_GENES:
                continue
            key = (gene, row["pos"], row["ref"], row["alt"])
            key_info.setdefault(key, row)
            isolate_vars[row["isolate_id"]].append(key)

    def is_non_syn(k: tuple) -> bool:
        return (key_info[k].get("mutation_type") or "").strip().lower() != "synonymous"

    def final_call(k: tuple) -> str:
        if k in in_server:
            return re_resolve(k[0], server_rows[k]["final_call"])
        return re_resolve(k[0], key_info[k]["resistance_label"])

    keys = {k for k in key_info if is_non_syn(k)}
    judged = {k for k in keys if k in in_server}
    uploaded_nj = {k for k in keys if k not in in_server and k in judge_keys}
    not_up = keys - judged - uploaded_nj

    def flows_by_source(src: set, call: str) -> int:
        return sum((1 for k in src if final_call(k) == call))

    r_triple = (
        flows_by_source(judged, "R"),
        flows_by_source(uploaded_nj, "R"),
        flows_by_source(not_up, "R"),
    )
    s_triple = (
        flows_by_source(judged, "S"),
        flows_by_source(uploaded_nj, "S"),
        flows_by_source(not_up, "S"),
    )
    i_triple = (
        flows_by_source(judged, "INDET"),
        flows_by_source(uploaded_nj, "INDET"),
        flows_by_source(not_up, "INDET"),
    )
    r_var = sum(r_triple)
    s_var = sum(s_triple)
    i_var = sum(i_triple)
    iso_calls: dict[str, str] = {}
    r_pairs: Counter[str] = Counter()
    for iso, ks in isolate_vars.items():
        calls = [final_call(k) for k in ks]
        if "R" in calls:
            iso_calls[iso] = "R"
            r_pairs.update((k[0] for k in ks if final_call(k) == "R"))
        elif "S" in calls:
            iso_calls[iso] = "S"
        else:
            iso_calls[iso] = "INDET"
    iso_dist = Counter(iso_calls.values())
    fs_by_gene = Counter()
    lof_pairs = Counter()
    for iso, ks in isolate_vars.items():
        for k in ks:
            if final_call(k) != "R":
                continue
            info = key_info[k]
            mt = (info.get("mutation_type") or "").strip().lower()
            aa = info.get("aa_change") or ""
            if k in in_server:
                lof_pairs[k[0]] += 1
            elif mt == "frameshift" or "fs" in aa:
                fs_by_gene[k[0]] += 1
    screening_var = {k[0] for k in keys - judged if final_call(k) == "R"}
    fs_var_gene = Counter()
    for k in keys - judged:
        if final_call(k) == "R":
            fs_var_gene[k[0]] += 1
    return {
        "n_judged": len(judged),
        "n_njudge": len(uploaded_nj),
        "n_notup": len(not_up),
        "n_var": len(keys),
        "r_var": r_var,
        "s_var": s_var,
        "i_var": i_var,
        "r_triple": r_triple,
        "s_triple": s_triple,
        "i_triple": i_triple,
        "iso": iso_dist,
        "n_carriers": len(iso_calls),
        "r_pairs": dict(r_pairs),
        "lof_pairs": dict(lof_pairs),
        "fs_pair_gene": dict(fs_by_gene),
        "screening_var_gene": dict(fs_var_gene),
        "server_r_gene": Counter((k[0] for k in judged if final_call(k) == "R")),
    }


def flow(ax, x0, w0, yt0, yb0, x1, w1, yt1, yb1, color, alpha=0.55):
    xm = (x0 + w0 + x1) / 2
    verts = [
        (x0 + w0, yt0),
        (xm, yt0),
        (x1, yt1),
        (x1, yb1),
        (xm, yb0),
        (x0 + w0, yb0),
        (0, 0),
    ]
    codes = [
        MplPath.MOVETO,
        MplPath.CURVE3,
        MplPath.CURVE3,
        MplPath.LINETO,
        MplPath.CURVE3,
        MplPath.CURVE3,
        MplPath.CLOSEPOLY,
    ]
    ax.add_patch(
        PathPatch(
            MplPath(verts, codes),
            facecolor=color,
            alpha=alpha,
            edgecolor="none",
            zorder=2,
        )
    )


def segment(
    ax,
    x,
    w,
    ytop,
    ybot,
    color,
    label,
    value,
    pct,
    sub_fractions=None,
    sub_colors=None,
    label_side="right",
    name=None,
):
    ax.add_patch(
        FancyBboxPatch(
            (x, ybot),
            w,
            ytop - ybot,
            boxstyle="round,pad=0,rounding_size=0.004",
            facecolor=color,
            edgecolor="white",
            linewidth=1.2,
            zorder=3,
        )
    )
    if sub_fractions:
        total = sum(sub_fractions)
        fracs = [f / total for f in sub_fractions]
        y = ybot
        for fr, c in zip(fracs, sub_colors):
            y1 = y + (ytop - ybot) * fr
            ax.plot(
                [x + 0.012, x + w - 0.012],
                [y1, y1],
                color=c,
                lw=1.0,
                ls=(0, (4, 2)),
                alpha=0.95,
                zorder=4,
            )
            y = y1
    ymid = (ytop + ybot) / 2
    h = ytop - ybot
    bbox_kw = dict(
        boxstyle="round,pad=0.22", fc="white", ec="#cccccc", lw=0.6, alpha=0.9
    )
    if h > 0.14:
        if name:
            ax.text(
                x + w / 2,
                ytop - 0.045,
                name,
                ha="center",
                va="center",
                fontsize=7.5,
                color="white",
                alpha=0.95,
                zorder=5,
            )
        ax.text(
            x + w / 2,
            ymid + 0.02,
            f"{value:,}",
            ha="center",
            va="center",
            fontsize=10,
            fontweight="bold",
            color="white",
            zorder=5,
        )
        ax.text(
            x + w / 2,
            ymid - 0.04,
            f"({pct:.1f}%)",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            alpha=0.9,
            zorder=5,
        )
    elif h > 0.075:
        if name:
            ax.text(
                x + w / 2,
                ymid + 0.022,
                name,
                ha="center",
                va="center",
                fontsize=7,
                color="white",
                alpha=0.95,
                zorder=5,
            )
        ax.text(
            x + w / 2,
            ymid - 0.003,
            f"{value:,}",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color="white",
            zorder=5,
        )
        ax.text(
            x + w / 2,
            ymid - 0.028,
            f"({pct:.1f}%)",
            ha="center",
            va="center",
            fontsize=7,
            color="white",
            alpha=0.9,
            zorder=5,
        )
    elif h > 0.035:
        ax.text(
            x + w / 2,
            ymid,
            f"{value:,}",
            ha="center",
            va="center",
            fontsize=9,
            fontweight="bold",
            color="white",
            zorder=5,
        )
        txt = f"{label} ({pct:.1f}%)"
        if label_side == "left":
            ax.text(
                x - 0.012,
                ymid,
                txt,
                ha="right",
                va="center",
                fontsize=7.5,
                color="#333333",
                zorder=7,
                bbox=bbox_kw,
            )
        else:
            ax.text(
                x + w + 0.012,
                ymid,
                txt,
                ha="left",
                va="center",
                fontsize=7.5,
                color="#333333",
                zorder=7,
                bbox=bbox_kw,
            )
    else:
        txt = f"{label} {value:,} ({pct:.1f}%)"
        if label_side == "left":
            ax.text(
                x - 0.012,
                ymid,
                txt,
                ha="right",
                va="center",
                fontsize=7.5,
                color="#333333",
                zorder=7,
                bbox=bbox_kw,
            )
        else:
            ax.text(
                x + w + 0.012,
                ymid,
                txt,
                ha="left",
                va="center",
                fontsize=7.5,
                color="#333333",
                zorder=7,
                bbox=bbox_kw,
            )


def main() -> None:
    d = load_flows()
    n_j, n_n, n_u = (d["n_judged"], d["n_njudge"], d["n_notup"])
    n_var = d["n_var"]
    r, s, i = (d["r_var"], d["s_var"], d["i_var"])
    iso = d["iso"]
    n_c = d["n_carriers"]
    fig, ax = plt.subplots(figsize=(13.6, 7.4))
    ax.set_xlim(-0.06, 1.04)
    ax.set_ylim(0.0, 1.14)
    ax.axis("off")
    x0, x1, x2 = (0.02, 0.46, 0.865)
    wcol = 0.1
    seg1 = [n_j, n_n, n_u]
    ys1 = []
    y = 0.0
    for v in seg1:
        ys1.append((y, y + v / n_var))
        y += v / n_var
    j_y, n_y, u_y = ys1
    seg2 = [r, s, i]
    ys2 = []
    y = 1.0
    for v in seg2:
        ys2.append((y - v / n_var, y))
        y -= v / n_var
    r_y, s_y, i_y = ys2
    rt, st, it = (d["r_triple"], d["s_triple"], d["i_triple"])

    def band(src_range, src_size, dst_range, dst_size, color, alpha=0.5):
        flow(
            ax,
            x0,
            wcol,
            src_range[0],
            src_range[1],
            x1,
            wcol,
            dst_range[0],
            dst_range[1],
            color,
            alpha,
        )

    j_rf, j_sf, j_if = (rt[0] / n_j, st[0] / n_j, it[0] / n_j)
    j_seg_r = (j_y[1], j_y[1] - (j_y[1] - j_y[0]) * j_rf)
    j_seg_s = (
        j_y[1] - (j_y[1] - j_y[0]) * j_rf,
        j_y[1] - (j_y[1] - j_y[0]) * (j_rf + j_sf),
    )
    j_seg_i = (j_y[1] - (j_y[1] - j_y[0]) * (j_rf + j_sf), j_y[0])
    n_rf, n_sf, n_if = (
        rt[1] / n_n if n_n else 0,
        st[1] / n_n if n_n else 0,
        it[1] / n_n if n_n else 0,
    )
    n_seg_r = (n_y[1], n_y[1] - (n_y[1] - n_y[0]) * n_rf)
    n_seg_s = (
        n_y[1] - (n_y[1] - n_y[0]) * n_rf,
        n_y[1] - (n_y[1] - n_y[0]) * (n_rf + n_sf),
    )
    n_seg_i = (n_y[1] - (n_y[1] - n_y[0]) * (n_rf + n_sf), n_y[0])
    u_rf, u_sf, u_if = (
        rt[2] / n_u if n_u else 0,
        st[2] / n_u if n_u else 0,
        it[2] / n_u if n_u else 0,
    )
    u_seg_r = (u_y[1], u_y[1] - (u_y[1] - u_y[0]) * u_rf)
    u_seg_s = (
        u_y[1] - (u_y[1] - u_y[0]) * u_rf,
        u_y[1] - (u_y[1] - u_y[0]) * (u_rf + u_sf),
    )
    u_seg_i = (u_y[1] - (u_y[1] - u_y[0]) * (u_rf + u_sf), u_y[0])

    def dst_bands(yseg, triple):
        total = sum(triple)
        ranges, ytop = ([], yseg[1])
        for v in triple:
            if v == 0:
                ranges.append((ytop, ytop))
                continue
            ybot = ytop - (yseg[1] - yseg[0]) * v / total
            ranges.append((ytop, ybot))
            ytop = ybot
        return ranges

    j_rb, n_rb, u_rb = dst_bands(r_y, rt)
    j_sb, n_sb, u_sb = dst_bands(s_y, st)
    j_ib, n_ib, u_ib = dst_bands(i_y, it)
    band(j_seg_r, n_j, j_rb, r, C_JUDGED, 0.5)
    band(j_seg_s, n_j, j_sb, s, C_JUDGED, 0.5)
    band(j_seg_i, n_j, j_ib, i, C_JUDGED, 0.5)
    if st[1]:
        band(n_seg_s, n_n, n_sb, s, C_NJUDGE, 0.7)
    band(n_seg_i, n_n, n_ib, i, C_NJUDGE, 0.7)
    if rt[2]:
        band(u_seg_r, n_u, u_rb, r, C_NOTUP, 0.8)
    if st[2]:
        band(u_seg_s, n_u, u_sb, s, C_NOTUP, 0.8)
    band(u_seg_i, n_u, u_ib, i, C_NOTUP, 0.8)
    segment(
        ax,
        x0,
        wcol,
        *j_y[::-1],
        C_JUDGED,
        "Full pipeline analysis",
        n_j,
        100 * n_j / n_var,
        name="Full pipeline\nanalysis",
    )
    segment(
        ax,
        x0,
        wcol,
        *n_y[::-1],
        C_NJUDGE,
        "Uploaded, not analysed",
        n_n,
        100 * n_n / n_var,
        name="Uploaded,\nnot analysed",
    )
    segment(
        ax,
        x0,
        wcol,
        *u_y[::-1],
        C_NOTUP,
        "Not uploaded",
        n_u,
        100 * n_u / n_var,
        name="Not uploaded",
    )
    tri_colors = [C_JUDGED, C_NJUDGE, C_NOTUP]
    segment(
        ax,
        x1,
        wcol,
        *r_y[::-1],
        C_R,
        "Resistant",
        r,
        100 * r / n_var,
        sub_fractions=list(reversed(rt)),
        sub_colors=list(reversed(tri_colors)),
        label_side="left",
    )
    segment(
        ax,
        x1,
        wcol,
        *s_y[::-1],
        C_S,
        "Susceptible",
        s,
        100 * s / n_var,
        sub_fractions=list(reversed(st)),
        sub_colors=list(reversed(tri_colors)),
        name="Susceptible",
    )
    segment(
        ax,
        x1,
        wcol,
        *i_y[::-1],
        C_I,
        "Indeterminate",
        i,
        100 * i / n_var,
        sub_fractions=list(reversed(it)),
        sub_colors=list(reversed(tri_colors)),
        name="Indeterminate",
    )
    iso_r, iso_s, iso_i = (iso["R"], iso["S"], iso["INDET"])
    y = 1.0
    ir_y = (y - iso_r / n_c, y)
    y -= iso_r / n_c
    is_y = (y - iso_s / n_c, y)
    y -= iso_s / n_c
    ii_y = (y - iso_i / n_c, y)
    flow(ax, x1, wcol, r_y[1], r_y[0], x2, wcol, ir_y[1], ir_y[0], C_R, 0.75)
    flow(ax, x1, wcol, s_y[1], s_y[0], x2, wcol, is_y[1], is_y[0], C_S, 0.75)
    flow(ax, x1, wcol, i_y[1], i_y[0], x2, wcol, ii_y[1], ii_y[0], C_I, 0.75)
    segment(
        ax,
        x2,
        wcol,
        *ir_y[::-1],
        C_R,
        "Resistant",
        iso_r,
        100 * iso_r / n_c,
        label_side="right",
    )
    segment(
        ax,
        x2,
        wcol,
        *is_y[::-1],
        C_S,
        "Susceptible",
        iso_s,
        100 * iso_s / n_c,
        name="Susceptible",
    )
    segment(
        ax,
        x2,
        wcol,
        *ii_y[::-1],
        C_I,
        "Indeterminate",
        iso_i,
        100 * iso_i / n_c,
        name="Indeterminate",
    )
    xgap = (x1 + wcol + x2) / 2
    for ya, yb, nv, ni in (
        (r_y, ir_y, r, iso_r),
        (s_y, is_y, s, iso_s),
        (i_y, ii_y, i, iso_i),
    ):
        ym = (ya[0] + ya[1] + (yb[0] + yb[1])) / 4
        ax.text(
            xgap,
            ym,
            f"{nv:,} variants  →  {ni:,} isolates",
            ha="center",
            va="center",
            fontsize=7.5,
            color="#333333",
            zorder=6,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none", alpha=0.8),
        )
    ax.text(
        x0 + wcol / 2,
        1.045,
        "Unique non-synonymous\nvariants",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        x0 + wcol / 2,
        1.005,
        f"(n = {n_var:,})",
        ha="center",
        va="bottom",
        fontsize=8,
        color="dimgrey",
    )
    ax.text(
        x1 + wcol / 2,
        1.045,
        "Variant-level calls",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        x1 + wcol / 2,
        1.005,
        f"(n = {n_var:,})",
        ha="center",
        va="bottom",
        fontsize=8,
        color="dimgrey",
    )
    ax.text(
        x2 + wcol / 2,
        1.045,
        "Isolate-level calls",
        ha="center",
        va="bottom",
        fontsize=10,
        fontweight="bold",
    )
    ax.text(
        x2 + wcol / 2,
        1.005,
        f"(n = {n_c:,} carriers)",
        ha="center",
        va="bottom",
        fontsize=8,
        color="dimgrey",
    )
    out = OUT_DIR / "fig_alluvial_application"
    for ext in ("png", "pdf", "svg"):
        fig.savefig(f"{out}.{ext}")
    plt.close(fig)
    print(f"Alluvial figure written to: {out}.png/.pdf/.svg")
    print(
        f"flows: judged={n_j} uploaded_nj={n_n} not_up={n_u} | R={r} S={s} INDET={i} | isolates R={iso_r} S={iso_s} INDET={iso_i} | triples R={rt} S={st} I={it}"
    )


if __name__ == "__main__":
    main()
