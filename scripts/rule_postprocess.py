from __future__ import annotations
import argparse
import csv
import sqlite3
import sys
from collections import defaultdict
from pathlib import Path

RV0678_UPGRADE_MISSENSE = {
    "Ile67Leu",
    "Leu40Met",
    "Leu40Phe",
    "Phe19Ser",
    "Arg109Pro",
    "Ala99Pro",
    "Ser68Asn",
}
RV0678_DOWNGRADE_VUS = {"Glu55Asp"}
LOF_CLASSES = {"frameshift", "stop_gained", "start_lost"}
AA3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}


def aa3_to_1(aa3: str) -> str:
    m = (aa3 or "").strip()
    if not m:
        return m
    import re

    mt = re.match("^([A-Za-z]{3})(\\d+)([A-Za-z]{3})$", m)
    if not mt:
        return m
    wt1 = AA3_TO_1.get(mt.group(1).capitalize())
    mt1 = AA3_TO_1.get(mt.group(3).capitalize())
    if wt1 is None or mt1 is None:
        return m
    return f"{wt1}{mt.group(2)}{mt1}"


def norm_aa(mutation: str) -> str:
    m = (mutation or "").strip()
    if m.startswith("p."):
        m = m[2:]
    return m


def parse_prediction(value: str) -> str:
    v = (value or "").strip().upper()
    if v in {"R", "RESISTANT"}:
        return "R"
    if v in {"S", "SENSITIVE", "SUSCEPTIBLE"}:
        return "S"
    return "UNKNOWN"


def load_evo_delta(
    db_path: str, gene: str, aa_list: set[str]
) -> dict[str, float | None]:
    result: dict[str, float | None] = {aa: None for aa in aa_list}
    try:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    except sqlite3.Error as exc:
        print(f"[rule_postprocess] WARNING: cannot open db {db_path}: {exc}")
        return result
    try:
        aa1_list = {aa3_to_1(aa) for aa in aa_list}
        placeholders = ",".join(("?" for _ in aa1_list))
        sql = f"SELECT m.aa_change, MIN(m.evo_delta) FROM mutations m JOIN genes g ON g.id = m.gene_id WHERE g.name = ? AND m.aa_change IN ({placeholders}) GROUP BY m.aa_change"
        aa1_to_3 = {aa3_to_1(aa): aa for aa in aa_list}
        for aa1, delta in conn.execute(sql, [gene, *sorted(aa1_list)]):
            if aa1 in aa1_to_3:
                result[aa1_to_3[aa1]] = delta
    finally:
        conn.close()
    return result


def apply_variant_rules(
    row: dict[str, str],
    *,
    enable_upgrade: bool = True,
    enable_downgrade_vus: bool = True,
    enable_mtra_dereg: bool = True,
    enable_pepq_reviewed: bool = True,
    gate_evo_delta: dict[str, float | None] | None = None,
    evo_threshold: float | None = None,
) -> tuple[str, str]:
    pred = parse_prediction(row.get("amr_hunter_prediction"))
    gene = (row.get("gene") or "").strip()
    mutation = (row.get("mutation") or "").strip()
    aa = norm_aa(mutation)
    vc = (row.get("variant_class") or "").strip()
    ev = (row.get("evidence_code") or "").strip()
    tag = ""
    if enable_upgrade and gene == "Rv0678" and (aa in RV0678_UPGRADE_MISSENSE):
        if gate_evo_delta is not None:
            delta = gate_evo_delta.get(aa)
            if delta is None:
                tag = "upgrade_missense:SKIPPED_NO_EVO2"
            elif evo_threshold is not None and delta > evo_threshold:
                tag = f"upgrade_missense:SKIPPED_EVO2_{delta:+.1f}"
            else:
                pred = "R"
                tag = "upgrade_missense"
        else:
            pred = "R"
            tag = "upgrade_missense"
    if (
        enable_downgrade_vus
        and gene == "Rv0678"
        and (aa in RV0678_DOWNGRADE_VUS)
        and (pred == "R")
    ):
        pred = "UNKNOWN"
        tag = "downgrade_vus"
    if enable_mtra_dereg and gene == "mtrB" and (ev == "MTRAB_DEREGULATION"):
        pred = "S"
        tag = "mtra_dereg"
    if (
        enable_pepq_reviewed
        and gene == "pepQ"
        and (ev == "NEGATIVE_LOF_REVIEWED")
        and (vc == "protein_substitution")
    ):
        pred = "UNKNOWN"
        tag = "pepq_reviewed"
    return (pred, tag)


def recompute_isolate_calls(variant_rows: list[dict[str, str]]) -> dict[str, str]:
    by_iso: dict[str, list[str]] = defaultdict(list)
    for r in variant_rows:
        by_iso[r["unique_id"]].append(parse_prediction(r["amr_hunter_prediction"]))
    calls: dict[str, str] = {}
    for uid, ps in by_iso.items():
        if "R" in ps:
            calls[uid] = "R"
        elif "S" in ps:
            calls[uid] = "S"
        else:
            calls[uid] = "UNKNOWN"
    return calls


def confusion(iso_rows: list[dict[str, str]], calls: dict[str, str]) -> dict:
    tp = tn = fp = fn = unc = 0
    for r in iso_rows:
        u = r["unique_id"]
        p = parse_prediction(r.get("bdq_binary_phenotype"))
        c = calls.get(u, "UNKNOWN")
        if p not in ("R", "S"):
            continue
        if c == "UNKNOWN":
            unc += 1
        elif p == "R":
            tp += 1 if c == "R" else 0
            fn += 0 if c == "R" else 1
        else:
            fp += 1 if c == "R" else 0
            tn += 0 if c == "R" else 1
    return dict(tp=tp, tn=tn, fp=fp, fn=fn, unc=unc)


def fmt_metrics(m: dict) -> str:
    sen = m["tp"] / (m["tp"] + m["fn"]) if m["tp"] + m["fn"] else 0.0
    spe = m["tn"] / (m["tn"] + m["fp"]) if m["tn"] + m["fp"] else 0.0
    ppv = m["tp"] / (m["tp"] + m["fp"]) if m["tp"] + m["fp"] else 0.0
    return f"TP={m['tp']:3d} FP={m['fp']:3d} FN={m['fn']:3d} TN={m['tn']:5d} UNK={m['unc']:3d} | sens={sen * 100:5.1f}% spec={spe * 100:6.2f}% ppv={ppv * 100:5.1f}%"


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description="CRyPTIC 判定后处理规则层")
    ap.add_argument("--variant-csv", required=True, help="变体级判定 CSV")
    ap.add_argument("--isolate-csv", required=True, help="株级判定 CSV")
    ap.add_argument("--out-dir", required=True, help="输出目录")
    ap.add_argument(
        "--skip",
        default="",
        help="逗号分隔要禁用的规则: upgrade,downgrade_vus,mtra_dereg,pepq_reviewed",
    )
    ap.add_argument("--gate-evo2", default="", help="BDQ-SatDB 路径（启用 Evo2 门槛）")
    ap.add_argument(
        "--evo-threshold",
        type=float,
        default=-8.0,
        help="Evo2 evo_delta 升级门槛（默认 -8）",
    )
    args = ap.parse_args(argv)
    skip = {s.strip() for s in args.skip.split(",") if s.strip()}
    enable = {
        "upgrade": "upgrade" not in skip,
        "downgrade_vus": "downgrade_vus" not in skip,
        "mtra_dereg": "mtra_dereg" not in skip,
        "pepq_reviewed": "pepq_reviewed" not in skip,
    }
    var_path = Path(args.variant_csv)
    iso_path = Path(args.isolate_csv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with var_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        var_rows = list(reader)
        var_fields = list(reader.fieldnames or [])
    with iso_path.open(newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        iso_rows = list(reader)
        iso_fields = list(reader.fieldnames or [])
    gate = None
    if args.gate_evo2 and enable["upgrade"]:
        gate = load_evo_delta(args.gate_evo2, "Rv0678", RV0678_UPGRADE_MISSENSE)
        missing = [aa for aa, d in gate.items() if d is None]
        if missing:
            print(
                f"[rule_postprocess] Evo2 门槛：{len(missing)}/{len(gate)} 个候选在数据库中无记录，将保守不升级: {sorted(missing)}"
            )
        for aa, d in sorted(gate.items()):
            if d is not None:
                status = "通过" if d <= args.evo_threshold else "未通过"
                print(f"[rule_postprocess] {aa}: evo_delta={d:+.2f} ({status})")
    tag_counts: dict[str, int] = defaultdict(int)
    new_variant_fields = (
        var_fields + ["postprocess_rule"]
        if "postprocess_rule" not in var_fields
        else var_fields
    )
    out_rows = []
    for r in var_rows:
        new_pred, tag = apply_variant_rules(
            r,
            enable_upgrade=enable["upgrade"],
            enable_downgrade_vus=enable["downgrade_vus"],
            enable_mtra_dereg=enable["mtra_dereg"],
            enable_pepq_reviewed=enable["pepq_reviewed"],
            gate_evo_delta=gate,
            evo_threshold=args.evo_threshold if args.gate_evo2 else None,
        )
        if tag:
            tag_counts[tag] += 1
        r = dict(r)
        r["amr_hunter_prediction"] = new_pred
        r["postprocess_rule"] = tag
        out_rows.append(r)
    baseline_variant = []
    with var_path.open(newline="", encoding="utf-8-sig") as f:
        baseline_variant = list(csv.DictReader(f))
    base_calls = recompute_isolate_calls(baseline_variant)
    new_calls = recompute_isolate_calls(out_rows)
    print(f"[rule_postprocess] 规则变更计数: {dict(tag_counts) or '无变化'}")
    print(f"[rule_postprocess] 基线:   {fmt_metrics(confusion(iso_rows, base_calls))}")
    print(f"[rule_postprocess] 应用后: {fmt_metrics(confusion(iso_rows, new_calls))}")
    out_var = out_dir / var_path.name
    with out_var.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=new_variant_fields)
        w.writeheader()
        for r in out_rows:
            w.writerow({k: r.get(k, "") for k in new_variant_fields})
    print(f"[rule_postprocess] 变体级输出: {out_var}")
    out_iso = out_dir / iso_path.name
    with out_iso.open("w", newline="", encoding="utf-8-sig") as f:
        reader2 = csv.DictReader(iso_path.open(newline="", encoding="utf-8-sig"))
        iso_fields2 = list(reader2.fieldnames or [])
        iso_rows2 = list(reader2)
    changed = 0
    for r in iso_rows2:
        old = parse_prediction(r.get("amr_hunter_isolate_prediction"))
        new = new_calls.get(r["unique_id"], "UNKNOWN")
        r = dict(r)
        r["amr_hunter_isolate_prediction"] = new
        r["postprocess_rule"] = f"{old}->{new}" if old != new else ""
        if old != new:
            changed += 1
    with out_iso.open("w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=iso_fields2 + ["postprocess_rule"])
        w.writeheader()
        for r in iso_rows2:
            w.writerow({k: r.get(k, "") for k in iso_fields2 + ["postprocess_rule"]})
    print(f"[rule_postprocess] 株级输出: {out_iso} (判定变化 {changed} 株)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
