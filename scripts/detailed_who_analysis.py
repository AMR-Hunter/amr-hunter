from __future__ import annotations
import re
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple
import pandas as pd
import yaml

DB_PATH = Path("/root/gpufree-data/db/amr_hunter.db")
XLSX_PATH = Path(
    "/root/gpufree-share/amr_hunter_workspace/amr_hunter_data/reference/who/WHO-UCN-TB-2023.5-eng.xlsx"
)
CONFIG_PATH = Path("/root/gpufree-share/amr_hunter_workspace/code/config/config.yaml")
DRUG = "Bedaquiline"
AA3 = {
    "A": "Ala",
    "R": "Arg",
    "N": "Asn",
    "D": "Asp",
    "C": "Cys",
    "Q": "Gln",
    "E": "Glu",
    "G": "Gly",
    "H": "His",
    "I": "Ile",
    "L": "Leu",
    "K": "Lys",
    "M": "Met",
    "F": "Phe",
    "P": "Pro",
    "S": "Ser",
    "T": "Thr",
    "W": "Trp",
    "Y": "Tyr",
    "V": "Val",
    "*": "*",
}
SIMPLE_CDS_RE = re.compile("c\\.(\\d+)([ACGT]+)>([ACGT]+)")
UPSTREAM_RE = re.compile("c\\.\\-(\\d+)([ACGT]+)>([ACGT]+)")
GENERIC_NON_EXEC = {"LoF"}
ALL_MODELED = {"atpE", "Rv0678", "pepQ", "mmpL5", "mmpS5"}


def grade_key(grading: str) -> int:
    g = str(grading).strip()
    m = re.match("(\\d)\\)", g)
    if m:
        return int(m.group(1))
    if "Assoc. w. R" in g:
        return 1
    if "Assoc. w. R - interim" in g:
        return 2
    return 99


def truth_label(grading: str) -> Optional[str]:
    g = grade_key(grading)
    if g in (1, 2):
        return "R"
    if g in (4, 5):
        return "S"
    return None


def is_strict_comparable(mutation: str, effect: str) -> bool:
    mutation = str(mutation).strip()
    effect = str(effect).strip()
    if effect == "upstream_gene_variant":
        return False
    if mutation in GENERIC_NON_EXEC:
        return False
    if mutation.startswith("p."):
        return all((tok not in mutation for tok in ("del", "ins", "dup", "fs", "?")))
    return bool(SIMPLE_CDS_RE.fullmatch(mutation))


def classify_support(gene: str, mutation: str, effect: str) -> Tuple[bool, str]:
    mutation = str(mutation).strip()
    effect = str(effect).strip()
    if mutation in GENERIC_NON_EXEC:
        return (False, "generic_LoF_category")
    if mutation.startswith("p.") and any(
        (tok in mutation for tok in ("del", "ins", "dup", "fs", "?"))
    ):
        return (False, "non-SNP_variant")
    if effect == "upstream_gene_variant":
        return (True, "")
    if mutation.startswith("p."):
        return (True, "")
    if SIMPLE_CDS_RE.fullmatch(mutation):
        return (True, "")
    return (False, "unsupported_hgvs")


def aa_change_to_mutations(change: str) -> Set[str]:
    if not change or len(change) < 3:
        return set()
    if change.startswith("p."):
        if change.endswith("Ter"):
            return {change, f"{change[:-3]}*"}
        if change.endswith("*"):
            return {change, f"{change[:-1]}Ter"}
        return {change}
    wt, mt, pos = (change[0], change[-1], change[1:-1])
    if wt not in AA3 or mt not in AA3 or (not pos.isdigit()):
        return set()
    mutation = f"p.{AA3[wt]}{int(pos)}{AA3[mt]}"
    if mt == "*":
        return {mutation, f"p.{AA3[wt]}{int(pos)}Ter"}
    return {mutation}


def load_who() -> pd.DataFrame:
    who = pd.read_excel(
        XLSX_PATH,
        sheet_name="Catalogue_master_file",
        header=2,
        usecols=["drug", "gene", "mutation", "effect", "FINAL CONFIDENCE GRADING"],
    )
    who = who[who["drug"] == DRUG].copy()
    who["grade"] = who["FINAL CONFIDENCE GRADING"].apply(grade_key)
    who["truth"] = who["FINAL CONFIDENCE GRADING"].apply(truth_label)
    who["mutation"] = who["mutation"].astype(str).str.strip()
    who["supported"], who["gap_reason"] = zip(
        *who.apply(
            lambda r: classify_support(r["gene"], r["mutation"], r["effect"]), axis=1
        )
    )
    return who


def load_db() -> Tuple[
    Dict[Tuple[str, str], List[Dict]],
    Dict[Tuple[str, str, int, str, str], List[Dict]],
]:
    conn = sqlite3.connect(DB_PATH)
    conn.row_factory = sqlite3.Row
    genes_list = list(ALL_MODELED)
    rows = conn.execute(
        "SELECT g.name AS gene, g.logic_type, m.pos, m.ref, m.alt, m.aa_change,\n                  m.region_type, m.region_name, m.status, m.final_interpretation,\n                  m.resistance_phenotype\n           FROM mutations m JOIN genes g ON g.id = m.gene_id\n           WHERE g.name IN ({})".format(
            ",".join(("?" for _ in genes_list))
        ),
        genes_list,
    ).fetchall()
    conn.close()
    cds: Dict[Tuple[str, str], List[Dict]] = defaultdict(list)
    reg: Dict[Tuple[str, str, int, str, str], List[Dict]] = defaultdict(list)
    for row in rows:
        gene = str(row["gene"])
        rt = str(row["region_type"] or "CDS").upper()
        rec = {
            "logic_type": str(row["logic_type"] or "POSITIVE").upper(),
            "status": str(row["status"] or ""),
            "final_interpretation": str(row["final_interpretation"] or ""),
            "resistance_phenotype": str(row["resistance_phenotype"] or ""),
        }
        if rt == "REGULATORY":
            reg[
                gene,
                str(row["region_name"] or "REG"),
                int(row["pos"]),
                str(row["ref"]),
                str(row["alt"]),
            ].append(rec)
        else:
            cds[gene, f"c.{row['pos']}{row['ref']}>{row['alt']}"].append(rec)
            for mut in aa_change_to_mutations(row["aa_change"]):
                cds[gene, mut].append(rec)
    return (cds, reg)


def derive_pred(records: List[Dict]) -> Optional[str]:
    labels: Set[str] = set()
    for r in records:
        sp = str(r.get("resistance_phenotype") or "").strip().upper()
        if sp == "RESISTANT":
            labels.add("R")
            continue
        if sp == "SENSITIVE":
            labels.add("S")
            continue
        if sp == "UNKNOWN":
            continue
        status = str(r.get("status") or "")
        interp = str(r.get("final_interpretation") or "").strip().lower()
        lt = str(r.get("logic_type") or "POSITIVE").upper()
        if status == "RE_SENSITIZED":
            labels.add("S")
            continue
        if interp.startswith("sensitive:"):
            labels.add("S")
            continue
        if status == "LOSS_OF_FUNCTION":
            labels.add("S" if lt == "STRUCTURAL" else "R")
            continue
        if status == "COMPLETED":
            if lt in ("NEGATIVE", "STRUCTURAL"):
                labels.add("S")
                continue
            if interp.startswith("resistant:"):
                labels.add("R")
                continue
    return next(iter(labels)) if len(labels) == 1 else None


def run():
    who = load_who()
    cds, reg = load_db()
    who_m = who[who["gene"].isin(ALL_MODELED)].copy()
    who_eval = who_m[who_m["truth"].notna()].copy()
    print("=" * 80)
    print("WHO BDQ 详细逐级分析报告")
    print("=" * 80)
    print(f"\nWHO BDQ 总条目: {len(who)}")
    print(f"5基因 (atpE,Rv0678,pepQ,mmpL5,mmpS5) 条目: {len(who_m)}")
    print(f"其中 G3(不确定) 条目: {len(who_m[who_m['truth'].isna()])}")
    print(f"可评估 (G1+G2+G4+G5) 条目: {len(who_eval)}")
    print("\n" + "=" * 80)
    print("一、各基因 WHO Grade 分布总览")
    print("=" * 80)
    for gene in sorted(ALL_MODELED):
        gdf = who[who["gene"] == gene]
        g1 = len(gdf[gdf["grade"] == 1])
        g2 = len(gdf[gdf["grade"] == 2])
        g3 = len(gdf[gdf["grade"] == 3])
        g4 = len(gdf[gdf["grade"] == 4])
        g5 = len(gdf[gdf["grade"] == 5])
        total = len(gdf)
        evaluable = g1 + g2 + g4 + g5
        print(f"\n  {gene}: WHO总={total}, G1={g1}, G2={g2}, G3={g3}, G4={g4}, G5={g5}")
        print(f"         可评估(G1+G2+G4+G5)={evaluable}")
    print("\n" + "=" * 80)
    print("二、各基因逐Grade 匹配详情")
    print("=" * 80)
    all_matched = []
    all_absent = []
    for gene in sorted(ALL_MODELED):
        gdf = who_eval[who_eval["gene"] == gene].copy()
        print(f"\n{'─' * 60}")
        print(f"  【{gene}】 可评估条目: {len(gdf)}")
        print(f"{'─' * 60}")
        for grade_num in [1, 2, 4, 5]:
            grade_label = f"G{grade_num}"
            g_rows = gdf[gdf["grade"] == grade_num]
            if len(g_rows) == 0:
                continue
            truth = "R" if grade_num <= 2 else "S"
            matched_correct = 0
            matched_wrong = 0
            absent_list = []
            matched_detail = []
            for _, row in g_rows.iterrows():
                mutation = str(row["mutation"])
                effect = str(row["effect"] or "")
                grade_str = str(row["FINAL CONFIDENCE GRADING"])
                records = cds.get((gene, mutation), [])
                if effect == "upstream_gene_variant" and (not records):
                    parsed = UPSTREAM_RE.fullmatch(mutation)
                    if parsed:
                        up_dist = int(parsed.group(1))
                        ref_nt = parsed.group(2)
                        alt_nt = parsed.group(3)
                        for (g, rn, pos, ref, alt), recs in reg.items():
                            if g == gene:
                                records.extend(recs)
                pred = derive_pred(records)
                if pred in ("R", "S"):
                    correct = pred == truth
                    if correct:
                        matched_correct += 1
                    else:
                        matched_wrong += 1
                    matched_detail.append(
                        {
                            "mutation": mutation,
                            "truth": truth,
                            "pred": pred,
                            "correct": correct,
                            "grade_str": grade_str,
                            "statuses": [r["status"] for r in records],
                        }
                    )
                else:
                    reason = "DB中无此突变记录"
                    is_supported, gap = classify_support(gene, mutation, effect)
                    if not is_supported:
                        reason = f"不支持({gap})"
                    elif records:
                        statuses = set((r["status"] for r in records))
                        reason = f"有记录但无最终判定(状态:{','.join(statuses)})"
                    absent_list.append(
                        {
                            "mutation": mutation,
                            "effect": effect,
                            "grade_str": grade_str,
                            "reason": reason,
                        }
                    )
            total_grade = len(g_rows)
            matched_total = matched_correct + matched_wrong
            print(
                f"\n    {grade_label} (真值={truth}): WHO条目={total_grade}, DB匹配={matched_total}, 匹配正确={matched_correct}, 匹配错误={matched_wrong}, 未匹配={len(absent_list)}"
            )
            if matched_detail:
                for d in matched_detail:
                    icon = "✅" if d["correct"] else "❌"
                    print(
                        f"      {icon} {d['mutation']} | truth={d['truth']} pred={d['pred']} | status={d['statuses']}"
                    )
            if absent_list:
                for a in absent_list:
                    print(f"      ⬜ {a['mutation']} ({a['effect']}) — {a['reason']}")
            all_matched.extend(matched_detail)
            all_absent.extend(absent_list)
    print("\n" + "=" * 80)
    print("三、各基因汇总统计")
    print("=" * 80)
    print(
        f"\n{'基因':<10} {'G1':>4} {'G2':>4} {'G4':>4} {'G5':>4} {'可评估':>6} {'匹配':>4} {'正确R':>5} {'正确S':>5} {'错误':>4} {'未匹配':>5} {'准确率':>8}"
    )
    print("-" * 70)
    grand_total_eval = 0
    grand_matched = 0
    grand_correct_r = 0
    grand_correct_s = 0
    grand_wrong = 0
    for gene in sorted(ALL_MODELED):
        gdf_all = who[who["gene"] == gene]
        gdf = who_eval[who_eval["gene"] == gene]
        g1 = len(gdf_all[gdf_all["grade"] == 1])
        g2 = len(gdf_all[gdf_all["grade"] == 2])
        g3 = len(gdf_all[gdf_all["grade"] == 3])
        g4 = len(gdf_all[gdf_all["grade"] == 4])
        g5 = len(gdf_all[gdf_all["grade"] == 5])
        total_eval = len(gdf)
        matched = 0
        correct_r = 0
        correct_s = 0
        wrong = 0
        for _, row in gdf.iterrows():
            mutation = str(row["mutation"])
            effect = str(row["effect"] or "")
            truth = str(row["truth"])
            records = cds.get((gene, mutation), [])
            pred = derive_pred(records)
            if pred in ("R", "S"):
                matched += 1
                if pred == truth:
                    if truth == "R":
                        correct_r += 1
                    else:
                        correct_s += 1
                else:
                    wrong += 1
        absent = total_eval - matched
        acc = f"{matched / matched * 100:.0f}%" if matched > 0 else "N/A"
        if correct_r + correct_s == matched:
            acc = "100%"
        print(
            f"{gene:<10} {g1:>4} {g2:>4} {g4:>4} {g5:>4} {total_eval:>6} {matched:>4} {correct_r:>5} {correct_s:>5} {wrong:>4} {absent:>5} {acc:>8}"
        )
        grand_total_eval += total_eval
        grand_matched += matched
        grand_correct_r += correct_r
        grand_correct_s += correct_s
        grand_wrong += wrong
    print("-" * 70)
    grand_absent = grand_total_eval - grand_matched
    print(
        f"{'合计':<10} {'':>4} {'':>4} {'':>4} {'':>4} {grand_total_eval:>6} {grand_matched:>4} {grand_correct_r:>5} {grand_correct_s:>5} {grand_wrong:>4} {grand_absent:>5} {('100%' if grand_wrong == 0 else f'{(grand_correct_r + grand_correct_s) / grand_matched * 100:.1f}%'):>8}"
    )
    print(
        f"\n覆盖率: {grand_matched}/{grand_total_eval} = {grand_matched / grand_total_eval * 100:.1f}%"
    )
    print("\n" + "=" * 80)
    print("四、未匹配条目原因分析")
    print("=" * 80)
    absent_by_gene: Dict[str, List] = defaultdict(list)
    for a in all_absent:
        absent_by_gene[
            a["mutation"].split("_")[0] if "_" in a["mutation"] else "?"
        ].append(a)
    absent_by_gene2: Dict[str, List] = defaultdict(list)
    for _, row in who_eval.iterrows():
        gene = str(row["gene"])
        mutation = str(row["mutation"])
        effect = str(row["effect"] or "")
        records = cds.get((gene, mutation), [])
        pred = derive_pred(records)
        if pred not in ("R", "S"):
            is_supported, gap = classify_support(gene, mutation, effect)
            absent_by_gene2[gene].append(
                {
                    "mutation": mutation,
                    "effect": effect,
                    "grade": grade_key(row["FINAL CONFIDENCE GRADING"]),
                    "reason": (
                        gap
                        if not is_supported
                        else "有记录无判定" if records else "DB无此突变"
                    ),
                }
            )
    for gene in sorted(absent_by_gene2.keys()):
        entries = absent_by_gene2[gene]
        reasons = Counter((e["reason"] for e in entries))
        print(f"\n  {gene}: {len(entries)} 条未匹配")
        for reason, count in reasons.most_common():
            print(f"    - {reason}: {count} 条")
        for e in entries[:5]:
            print(f"      例: {e['mutation']} (G{e['grade']}, {e['effect']})")
    print("\n" + "=" * 80)
    print("分析完成。")
    print("=" * 80)


if __name__ == "__main__":
    run()
