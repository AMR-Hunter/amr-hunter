import csv
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Set, Optional, Tuple

BDQ_GENES = {
    "atpE": "POSITIVE",
    "Rv0678": "NEGATIVE",
    "mmpL5": "STRUCTURAL",
    "mmpS5": "STRUCTURAL",
    "pepQ": "NEGATIVE",
    "mtrA": "STRUCTURAL",
    "mtrB": "STRUCTURAL",
    "Rv1979c": "NEGATIVE",
    "glpK": "NEGATIVE",
    "lpqB": "STRUCTURAL",
}
LOF_TYPES = {
    "frameshift",
    "nonsense",
    "start_lost",
    "Insert_frameshift",
    "Delete_frameshift",
    "stop_gained",
}
KNOWN_SENSITIVE_MARKERS = {
    "mmpL5:p.Ile948Val",
    "mmpL5:p.Thr794Ile",
    "mmpL5:p.Asp767Asn",
    "mtrB:p.Met517Leu",
    "mtrB:p.Pro18Ser",
}
BDQ_LOF_NEUTRAL_GENES = {"glpK", "Rv1979c"}
BDQ_RV0678_NEUTRAL_FS_POSITIONS = {64, 71, 78, 10, 44, 124, 140, 164}


def normalize_mutation_type(raw: str) -> str:
    raw = str(raw or "").strip().lower()
    if "frameshift" in raw or "fs" in raw:
        return "frameshift"
    if raw in {"nonsense", "stop_gained", "stop-gained"}:
        return "nonsense"
    if raw in {"start_lost", "start-lost"}:
        return "start_lost"
    if raw == "synonymous":
        return "synonymous"
    if raw == "nonsynonymous":
        return "nonsynonymous"
    if raw == "intergenic":
        return "intergenic"
    return raw


def normalize_aa_change(raw: str) -> str:
    raw = str(raw or "").strip()
    if not raw or raw == "-":
        return ""
    if raw.startswith("p."):
        return raw
    parts = raw.split("-")
    if len(parts) == 2 and len(parts[0]) == 1 and (len(parts[1]) == 1):
        return f"p.{parts[0]}{parts[1]}"
    return raw


def extract_aa_position(aa_change: str) -> Optional[int]:
    match = re.search("(\\d+)", str(aa_change))
    if match:
        return int(match.group(1))
    return None


def extract_isolate_ids(info_str: str) -> Dict[str, Tuple[int, int, float]]:
    isolates = {}
    if not info_str or info_str == "-":
        return isolates
    parts = re.split('","|\\",\\"', info_str)
    for part in parts:
        part = part.strip().strip('"')
        if not part:
            continue
        match = re.match("([\\w\\-\\.]+):(\\d+):(\\d+),(\\d+):([\\d\\.]+)", part)
        if match:
            iso_id = match.group(1)
            depth = int(match.group(2))
            ref_reads = int(match.group(3))
            alt_reads = int(match.group(4))
            freq = float(match.group(5))
            isolates[iso_id] = (ref_reads, alt_reads, freq)
            continue
        match = re.match("([\\w\\-\\.]+):(\\d+):(\\d+):([\\d\\.]+)", part)
        if match:
            iso_id = match.group(1)
            ref = int(match.group(3))
            alt = 0
            freq = float(match.group(4))
            isolates[iso_id] = (ref, alt, freq)
    return isolates


def parse_mutation_classification(raw: str) -> str:
    return str(raw or "").strip()


def infer_resistance_label(
    gene: str,
    logic_type: str,
    mutation_type: str,
    aa_change: str,
    aa_pos: Optional[int],
) -> Tuple[str, str, str]:
    gene = str(gene or "").strip()
    logic_type = str(logic_type or "").strip().upper()
    if mutation_type == "synonymous":
        return ("SENSITIVE", "SYNONYMOUS_DEFAULT", "同义突变不改变蛋白序列")
    mutation_key = f"{gene}:{aa_change}"
    if mutation_key in KNOWN_SENSITIVE_MARKERS:
        return ("SENSITIVE", "KNOWN_SENSITIVE_MARKER", f"已知敏感标记: {mutation_key}")
    if mutation_type in {"intergenic", "regulatory"}:
        return ("INDETERMINATE", "NEEDS_EVO_REGULATORY", "调控区突变需要Evo2评分")
    if mutation_type in LOF_TYPES:
        if gene in BDQ_LOF_NEUTRAL_GENES:
            return (
                "INDETERMINATE",
                "EXPLORATORY_ONLY",
                f"{gene} LoF经CRyPTIC校准为BDQ中性，方向未定",
            )
        if gene == "Rv0678" and mutation_type == "frameshift":
            if aa_pos is not None and aa_pos in BDQ_RV0678_NEUTRAL_FS_POSITIONS:
                return (
                    "SENSITIVE",
                    "RV0678_NEUTRAL_FS",
                    f"Rv0678 aa{aa_pos}fs经CRyPTIC校准为中性",
                )
        if logic_type == "NEGATIVE":
            return ("RESISTANT", "NEGATIVE_LOF_SIGNATURE", f"{gene} LoF(NEGATIVE)→耐药")
        elif logic_type == "STRUCTURAL":
            return (
                "SENSITIVE",
                "STRUCTURAL_LOF_SIGNATURE",
                f"{gene} LoF(STRUCTURAL)→敏感",
            )
        elif logic_type == "POSITIVE":
            return (
                "SENSITIVE",
                "POSITIVE_TARGET_COLLAPSE",
                f"{gene} LoF(POSITIVE)→靶标崩溃敏感",
            )
    if mutation_type == "nonsynonymous":
        return (
            "INDETERMINATE",
            "NEEDS_EVO_BOLTZ",
            f"{gene}错义突变需要Evo2评分+Boltz结构复核",
        )
    return ("INDETERMINATE", "UNKNOWN", f"无法判定: {mutation_type}")


def main():
    input_csv = Path("c:/Users/zhang/Desktop/merged_all_snps.csv")
    output_dir = Path("c:/Users/zhang/Desktop/amr_hunter_local_analysis")
    output_dir.mkdir(exist_ok=True)
    isolate_variants_file = output_dir / "per_isolate_variants.csv"
    isolate_summary_file = output_dir / "isolate_summary.csv"
    gene_summary_file = output_dir / "gene_summary.csv"
    print(f"Reading: {input_csv}")
    bdq_snps = []
    all_isolates = set()
    encoding_attempts = ["gbk", "gb2312", "gb18030", "utf-8", "utf-8-sig", "latin-1"]
    lines = None
    used_encoding = None
    for enc in encoding_attempts:
        try:
            with open(input_csv, "r", encoding=enc) as f:
                lines = f.readlines()
            used_encoding = enc
            break
        except (UnicodeDecodeError, UnicodeError):
            continue
    if lines is None:
        print("ERROR: Cannot decode CSV file with any encoding")
        sys.exit(1)
    print(f"Using encoding: {used_encoding}, total lines: {len(lines)}")
    COL_POS = 0
    COL_REF = 1
    COL_ALT = 2
    COL_GENE = 3
    COL_RV = 4
    COL_GENE_POS = 5
    COL_AA_POS = 6
    COL_CODON = 7
    COL_AA_CHANGE = 8
    COL_MUT_TYPE = 9
    COL_PRODUCT = 10
    COL_FUNC = 11
    COL_ISOLATE_COUNT = 12
    COL_ISOLATE_INFO = 13
    RV_TO_GENE = {
        "Rv0676c": "mmpL5",
        "Rv1305": "atpE",
        "Rv0678": "Rv0678",
        "Rv1979c": "Rv1979c",
        "Rv2535c": "pepQ",
        "Rv2333c": "mmpS5",
        "Rv2334": "mmpL5",
        "Rv3245c": "mtrB",
        "Rv3246c": "mtrA",
        "Rv3696c": "glpK",
        "Rv3240c": "lpqB",
    }
    reader = csv.reader(lines[1:])
    for row in reader:
        if len(row) < 14:
            continue
        gene_name = row[COL_GENE].strip()
        rv_id = row[COL_RV].strip()
        resolved_gene = gene_name
        if resolved_gene not in BDQ_GENES:
            resolved_gene = RV_TO_GENE.get(rv_id, "")
            if resolved_gene not in BDQ_GENES:
                if rv_id in BDQ_GENES:
                    resolved_gene = rv_id
                else:
                    for g in BDQ_GENES:
                        if g in rv_id or g in gene_name:
                            resolved_gene = g
                            break
        if resolved_gene not in BDQ_GENES:
            continue
        mutation_type_raw = row[COL_MUT_TYPE].strip()
        mutation_type = normalize_mutation_type(mutation_type_raw)
        aa_change_raw = row[COL_AA_CHANGE].strip()
        aa_change = normalize_aa_change(aa_change_raw)
        aa_pos_str = row[COL_AA_POS].strip()
        aa_pos = None
        try:
            aa_pos = int(aa_pos_str) if aa_pos_str and aa_pos_str != "-" else None
        except ValueError:
            aa_pos = extract_aa_position(aa_change)
        isolate_info = row[COL_ISOLATE_INFO].strip()
        isolates = extract_isolate_ids(isolate_info)
        if not isolates:
            continue
        for iso_id in isolates:
            all_isolates.add(iso_id)
        bdq_snps.append(
            {
                "gene": resolved_gene,
                "rv_id": rv_id,
                "pos": row[COL_POS].strip(),
                "ref": row[COL_REF].strip(),
                "alt": row[COL_ALT].strip(),
                "gene_pos": row[COL_GENE_POS].strip(),
                "aa_change": aa_change,
                "aa_pos": aa_pos,
                "mutation_type": mutation_type,
                "mutation_type_raw": mutation_type_raw,
                "gene_product": (
                    row[COL_PRODUCT].strip() if len(row) > COL_PRODUCT else ""
                ),
                "isolates": isolates,
            }
        )
    print(f"BDQ gene SNPs: {len(bdq_snps)}")
    print(f"Unique isolates with BDQ variants: {len(all_isolates)}")
    isolate_variants = defaultdict(list)
    for snp in bdq_snps:
        gene = snp["gene"]
        logic_type = BDQ_GENES.get(gene, "UNKNOWN")
        for iso_id, (ref_r, alt_r, freq) in snp["isolates"].items():
            label, evidence, rationale = infer_resistance_label(
                gene=gene,
                logic_type=logic_type,
                mutation_type=snp["mutation_type"],
                aa_change=snp["aa_change"],
                aa_pos=snp["aa_pos"],
            )
            isolate_variants[iso_id].append(
                {
                    "gene": gene,
                    "logic_type": logic_type,
                    "aa_change": snp["aa_change"],
                    "mutation_type": snp["mutation_type"],
                    "pos": snp["pos"],
                    "ref": snp["ref"],
                    "alt": snp["alt"],
                    "freq": freq,
                    "resistance_label": label,
                    "evidence_code": evidence,
                    "rationale": rationale,
                }
            )
    variant_fields = [
        "isolate_id",
        "gene",
        "logic_type",
        "aa_change",
        "mutation_type",
        "pos",
        "ref",
        "alt",
        "freq",
        "resistance_label",
        "evidence_code",
        "rationale",
    ]
    with open(isolate_variants_file, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=variant_fields)
        writer.writeheader()
        for iso_id in sorted(isolate_variants.keys()):
            for var in isolate_variants[iso_id]:
                row = {"isolate_id": iso_id, **var}
                writer.writerow(row)
    print(f"Per-isolate variants: {isolate_variants_file}")
    summary_fields = [
        "isolate_id",
        "amr_hunter_label",
        "num_r",
        "num_likely_r",
        "num_indeterminate",
        "num_s",
        "num_variants",
        "r_genes",
        "indeterminate_genes",
        "top_rationale",
    ]
    with open(isolate_summary_file, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=summary_fields)
        writer.writeheader()
        label_counts = defaultdict(int)
        for iso_id in sorted(isolate_variants.keys()):
            variants = isolate_variants[iso_id]
            r_variants = [v for v in variants if v["resistance_label"] == "RESISTANT"]
            s_variants = [v for v in variants if v["resistance_label"] == "SENSITIVE"]
            indet_variants = [
                v for v in variants if v["resistance_label"] == "INDETERMINATE"
            ]
            if r_variants:
                agg_label = "R"
            elif indet_variants:
                agg_label = "INDETERMINATE"
            elif s_variants:
                agg_label = "S"
            else:
                agg_label = "NO_VARIANT"
            label_counts[agg_label] += 1
            r_genes = sorted(set((v["gene"] for v in r_variants)))
            indet_genes = sorted(set((v["gene"] for v in indet_variants)))
            top_rationale = (
                r_variants[0]["rationale"]
                if r_variants
                else indet_variants[0]["rationale"] if indet_variants else ""
            )
            writer.writerow(
                {
                    "isolate_id": iso_id,
                    "amr_hunter_label": agg_label,
                    "num_r": len(r_variants),
                    "num_likely_r": 0,
                    "num_indeterminate": len(indet_variants),
                    "num_s": len(s_variants),
                    "num_variants": len(variants),
                    "r_genes": ";".join(r_genes),
                    "indeterminate_genes": ";".join(indet_genes),
                    "top_rationale": top_rationale[:200],
                }
            )
    print(f"Isolate summary: {isolate_summary_file}")
    gene_stats = defaultdict(
        lambda: {
            "total_snps": 0,
            "isolates_affected": 0,
            "r_count": 0,
            "s_count": 0,
            "indet_count": 0,
            "isolate_set": set(),
        }
    )
    for snp in bdq_snps:
        gene = snp["gene"]
        gene_stats[gene]["total_snps"] += 1
        for iso_id in snp["isolates"]:
            gene_stats[gene]["isolate_set"].add(iso_id)
    for iso_id in isolate_variants:
        for var in isolate_variants[iso_id]:
            gene = var["gene"]
            if var["resistance_label"] == "RESISTANT":
                gene_stats[gene]["r_count"] += 1
            elif var["resistance_label"] == "SENSITIVE":
                gene_stats[gene]["s_count"] += 1
            elif var["resistance_label"] == "INDETERMINATE":
                gene_stats[gene]["indet_count"] += 1
    gene_fields = [
        "gene",
        "logic_type",
        "total_snps",
        "isolates_affected",
        "r_calls",
        "s_calls",
        "indet_calls",
    ]
    with open(gene_summary_file, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=gene_fields)
        writer.writeheader()
        for gene in sorted(gene_stats.keys()):
            gs = gene_stats[gene]
            writer.writerow(
                {
                    "gene": gene,
                    "logic_type": BDQ_GENES.get(gene, ""),
                    "total_snps": gs["total_snps"],
                    "isolates_affected": len(gs["isolate_set"]),
                    "r_calls": gs["r_count"],
                    "s_calls": gs["s_count"],
                    "indet_calls": gs["indet_count"],
                }
            )
    print(f"Gene summary: {gene_summary_file}")
    print("\n===== AMR-Hunter 本地菌株判定摘要 =====")
    print(f"总菌株数: {len(all_isolates)}")
    print(f"BDQ 范围 SNP 位点: {len(bdq_snps)}")
    print(f"有 BDQ 变异的菌株: {len(isolate_variants)}")
    print(f"\n菌株级聚合判定分布:")
    for label in ["R", "INDETERMINATE", "S", "NO_VARIANT"]:
        count = label_counts.get(label, 0)
        print(f"  {label}: {count}")
    r_isolates = [
        iso_id
        for iso_id in isolate_variants
        if any((v["resistance_label"] == "RESISTANT" for v in isolate_variants[iso_id]))
    ]
    print(f"\n高置信耐药菌株（建议优先 MIC 验证）: {len(r_isolates)} 株")
    if r_isolates:
        print("Top 10 耐药菌株:")
        for iso_id in r_isolates[:10]:
            r_vars = [
                v
                for v in isolate_variants[iso_id]
                if v["resistance_label"] == "RESISTANT"
            ]
            genes = sorted(set((v["gene"] for v in r_vars)))
            aas = "; ".join((f"{v['gene']}:{v['aa_change']}" for v in r_vars))
            print(f"  {iso_id}: {aas} [{', '.join(genes)}]")
    print(f"\n输出目录: {output_dir}")


if __name__ == "__main__":
    main()
