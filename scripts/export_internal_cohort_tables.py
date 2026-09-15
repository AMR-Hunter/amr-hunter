from __future__ import annotations
import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.stdout.reconfigure(encoding="utf-8", errors="replace")
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
LOF_R_GENES = {"Rv0678", "pepQ"}
STRUCTURAL_GENES = {"mtrA", "mtrB", "mmpS5", "mmpL5", "lpqB"}
EXPLORATORY_GENES = {"glpK", "Rv1979c"}
DEFAULT_ANALYSIS_DIR = Path(
    "C:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\内部数据\\outdir\\analysis"
)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, fieldnames: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(fieldnames)
        writer.writerows(rows)


def server_direction(gene: str, db_row: dict[str, str] | None) -> str:
    if db_row is None:
        return "INDET"
    if (db_row["status"] or "").strip() == "BOLTZ_READY":
        return "S" if gene in STRUCTURAL_GENES else "INDET"
    phenotype = (db_row["resistance_phenotype"] or "").strip().upper()
    if phenotype == "RESISTANT":
        return "R"
    if phenotype == "SENSITIVE":
        return "S"
    return "INDET"


def classify_variant(
    gene: str, variant_type: str, db_row: dict[str, str] | None
) -> str:
    if variant_type.startswith("Intergenic"):
        return "INDET"
    if variant_type.startswith("Synonymous"):
        return "S"
    is_coding_consequence = (
        variant_type.startswith("Nonsynonymous") or "frameshift" in variant_type
    )
    if not is_coding_consequence:
        return "INDET"
    if gene in STRUCTURAL_GENES:
        return "S"
    if gene in EXPLORATORY_GENES:
        return "INDET"
    if "frameshift" in variant_type:
        return "R" if gene in LOF_R_GENES else "INDET"
    raw = server_direction(gene, db_row)
    if gene in STRUCTURAL_GENES and raw == "R":
        return "S"
    if gene in EXPLORATORY_GENES and raw == "R":
        return "INDET"
    return raw


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--analysis-dir", type=Path, default=DEFAULT_ANALYSIS_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_ANALYSIS_DIR)
    parser.add_argument("--variant-output", default="variant_calls_final.csv")
    parser.add_argument("--isolate-output", default="isolate_calls_final.csv")
    args = parser.parse_args()
    per_sample_rows = read_csv(args.analysis_dir / "bdq_variants_per_sample.csv")
    db_rows = read_csv(args.analysis_dir / "db_match.csv")
    db = {(r["gene"], r["pos"], r["ref"], r["alt"]): r for r in db_rows}
    rows = [
        r
        for r in per_sample_rows
        if r["if_phylogeny"] == "1" and r["gene"] in TARGET_GENES
    ]
    dropped = Counter(
        (r["gene"] for r in per_sample_rows if r["gene"] not in TARGET_GENES)
    )
    print(f"主队列样本-变异对: {len(rows)} | 非目标基因已剔除: {dict(dropped)}")
    variants: dict[tuple[str, str, str, str], str] = {}
    for r in rows:
        key = (r["gene"], r["pos"], r["ref"], r["alt"])
        variants.setdefault(key, r["variant_type"])
    calls = {
        key: classify_variant(key[0], vtype, db.get(key))
        for key, vtype in variants.items()
    }
    variant_counts = Counter(calls.values())
    print(f"唯一变异: {len(calls)} | 判定分布: {dict(variant_counts)}")
    print(
        f"  耐药判定按基因: {dict(Counter((k[0] for k, c in calls.items() if c == 'R')))}"
    )
    write_csv(
        args.output_dir / args.variant_output,
        ["gene", "pos", "ref", "alt", "variant_type", "call"],
        [[k[0], k[1], k[2], k[3], variants[k], calls[k]] for k in sorted(variants)],
    )
    per_sample: dict[str, list[str]] = {}
    for r in rows:
        key = (r["gene"], r["pos"], r["ref"], r["alt"])
        per_sample.setdefault(r["sample_id"], []).append(calls[key])
    isolate_calls: dict[str, str] = {}
    for sample_id, cs in per_sample.items():
        if "R" in cs:
            isolate_calls[sample_id] = "R"
        elif "S" in cs:
            isolate_calls[sample_id] = "S"
        else:
            isolate_calls[sample_id] = "INDET"
    isolate_counts = Counter(isolate_calls.values())
    print(f"菌株数: {len(isolate_calls)} | 判定分布: {dict(isolate_counts)}")
    carriers_of_indet = sum(
        (1 for cs in per_sample.values() if "INDET" in cs and "R" not in cs)
    )
    print(f"  携带不确定变异但无耐药判定的菌株: {carriers_of_indet}")
    write_csv(
        args.output_dir / args.isolate_output,
        ["sample_id", "call"],
        [[sid, isolate_calls[sid]] for sid in sorted(isolate_calls)],
    )
    print(f"written: {args.output_dir / args.variant_output}")
    print(f"written: {args.output_dir / args.isolate_output}")


if __name__ == "__main__":
    main()
