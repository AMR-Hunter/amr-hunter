from __future__ import annotations
import csv
from collections import Counter, defaultdict
from pathlib import Path

BASE = Path("c:\\Users\\zhang\\Desktop\\深度学习\\3、文章\\内部数据")
PER_ISOLATE = BASE / "amr_hunter_local_analysis" / "per_isolate_variants.csv"
SERVER = BASE / "final_calls_server.csv"
OUT_VARIANT = BASE / "final_calls_v3_five_category.csv"
OUT_ISOLATE = BASE / "isolate_final_calls_v3.csv"
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


def normalize(call: str) -> str:
    call = call.strip().upper()
    return LABEL_MAP.get(call, call)


def re_resolve(gene: str, call: str) -> str:
    call = normalize(call)
    if gene in STRUCTURAL_GENES and call == "R":
        return "S"
    if gene in EXPLORATORY_GENES and call == "R":
        return "INDET"
    return call


def main() -> None:
    server_rows = {}
    with SERVER.open(encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            key = (row["gene"], row["pos"], row["ref"], row["alt"])
            server_rows[key] = row
    variant_calls: dict[tuple, dict] = {}
    for key, row in server_rows.items():
        gene = row["gene"]
        if gene not in TARGET_GENES:
            continue
        new_call = re_resolve(gene, row["final_call"])
        variant_calls[key] = {
            "gene": gene,
            "pos": row["pos"],
            "ref": row["ref"],
            "alt": row["alt"],
            "region": row["region"],
            "old_call": row["final_call"],
            "status": row["status"],
            "old_pheno": row["pheno"],
            "new_call": new_call,
        }
    isolate_rows = defaultdict(list)
    with PER_ISOLATE.open(encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            gene = row["gene"]
            if gene not in TARGET_GENES:
                continue
            key = (gene, row["pos"], row["ref"], row["alt"])
            if key in variant_calls:
                call = variant_calls[key]["new_call"]
            else:
                call = re_resolve(gene, row["resistance_label"])
            isolate_rows[row["isolate_id"]].append(
                (
                    gene,
                    row["aa_change"],
                    row["mutation_type"],
                    call,
                    row["pos"],
                    row["ref"],
                    row["alt"],
                )
            )
    with OUT_VARIANT.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "gene",
                "pos",
                "ref",
                "alt",
                "region",
                "old_call",
                "status",
                "old_pheno",
                "new_call",
            ],
        )
        writer.writeheader()
        for row in sorted(
            variant_calls.values(), key=lambda r: (r["gene"], int(r["pos"]))
        ):
            writer.writerow(row)
    isolate_calls = {}
    for iso, vars_ in isolate_rows.items():
        calls = [v[3] for v in vars_]
        if "R" in calls:
            isolate_calls[iso] = ("R", sum((1 for c in calls if c == "R")))
        elif "S" in calls:
            isolate_calls[iso] = ("S", 0)
        else:
            isolate_calls[iso] = ("INDET", 0)
    with OUT_ISOLATE.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["isolate_id", "final_call", "n_R"])
        writer.writeheader()
        for iso in sorted(isolate_calls):
            call, n_r = isolate_calls[iso]
            writer.writerow({"isolate_id": iso, "final_call": call, "n_R": n_r})
    total_carriers = len(isolate_calls)
    iso_counts = Counter((c for c, _ in isolate_calls.values()))
    print(f"变异携带株（目标基因）: {total_carriers}")
    print(f"菌株级判定: {dict(iso_counts)}")
    print(
        f"R 株中 n_R 分布: {dict(Counter((n for c, n in isolate_calls.values() if c == 'R')))}"
    )
    vc = Counter((v["new_call"] for v in variant_calls.values()))
    print(f"服务器判定变异（目标基因）: {len(variant_calls)}，重判定: {dict(vc)}")
    per_gene = defaultdict(Counter)
    for v in variant_calls.values():
        if v["new_call"] != "R":
            continue
        gene = v["gene"]
        if v["status"] == "LETHAL_SKIP":
            per_gene[gene]["lethal"] += 1
        elif v["region"] == "REG":
            per_gene[gene]["regulatory"] += 1
        elif v["status"] in {"LOSS_OF_FUNCTION"}:
            per_gene[gene]["lof"] += 1
        else:
            per_gene[gene]["other"] += 1
    print("R 变异按基因与证据:")
    for g in sorted(per_gene, key=lambda g: -sum(per_gene[g].values())):
        print(f"  {g}: {dict(per_gene[g])} 合计 {sum(per_gene[g].values())}")
    fs_by_gene = Counter()
    for iso, vars_ in isolate_rows.items():
        for gene, aa, mtype, call, *_ in vars_:
            if call == "R" and ("fs" in (aa or "") or "frameshift" in (mtype or "")):
                fs_by_gene[gene] += 1
    print("R 且带移码的菌株-变异对:", dict(fs_by_gene))
    pair_by_gene = Counter()
    r_iso_by_gene = defaultdict(set)
    for iso, vars_ in isolate_rows.items():
        for gene, aa, mtype, call, *_ in vars_:
            if call == "R":
                pair_by_gene[gene] += 1
                r_iso_by_gene[gene].add(iso)
    print("R 菌株-变异对按基因:", dict(pair_by_gene))
    print(
        "R 菌株按基因（主要耐药变异）:",
        {
            g: len(s)
            for g, s in sorted(r_iso_by_gene.items(), key=lambda kv: -len(kv[1]))
        },
    )
    screening_r = Counter()
    with PER_ISOLATE.open(encoding="utf-8-sig") as f:
        for row in csv.DictReader(f):
            gene = row["gene"]
            if gene not in TARGET_GENES:
                continue
            key = (gene, row["pos"], row["ref"], row["alt"])
            if key in variant_calls:
                continue
            if re_resolve(gene, row["resistance_label"]) == "R":
                screening_r[gene] += 1
    print("筛选阶段直接判 R（未上服务器）:", dict(screening_r))


if __name__ == "__main__":
    main()
