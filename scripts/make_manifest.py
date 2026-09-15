import csv

rows = list(
    csv.DictReader(
        open(
            "d:\\python_work\\amr_hunter\\amr-hunter\\reports\\publication\\external_bdq_phenotype_validation_isolates.csv",
            encoding="utf-8-sig",
        )
    )
)
print("isolates:", len(rows))
with open(
    "C:\\Users\\zhang\\Desktop\\cryptic_manifest.tsv", "w", newline="", encoding="utf-8"
) as f:
    w = csv.writer(f, delimiter="\t")
    w.writerow(["unique_id", "vcf_path"])
    for r in rows:
        p = r["vcf_path"]
        if p.startswith("../"):
            p = p[3:]
        w.writerow([r["unique_id"], p])
print("manifest written")
