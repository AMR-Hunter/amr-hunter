from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_fn_taxonomy import (
    build_category_rows,
    build_fn_rows,
    build_mic_rows,
)


def test_fn_taxonomy_and_mic_stratification():
    isolates = [
        {
            "unique_id": "fn-core",
            "ena_sample": "S1",
            "bdq_binary_phenotype": "R",
            "bdq_mic": "0.5",
            "amr_hunter_isolate_prediction": "S",
            "correct": "false",
        },
        {
            "unique_id": "fn-mtrb-double",
            "ena_sample": "S2",
            "bdq_binary_phenotype": "R",
            "bdq_mic": ">2",
            "amr_hunter_isolate_prediction": "S",
            "correct": "false",
        },
        {
            "unique_id": "tp",
            "ena_sample": "S3",
            "bdq_binary_phenotype": "R",
            "bdq_mic": "1.0",
            "amr_hunter_isolate_prediction": "R",
            "correct": "true",
        },
        {
            "unique_id": "tn",
            "ena_sample": "S4",
            "bdq_binary_phenotype": "S",
            "bdq_mic": "0.06",
            "amr_hunter_isolate_prediction": "S",
            "correct": "true",
        },
    ]
    variants_by_isolate = {
        "fn-core": [
            {
                "gene": "mmpL5",
                "mutation": "p.Ile948Val",
                "amr_hunter_prediction": "S",
                "call_source": "direct_record_match",
            },
            {
                "gene": "Rv0678",
                "mutation": "p.Arg90Cys",
                "amr_hunter_prediction": "UNKNOWN",
                "call_source": "unsupported_variant",
            },
        ],
        "fn-mtrb-double": [
            {
                "gene": "mtrB",
                "mutation": "p.Met517Leu",
                "amr_hunter_prediction": "UNKNOWN",
                "call_source": "unsupported_variant",
            },
            {
                "gene": "mtrB",
                "mutation": "p.Pro18Ser",
                "amr_hunter_prediction": "UNKNOWN",
                "call_source": "unsupported_variant",
            },
        ],
    }
    fn_rows = build_fn_rows(isolates, variants_by_isolate)
    by_id = {row["unique_id"]: row for row in fn_rows}
    assert by_id["fn-core"]["taxonomy"] == "core_gene_unsupported"
    assert by_id["fn-core"]["core_unsupported"] == "Rv0678:p.Arg90Cys"
    assert by_id["fn-mtrb-double"]["taxonomy"] == "mtrB_double_hit_context"
    assert by_id["fn-mtrb-double"]["mic_bin"] == ">2"
    category_rows = build_category_rows(fn_rows, isolates)
    category_counts = {row["taxonomy"]: row["fn_isolates"] for row in category_rows}
    assert category_counts == {
        "core_gene_unsupported": "1",
        "mtrB_double_hit_context": "1",
    }
    mic_rows = build_mic_rows(isolates)
    by_mic = {(row["phenotype"], row["mic_bin"]): row for row in mic_rows}
    assert by_mic["R", "0.5"]["primary_metric_value"] == "0.0000"
    assert by_mic["R", "1.0"]["primary_metric_value"] == "1.0000"
    assert by_mic["R", ">2"]["primary_metric_value"] == "0.0000"
    assert by_mic["S", "<=0.25"]["primary_metric_value"] == "1.0000"


if __name__ == "__main__":
    test_fn_taxonomy_and_mic_stratification()
    print("cryptic FN taxonomy tests passed")
