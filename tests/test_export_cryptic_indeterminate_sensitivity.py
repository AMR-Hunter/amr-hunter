from pathlib import Path
import sys

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))
from scripts.export_cryptic_indeterminate_sensitivity import (
    SCENARIOS,
    summarize_predictions,
)


def scenario_by_name(name: str) -> dict[str, object]:
    return next((scenario for scenario in SCENARIOS if scenario["scenario"] == name))


def test_indeterminate_scenarios_reclassify_only_triggered_sensitive_calls():
    isolates = [
        {
            "unique_id": "fn-core",
            "ena_sample": "S1",
            "bdq_binary_phenotype": "R",
            "bdq_mic": "0.5",
            "amr_hunter_isolate_prediction": "S",
            "call_source_summary": "direct_record_match:1|unsupported_variant:1",
        },
        {
            "unique_id": "fn-mtrb",
            "ena_sample": "S2",
            "bdq_binary_phenotype": "R",
            "bdq_mic": ">2",
            "amr_hunter_isolate_prediction": "S",
            "call_source_summary": "direct_record_match:1|unsupported_variant:1",
        },
        {
            "unique_id": "tp",
            "ena_sample": "S3",
            "bdq_binary_phenotype": "R",
            "bdq_mic": "1.0",
            "amr_hunter_isolate_prediction": "R",
            "call_source_summary": "direct_record_match:1",
        },
        {
            "unique_id": "tn",
            "ena_sample": "S4",
            "bdq_binary_phenotype": "S",
            "bdq_mic": "0.06",
            "amr_hunter_isolate_prediction": "S",
            "call_source_summary": "direct_record_match:1",
        },
    ]
    variants_by_isolate = {
        "fn-core": [
            {
                "gene": "Rv0678",
                "mutation": "p.Arg90Cys",
                "call_source": "unsupported_variant",
            }
        ],
        "fn-mtrb": [
            {
                "gene": "mtrB",
                "mutation": "p.Met517Leu",
                "call_source": "unsupported_variant",
            }
        ],
        "tn": [
            {
                "gene": "mmpL5",
                "mutation": "p.Ile948Val",
                "call_source": "direct_record_match",
            }
        ],
    }
    baseline, baseline_changed = summarize_predictions(
        isolates, variants_by_isolate, scenario_by_name("baseline")
    )
    assert baseline["tp"] == "1"
    assert baseline["tn"] == "1"
    assert baseline["fn"] == "2"
    assert baseline["coverage"] == "1.0000"
    assert baseline_changed == []
    core, core_changed = summarize_predictions(
        isolates, variants_by_isolate, scenario_by_name("core_unsupported")
    )
    assert core["fn"] == "1"
    assert core["r_indeterminate"] == "1"
    assert core["pred_INDETERMINATE"] == "1"
    assert core["coverage"] == "0.7500"
    assert core_changed[0]["trigger_variants"] == "Rv0678:p.Arg90Cys"
    mtrb, mtrb_changed = summarize_predictions(
        isolates, variants_by_isolate, scenario_by_name("mtrB_context")
    )
    assert mtrb["fn"] == "1"
    assert mtrb["r_indeterminate"] == "1"
    assert mtrb["pred_INDETERMINATE"] == "1"
    assert mtrb_changed[0]["trigger_variants"] == "mtrB:p.Met517Leu"


if __name__ == "__main__":
    test_indeterminate_scenarios_reclassify_only_triggered_sensitive_calls()
    print("cryptic indeterminate sensitivity tests passed")
