from core.analyzer import resolve_analysis_decision


def _decision(mutation: dict) -> tuple:
    decision = resolve_analysis_decision(mutation)
    return (
        decision.resistance_phenotype,
        decision.functional_state,
        decision.evidence_code,
    )


def test_exploratory_gene_lof_is_unknown():
    for gene in ("Rv1979c", "glpK"):
        assert _decision(
            {
                "gene_name": gene,
                "logic_type": "NEGATIVE",
                "status": "LOSS_OF_FUNCTION",
                "region_type": "CDS",
                "aa_change": "100*",
            }
        ) == ("UNKNOWN", "LOSS_OF_FUNCTION", "EXPLORATORY_ONLY")


def test_exploratory_gene_completed_is_unknown():
    for gene in ("Rv1979c", "glpK"):
        assert _decision(
            {
                "gene_name": gene,
                "logic_type": "NEGATIVE",
                "status": "COMPLETED",
                "region_type": "CDS",
                "aa_change": "p.Ala100Val",
            }
        ) == ("UNKNOWN", "INTACT", "EXPLORATORY_ONLY")


def test_exploratory_gene_overrides_stored_resistant():
    assert _decision(
        {
            "gene_name": "glpK",
            "logic_type": "NEGATIVE",
            "status": "LOSS_OF_FUNCTION",
            "region_type": "CDS",
            "aa_change": "100*",
            "resistance_phenotype": "RESISTANT",
            "functional_state": "LOSS_OF_FUNCTION",
            "evidence_code": "NEGATIVE_LOF_SIGNATURE",
        }
    ) == ("UNKNOWN", "LOSS_OF_FUNCTION", "EXPLORATORY_ONLY")


def test_gene_role_field_triggers_exploratory_short_circuit():
    assert _decision(
        {
            "gene_name": "custom_gene",
            "gene_role": "exploratory",
            "logic_type": "NEGATIVE",
            "status": "COMPLETED",
            "region_type": "CDS",
            "aa_change": "p.Ala100Val",
        }
    ) == ("UNKNOWN", "INTACT", "EXPLORATORY_ONLY")


def test_mtra_lof_now_sensitive():
    assert _decision(
        {
            "gene_name": "mtrA",
            "logic_type": "STRUCTURAL",
            "status": "LOSS_OF_FUNCTION",
            "region_type": "CDS",
            "aa_change": "229*",
        }
    ) == ("SENSITIVE", "LOSS_OF_FUNCTION", "STRUCTURAL_LOF_SIGNATURE")


def test_mtrb_stored_negative_decision_refreshed_under_structural():
    assert _decision(
        {
            "gene_name": "mtrB",
            "logic_type": "STRUCTURAL",
            "status": "LOSS_OF_FUNCTION",
            "region_type": "CDS",
            "aa_change": "310*",
            "resistance_phenotype": "RESISTANT",
            "functional_state": "LOSS_OF_FUNCTION",
            "evidence_code": "NEGATIVE_LOF_SIGNATURE",
        }
    ) == ("SENSITIVE", "LOSS_OF_FUNCTION", "STRUCTURAL_LOF_SIGNATURE")


def test_mtrb_synonymous_sensitive():
    assert _decision(
        {
            "gene_name": "mtrB",
            "logic_type": "STRUCTURAL",
            "status": "COMPLETED",
            "region_type": "CDS",
            "aa_change": "p.Met517Met",
        }
    ) == ("SENSITIVE", "INTACT", "SYNONYMOUS_NO_PROTEIN_CHANGE")


def test_mtrb_deregulation_tag_no_longer_resistant():
    assert _decision(
        {
            "gene_name": "mtrA",
            "logic_type": "STRUCTURAL",
            "status": "COMPLETED",
            "region_type": "CDS",
            "aa_change": "p.Ala100Val",
            "synergy_tags": "MTRAB_DEREGULATION|SYNERGIC_RESISTANCE",
        }
    ) == ("SENSITIVE", "INTACT", "STRUCTURAL_FUNCTION_RETAINED")


def test_rv0678_lof_still_resistant():
    assert _decision(
        {
            "gene_name": "Rv0678",
            "logic_type": "NEGATIVE",
            "status": "LOSS_OF_FUNCTION",
            "region_type": "CDS",
            "aa_change": "R132*",
        }
    ) == ("RESISTANT", "LOSS_OF_FUNCTION", "NEGATIVE_LOF_SIGNATURE")
