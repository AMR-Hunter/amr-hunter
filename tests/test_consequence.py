from core.consequence import (
    is_loss_of_function_consequence,
    resolve_consequence_level_decision,
)


def test_loss_of_function_consequence_detection() -> None:
    assert is_loss_of_function_consequence("p.Asp47fs")
    assert is_loss_of_function_consequence("p.Met1?")
    assert is_loss_of_function_consequence("LoF")
    assert is_loss_of_function_consequence("start_lost")
    assert is_loss_of_function_consequence("p.Arg132*")
    assert not is_loss_of_function_consequence("p.Gly121Arg")


def test_negative_logic_lof_maps_to_resistant() -> None:
    decision = resolve_consequence_level_decision(
        {"logic_type": "NEGATIVE", "mutation": "p.Asp47fs"}
    )
    assert decision is not None
    assert decision.resistance_phenotype == "RESISTANT"
    assert decision.functional_state == "LOSS_OF_FUNCTION"
    assert decision.evidence_code == "NEGATIVE_LOF_SIGNATURE"


def test_structural_and_positive_logic_lof_map_to_sensitive() -> None:
    structural = resolve_consequence_level_decision(
        {"logic_type": "STRUCTURAL", "mutation_class": "frameshift"}
    )
    positive = resolve_consequence_level_decision(
        {"logic_type": "POSITIVE", "mutation": "p.Met1?"}
    )
    assert structural is not None
    assert structural.resistance_phenotype == "SENSITIVE"
    assert structural.evidence_code == "STRUCTURAL_LOF_SIGNATURE"
    assert positive is not None
    assert positive.resistance_phenotype == "SENSITIVE"
    assert positive.evidence_code == "POSITIVE_TARGET_COLLAPSE"
