//! Tests for the direct-construction result / diagnostic types.

use super::*;

#[test]
fn validation_level_is_copyable_and_comparable() {
    let a = ValidationLevel::Standard;
    let b = a;
    assert_eq!(a, b);
    assert_ne!(a, ValidationLevel::None);
    assert_ne!(a, ValidationLevel::Full);
}

#[test]
fn corner_kind_distinct() {
    assert_ne!(CornerKind::Central, CornerKind::OuterHex);
}

#[test]
fn formula_family_tags_are_distinct() {
    let families = [
        FormulaFamily::ChiyokuraKimuraStandard,
        FormulaFamily::ChiyokuraKimuraTangentCorner,
        FormulaFamily::KatoKimuraChiyokura1991,
        FormulaFamily::ShirmanSequinEnclosure,
        FormulaFamily::PetersTwoStripeRow1,
    ];
    // All distinct pairwise.
    for (i, a) in families.iter().enumerate() {
        for b in &families[i + 1..] {
            assert_ne!(a, b, "family duplicate");
        }
    }
}

#[test]
fn feasibility_warning_variants_distinguishable() {
    let w1 = FeasibilityWarning::TangentCornerVariant { corner_idx: 0 };
    let w2 = FeasibilityWarning::ValenceRequiresPolar { n: 5 };
    // Make sure Display / Debug don't panic and produce distinct output.
    assert_ne!(format!("{w1:?}"), format!("{w2:?}"));
}
