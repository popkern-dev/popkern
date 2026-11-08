use crate::error::ValidationError;
use crate::topology::BRep;
use crate::topology::radial_order::book_test_helpers::{DEG, book, tol};

/// A ring linked in radial order, and every ring of two, passes.
#[test]
fn rings_in_radial_order_pass() {
    let (mut pages, spine, _) = book(&[0.0, 200.0 * DEG, 90.0 * DEG], &[true, false, true]);
    pages.order_ring(spine, &tol()).expect("ordered");
    assert!(pages.validate_ring_order_all().is_empty());
    assert!(
        BRep::unit_cube(tol())
            .expect("cube")
            .validate_ring_order_all()
            .is_empty()
    );
}

/// Negative control: the same pages linked out of order are reported
/// with both orders.
#[test]
fn a_ring_out_of_order_is_reported() {
    let (pages, spine, fins) = book(&[0.0, 200.0 * DEG, 90.0 * DEG], &[true, false, true]);
    let errors = pages.validate_ring_order_all();
    assert!(
        matches!(errors.as_slice(), [ValidationError::RingOrder { edge, stored, radial }]
            if *edge == spine
                && *stored == fins
                && *radial == vec![fins[0], fins[2], fins[1]]),
        "{errors:?}"
    );
}

/// Negative control: coincident pages leave the order undecided.
#[test]
fn an_undecided_ring_is_reported() {
    let (pages, spine, _) = book(&[0.0, 90.0 * DEG, 90.0 * DEG], &[true, true, true]);
    let errors = pages.validate_ring_order_all();
    assert!(
        matches!(errors.as_slice(), [ValidationError::RingOrderUndecided { edge, .. }] if *edge == spine),
        "{errors:?}"
    );
}
