use super::*;

/// The table hands out exactly the products it was built from, so
/// building it once changes no result.
#[test]
fn the_table_equals_the_products_bit_for_bit() {
    let bits = |p: &Poly27| p.map(f64::to_bits);
    for a in 0..8 {
        assert_eq!(bits(&shape(a)), bits(&shape_product(a)), "corner {a}");
        for axis in 0..3 {
            assert_eq!(
                bits(&shape_derivative(a, axis)),
                bits(&derivative_product(a, axis)),
                "corner {a} axis {axis}"
            );
        }
    }
    assert_eq!(shape(8), [0.0; 27]);
    assert_eq!(shape_derivative(0, 3), [0.0; 27]);
}
