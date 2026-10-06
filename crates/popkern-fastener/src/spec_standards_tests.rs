use crate::spec::ThreadSpec;
use crate::standards::{self, ClearanceFit};

/// The same M6 at millimetre magnitudes and in SI metres.
fn both() -> (ThreadSpec, ThreadSpec) {
    let mm = match ThreadSpec::parse("M6") {
        Ok(s) => s,
        Err(e) => panic!("{e}"),
    };
    let si =
        ThreadSpec::custom(mm.major_diameter() * 1e-3, mm.pitch() * 1e-3).with_millimetre(1e-3);
    (mm, si)
}

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 * a.abs().max(b.abs())
}

#[test]
fn an_si_m6_gets_the_iso_head_not_derived_proportions() {
    let (mm, si) = both();
    // The ISO 4017 row (s = 10, k = 4), not 1.6 d = 9.6.
    assert_eq!(mm.iso_hex_head().across_flats, 10.0);
    assert!(close(si.iso_hex_head().across_flats, 0.010));
    assert!(close(
        si.iso_hex_head().height,
        mm.iso_hex_head().height * 1e-3
    ));
    assert!(close(
        si.iso_socket_head().diameter,
        mm.iso_socket_head().diameter * 1e-3
    ));
    assert!(close(
        si.iso_countersunk_head_diameter(),
        mm.iso_countersunk_head_diameter() * 1e-3
    ));
    assert!(close(
        si.iso_hex_nut_thickness(),
        mm.iso_hex_nut_thickness() * 1e-3
    ));
    for fit in [
        ClearanceFit::Close,
        ClearanceFit::Normal,
        ClearanceFit::Loose,
    ] {
        assert!(close(
            si.iso_clearance_diameter(fit),
            mm.iso_clearance_diameter(fit) * 1e-3
        ));
    }
    // Without the scale, the SI major misses every row.
    assert_eq!(standards::hex_head(0.006).across_flats, 1.6 * 0.006);
}

#[test]
fn a_custom_clearance_diameter_is_already_in_the_spec_unit() {
    let (_, si) = both();
    assert_eq!(
        si.iso_clearance_diameter(ClearanceFit::Custom(0.0065)),
        0.0065
    );
}
