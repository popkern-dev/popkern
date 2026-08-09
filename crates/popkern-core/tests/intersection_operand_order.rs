//! Operand-order equivalence of the surface-surface intersection
//! solver on the measured grazing fillet pair.
//!
//! The fixture is the cylinder strip × NURBS cap sector from a served
//! mixed-radius pyramid corner (captured 2026-08-08), the pair whose
//! crossing the solver historically found cylinder-first and lost
//! free-form-first: the seed Newton anchored its advancing plane on
//! operand 1 only, tangent-classified seeds were dropped when no
//! coincident region covered them, and a seed converging slightly
//! outside a non-periodic window terminated the march at its first
//! step. All three are fixed (dual-anchor seeding, angular-floor
//! tangent demotion, non-periodic march-window entry band); this pin
//! holds the solver to order-equivalent results on the real geometry.

use popkern_core::geometry::intersection::degeneracy::DefaultDegeneracyHandler;
use popkern_core::geometry::intersection::{
    SurfaceSurfaceIntersection, intersect_surfaces_with_handler,
};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;

struct TrimmedSurface {
    geometry: SurfaceGeometry,
    u_range: (f64, f64),
    v_range: (f64, f64),
}

fn load_pair() -> (TrimmedSurface, TrimmedSurface) {
    let raw = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/fixtures/grazing_strip_cap_pair.json"
    ))
    .expect("read grazing pair fixture");
    let v: serde_json::Value = serde_json::from_str(&raw).expect("parse fixture");
    let cyl: SurfaceGeometry = serde_json::from_value(v["cylinder"].clone()).expect("cylinder");
    let ff: SurfaceGeometry = serde_json::from_value(v["free_form"].clone()).expect("free_form");
    let range = |key: &str| -> (f64, f64) {
        (
            v[key][0].as_f64().expect("range lo"),
            v[key][1].as_f64().expect("range hi"),
        )
    };
    (
        TrimmedSurface {
            geometry: cyl,
            u_range: range("cylinder_u"),
            v_range: range("cylinder_v"),
        },
        TrimmedSurface {
            geometry: ff,
            u_range: range("free_form_u"),
            v_range: range("free_form_v"),
        },
    )
}

#[test]
fn grazing_strip_cap_pair_intersects_equivalently_in_both_operand_orders() {
    let (cyl, ff) = load_pair();
    let (cu, cv) = (cyl.u_range, cyl.v_range);
    let (fu, fv) = (ff.u_range, ff.v_range);
    let (cyl, ff) = (cyl.geometry, ff.geometry);
    // The self-intersection validator's coarse tolerance: 2× linear.
    let tol = Tolerance::default();
    let coarse = Tolerance::default()
        .with_linear(tol.linear() * 2.0)
        .with_angular(tol.angular())
        .with_parametric(tol.parametric_base());

    let analytic_first = intersect_surfaces_with_handler(
        &cyl,
        cu,
        cv,
        &ff,
        fu,
        fv,
        &coarse,
        &DefaultDegeneracyHandler,
    )
    .expect("analytic-first intersect");
    let free_form_first = intersect_surfaces_with_handler(
        &ff,
        fu,
        fv,
        &cyl,
        cu,
        cv,
        &coarse,
        &DefaultDegeneracyHandler,
    )
    .expect("free-form-first intersect");

    let count = |r: &popkern_core::geometry::intersection::SurfaceSurfaceIntersectionResult| {
        let curves = r
            .intersections
            .iter()
            .filter(|i| matches!(i, SurfaceSurfaceIntersection::Curve(_)))
            .count();
        (curves, r.intersections.len() - curves)
    };
    let (af_curves, af_regions) = count(&analytic_first);
    let (ff_curves, ff_regions) = count(&free_form_first);

    assert!(
        af_curves > 0,
        "analytic-first must find the strip×cap contact"
    );
    assert!(
        ff_curves > 0,
        "free-form-first must find the strip×cap contact"
    );
    assert_eq!(
        (af_curves, af_regions),
        (ff_curves, ff_regions),
        "operand order must not change the intersection set"
    );
}
