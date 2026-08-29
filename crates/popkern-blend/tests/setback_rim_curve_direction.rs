//! PK-303: setback rim/seam edge-curve direction and served pcurve
//! coherence on the 12-edge G2 cube.
//!
//! Root cause pinned here: `refit_cap_rim` authored each rim edge's
//! 3D curve along the CAP-side half-edge, but the kernel convention
//! puts the CANONICAL half-edge's origin at `t_start` — whenever the
//! cap side was the twin, the curve was stored reversed. Geometry was
//! exact (curve on surface at 1e-15, pcurve images exact), so no
//! quality gate caught it, but every fraction-sampling consumer (the
//! pcurve deviation metric, mesher boundary conformity, STEP trace
//! sync, declared accuracy) read the curve backwards: 300 of 528
//! half-edges measured "stale" at exactly the rim chord length
//! (4.96e-2), and both whole-model `fix_stale_pcurves` passes
//! re-cleared and re-fit them every run without ever repairing them.

use popkern_blend::{ContinuityLevel, FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::curve::CurveEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn served_g2_cube() -> BRep {
    let tol = Tolerance::default();
    let cube = BRep::unit_cube(tol).expect("unit cube builds");
    let specs: Vec<FilletSpec> = cube
        .edges
        .keys()
        .take(12)
        .map(|e| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.1),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G2,
        ..FilletOptions::default()
    };
    fillet_edges(&cube, &specs, &options, tol, None)
        .expect("12-edge G2 cube fillet serves")
        .brep
}

/// Every edge curve on the served body runs in the canonical
/// direction: the canonical half-edge's origin sits at `t_start`.
/// (An edge whose endpoints coincide within tolerance — closed or
/// near-closed — carries no direction evidence and is skipped.)
#[test]
fn served_g2_cube_edge_curves_run_canonical_forward() {
    let brep = served_g2_cube();
    let mut checked = 0usize;
    for (eh, edge) in brep.edges.iter() {
        let Some(curve) = edge.curve.as_ref() else {
            continue;
        };
        let he = &brep.half_edges[edge.half_edge];
        let origin = brep.vertices[he.origin].position;
        let dest = brep.vertices[brep.half_edges[he.twin].origin].position;
        let separation = origin.distance_to(&dest);
        if separation < 1.0e-6 {
            continue;
        }
        let start = curve.geometry.point_at(curve.t_start);
        let d_origin = start.distance_to(&origin);
        let d_dest = start.distance_to(&dest);
        assert!(
            d_origin < d_dest,
            "edge {eh:?}: curve(t_start) is nearer the canonical DEST \
             ({d_dest:.3e}) than the canonical ORIGIN ({d_origin:.3e}) — \
             the curve is stored reversed against the canonical half-edge \
             convention (PK-303 regression)"
        );
        checked += 1;
    }
    assert!(
        checked > 200,
        "direction check must cover the fan/rim edge population (checked {checked})"
    );
}

/// The served body carries ZERO stale pcurves: with rim curves in
/// canonical direction, the fan population measures inside the
/// staleness gate and the whole-model repair passes have nothing to
/// re-clear. Before the PK-303 root fix this was 300 of 528, every
/// run, unrepairable.
#[test]
fn served_g2_cube_has_no_stale_pcurves() {
    let brep = served_g2_cube();
    let stale = brep.stale_pcurve_handles();
    assert!(
        stale.is_empty(),
        "served 12-edge G2 cube carries {} stale pcurve(s) (first: {:?}) — \
         a fan/rim authoring or repair regression (PK-303)",
        stale.len(),
        stale.first(),
    );
}
