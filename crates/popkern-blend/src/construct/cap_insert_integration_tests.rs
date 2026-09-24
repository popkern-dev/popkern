//! §8.3 end-to-end integration tests for cap insertion on open-shell
//! topology (`BRep::open_box`). Split out of `cap_insert_tests.rs`
//! to stay under the 500-line structural cap.
//!
//! Unit-level tests for `endpoint_is_closed`, `insert_cap_face`, and
//! the closed-body cap-warning paths live in the sibling file.

use super::*;
use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::{CapKind, fillet_edges};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

fn tol() -> Tolerance {
    Tolerance::default()
}

/// Return all interior edges of `open_box` with at least one
/// endpoint on the boundary, sorted by handle for determinism.
fn boundary_adjacent_interior_edges(ob: &BRep) -> Vec<popkern_core::topology::EdgeHandle> {
    let mut edges: Vec<_> = ob
        .edges
        .iter()
        .filter(|(eh, _)| {
            if ob.boundary_edges.contains(eh) {
                return false;
            }
            let he = ob.edges[*eh].half_edge;
            let twin = ob.half_edges[he].twin;
            if twin == popkern_core::topology::HalfEdgeHandle::default() {
                return false;
            }
            let v_start = ob.half_edges[he].origin;
            let v_end = ob.half_edges[twin].origin;
            let touches_boundary = |v: popkern_core::topology::VertexHandle| -> bool {
                ob.half_edges
                    .values()
                    .any(|h| h.origin == v && ob.boundary_edges.contains(&h.edge))
            };
            touches_boundary(v_start) || touches_boundary(v_end)
        })
        .map(|(h, _)| h)
        .collect();
    edges.sort();
    edges
}

/// The edges filleted by the open-box fixtures: the first `n`
/// boundary-adjacent interior edges, at r = 0.1.
fn open_box_specs(ob: &BRep, n: usize) -> Vec<FilletSpec> {
    let edges = boundary_adjacent_interior_edges(ob);
    assert!(edges.len() >= n, "open_box has {} such edges", edges.len());
    edges
        .into_iter()
        .take(n)
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.1),
        })
        .collect()
}

fn fillet_open_box(n_edges: usize, kind: CapKind) -> Result<crate::BlendResult, crate::BlendError> {
    let ob = BRep::open_box(tol()).unwrap();
    let opts = FilletOptions {
        cap_kind: kind,
        ..FilletOptions::default()
    };
    fillet_edges(&ob, &open_box_specs(&ob, n_edges), &opts, tol(), None)
}

/// Open cap on open_box: the strip runs out at the sheet boundary and
/// its free end settles onto the strip's section by the rim plane
/// (PK-530). The rim edges split at the fillet's anchors keep their
/// curves. The whole body holds check 9.
#[test]
fn fillet_open_box_interior_edge_with_open_cap_validates() {
    let result = fillet_open_box(1, CapKind::Open).expect("open-cap fillet serves");
    let cap_warns = result
        .edge_warnings
        .iter()
        .filter(|(_, w)| matches!(w, BlendWarning::CapApplied { .. }))
        .count();
    assert_eq!(
        cap_warns, 0,
        "Open cap must not emit CapApplied warnings on open_box"
    );
    let errors = result.brep.validate_all();
    assert!(errors.is_empty(), "validate_all: {errors:?}");
    assert!(
        result.brep.edges.values().all(|e| e.curve.is_some()),
        "every edge carries a curve"
    );
}

/// Inserted caps on an open sheet author exact boundaries (PK-532,
/// PK-536). The strip | cap edge is the strip's end section. A flat
/// cap keeps its rim as lines in its plane. Rounded and tangent-match
/// rims become two meridians meeting at the outward pole. Single, two-edge and
/// all-boundary-edge fillets serve `validate_all` clean with every
/// cap inserted. Before, the strip | cap edge had no curve and check 9
/// silently skipped it.
#[test]
fn inserted_open_sheet_caps_serve_exact_boundaries() {
    use popkern_core::geometry::surface::SurfaceGeometry;
    let ob = BRep::open_box(tol()).unwrap();
    let all = boundary_adjacent_interior_edges(&ob).len();
    for (kind, n) in [
        (CapKind::Flat, 1),
        (CapKind::Flat, 2),
        (CapKind::Flat, all),
        (CapKind::Rounded, 1),
        (CapKind::Rounded, 2),
        (CapKind::Rounded, all),
        (CapKind::TangentMatch, 1),
        (CapKind::TangentMatch, 2),
        (CapKind::TangentMatch, all),
    ] {
        let result =
            fillet_open_box(n, kind).unwrap_or_else(|e| panic!("{kind:?} x{n}: must serve: {e}"));
        let errors = result.brep.validate_all();
        assert!(errors.is_empty(), "{kind:?} x{n}: validate_all: {errors:?}");
        let inserted = result
            .edge_warnings
            .iter()
            .filter(|(_, w)| {
                matches!(
                    w,
                    BlendWarning::CapApplied {
                        topology_inserted: true,
                        ..
                    }
                )
            })
            .count();
        assert_eq!(inserted, n, "{kind:?} x{n}: every free end capped");
        if kind == CapKind::Rounded {
            let spheres = result
                .brep
                .faces
                .values()
                .filter(|f| {
                    f.surface
                        .as_ref()
                        .is_some_and(|sof| matches!(sof.geometry, SurfaceGeometry::Sphere { .. }))
                })
                .count();
            assert_eq!(spheres, n, "{kind:?} x{n}: one quarter-ball per cap");
        }
    }
}

/// The tangent-match cap on the served open box (PK-536): along the
/// strip | cap section it shares the cylinder strip's tangent plane
/// (G1), and its normal curvature across the section matches the
/// strip's straight ruling (G2, the §7 relative bar 1e-2).
#[test]
fn tangent_match_cap_is_curvature_continuous_across_the_section() {
    use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
    use popkern_core::geometry::{CurveEval, Vec3};
    let result = fillet_open_box(1, CapKind::TangentMatch).expect("serves");
    let brep = &result.brep;
    let face_surface = |f| brep.faces.get(f).and_then(|x| x.surface.as_ref());
    let (cap, cylinder, section) = brep
        .edges
        .values()
        .find_map(|e| {
            let he = brep.half_edges.get(e.half_edge)?;
            let faces = [
                brep.loops.get(he.loop_)?.face,
                brep.loops.get(brep.half_edges.get(he.twin)?.loop_)?.face,
            ];
            let (mut cap, mut cyl) = (None, None);
            for f in faces {
                match &face_surface(f)?.geometry {
                    SurfaceGeometry::Nurbs(n) => cap = Some(n.clone()),
                    SurfaceGeometry::Cylinder {
                        origin,
                        axis,
                        radius,
                    } => cyl = Some((*origin, *axis, *radius)),
                    _ => {}
                }
            }
            Some((cap?, cyl?, e.curve.clone()?))
        })
        .expect("a strip | cap section edge");
    let (origin, axis, r) = cylinder;
    let axis = axis.normalized().expect("axis");
    let (mut worst_g1, mut worst_g2) = (0.0_f64, 0.0_f64);
    for k in 0..=16 {
        let t = section.t_start + (section.t_end - section.t_start) * k as f64 / 16.0;
        let p = section.geometry.point_at(t);
        let d = p - origin;
        let radial = (d - axis * d.dot(&axis)).normalized().expect("off axis");
        let (u, _) = popkern_core::geometry::uv_projection::project_point_to_uv(
            &SurfaceGeometry::Nurbs(cap.clone()),
            &p,
        )
        .expect("projects onto the cap");
        let n = cap.normal_at(u, 0.0).normalized().expect("regular");
        worst_g1 = worst_g1.max(n.cross(&radial).length());
        // Normal curvature along the meridian (v) against the ruling's 0.
        let h = 1e-4;
        let (s0, s1, s2) = (
            cap.point_at(u, 0.0),
            cap.point_at(u, h),
            cap.point_at(u, 2.0 * h),
        );
        let d1: Vec3 = s1 - s0;
        let kappa = ((s2 - s1) - d1).dot(&n) / d1.dot(&d1);
        worst_g2 = worst_g2.max(kappa.abs() * r);
    }
    assert!(
        worst_g1 < 1e-9,
        "G1 across the section: |n x n| {worst_g1:.3e}"
    );
    assert!(
        worst_g2 < 1e-2,
        "G2 across the section: relative dk {worst_g2:.3e}"
    );
}
