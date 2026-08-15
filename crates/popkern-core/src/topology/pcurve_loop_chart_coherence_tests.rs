//! Tests for per-loop chart coherence (PK-286): the invariant is
//! structural — enforced by `recompute_face_pcurves` itself — and
//! holds along every periodic chart direction.

use crate::geometry::pcurve::{Pcurve, PcurveGeometry};
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle};
use std::f64::consts::TAU;

fn lateral_face(brep: &BRep) -> FaceHandle {
    brep.faces
        .iter()
        .find(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Cylinder { .. })
            )
        })
        .map(|(fh, _)| fh)
        .expect("cylinder has a lateral face")
}

fn outer_loop_half_edges(brep: &BRep, fh: FaceHandle) -> Vec<HalfEdgeHandle> {
    let start = brep.loops[brep.faces[fh].outer_loop].half_edge;
    let mut hes = vec![start];
    let mut cursor = brep.half_edges[start].next;
    while cursor != start {
        hes.push(cursor);
        cursor = brep.half_edges[cursor].next;
    }
    hes
}

/// Endpoint coordinates of a half-edge's pcurve in traversal order.
fn pcurve_span(brep: &BRep, heh: HalfEdgeHandle) -> (Point3, Point3) {
    let pc = brep.half_edges[heh].pcurve.as_ref().expect("pcurve");
    let (t0, t1) = pc.parameter_range();
    (pc.point_at(t0), pc.point_at(t1))
}

/// Max |prev end − next start| gap around a loop, per UV coordinate.
fn loop_chain_gaps(brep: &BRep, fh: FaceHandle) -> (f64, f64) {
    let hes = outer_loop_half_edges(brep, fh);
    let mut max_u = 0.0_f64;
    let mut max_v = 0.0_f64;
    for i in 0..hes.len() {
        let (_, prev_end) = pcurve_span(brep, hes[i]);
        let (next_start, _) = pcurve_span(brep, hes[(i + 1) % hes.len()]);
        max_u = max_u.max((prev_end.x - next_start.x).abs());
        max_v = max_v.max((prev_end.y - next_start.y).abs());
    }
    (max_u, max_v)
}

/// A DIRECT per-face recompute must leave the loop branch-coherent:
/// the closure is part of `recompute_face_pcurves`' postcondition, not
/// a batch-wrapper courtesy. (Measured before the fix: the PK-159
/// parallel-seam completion recomputed two strip faces directly and
/// meshed 12 cracked weld edges when one rim arc landed a period
/// below its loop-mates.)
#[test]
fn direct_face_recompute_restores_loop_branch_coherence() {
    let mut brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let fh = lateral_face(&brep);
    let hes = outer_loop_half_edges(&brep, fh);
    // Push one loop member a whole period out of its loop's branch.
    // Its trace is identical (the chart is 2π-periodic), so the
    // survival rule alone would keep it forever.
    let victim = hes[hes.len() / 2];
    let shifted = brep.half_edges[victim]
        .pcurve
        .as_ref()
        .expect("cylinder lateral pcurves exist")
        .translated_uv(TAU, 0.0);
    brep.half_edges[victim].pcurve = Some(shifted);

    brep.recompute_face_pcurves(fh).unwrap();

    let (gap_u, _) = loop_chain_gaps(&brep, fh);
    assert!(
        gap_u < 1e-6,
        "loop must be branch-chained after a direct per-face recompute; max u gap {gap_u:.3e}"
    );
    // Principal cell: the loop as a whole sits within [0, 2π].
    for &heh in &hes {
        let (s, e) = pcurve_span(&brep, heh);
        for u in [s.x, e.x] {
            assert!(
                (-1e-6..=TAU + 1e-6).contains(&u),
                "loop member u={u} left the principal cell"
            );
        }
    }
}

/// The closure holds along V too (torus charts): the branch-honest
/// quarter-point chain aligns v on v-periodic surfaces, so a member
/// landing a period away in v must be re-chained exactly like u.
/// Unit-level: the surface is swapped for a torus and the loop's
/// pcurves are hand-authored — the closure reads only pcurves and
/// periodicity, never 3-D consistency.
#[test]
fn v_axis_closure_rechains_torus_loop() {
    let mut brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let fh = lateral_face(&brep);
    if let Some(sof) = brep.faces[fh].surface.as_mut() {
        sof.geometry = SurfaceGeometry::Torus {
            center: Point3::new(0.0, 0.0, 0.0),
            axis: Vec3::new(0.0, 0.0, 1.0),
            major_radius: 2.0,
            minor_radius: 0.5,
        };
    }
    let hes = outer_loop_half_edges(&brep, fh);
    assert_eq!(hes.len(), 4, "cylinder lateral loop is rim/seam/rim/seam");
    let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    for (i, &heh) in hes.iter().enumerate() {
        let a = corners[i];
        let b = corners[(i + 1) % corners.len()];
        brep.half_edges[heh].pcurve = Some(Pcurve {
            geometry: PcurveGeometry::Linear { start: a, end: b },
            t_start: 0.0,
            t_end: 1.0,
        });
    }
    let victim = hes[2];
    let shifted = brep.half_edges[victim]
        .pcurve
        .as_ref()
        .unwrap()
        .translated_uv(0.0, TAU);
    brep.half_edges[victim].pcurve = Some(shifted);

    brep.close_periodic_loop_winding(fh);

    let (gap_u, gap_v) = loop_chain_gaps(&brep, fh);
    assert!(
        gap_u < 1e-9 && gap_v < 1e-9,
        "torus loop must re-chain in v as in u; gaps u {gap_u:.3e} v {gap_v:.3e}"
    );
}

/// The mekr-hover distillation: an axial seam edge whose endpoint
/// projections straddle the chart wrap by ULP noise (atan2(±1e-12, ·)
/// picks the side) must classify as the iso line it is — the
/// branch-blind comparison read the straddle as a full-period sweep
/// and took the dense path, so a 7e-15 upstream change flipped the
/// authored pcurve degree-1 ↔ degree-3 (double_boss_plate).
#[test]
fn seam_hover_axial_edge_classifies_iso_line() {
    let mut brep = BRep::cylinder(1.0, 2.0, Tolerance::default()).unwrap();
    let fh = lateral_face(&brep);
    // The seam edge bounds the lateral face on both sides.
    let seam = brep
        .edges
        .iter()
        .find(|(_, e)| {
            let he = &brep.half_edges[e.half_edge];
            let twin = &brep.half_edges[he.twin];
            brep.loops[he.loop_].face == brep.loops[twin.loop_].face
        })
        .map(|(eh, _)| eh)
        .expect("cylinder has a seam edge");
    // Nudge one seam vertex across the wrap: its projection lands at
    // 2π − 1e-12 while the other endpoint and every interior sample
    // stay at +0.
    let (he, twin) = {
        let e = &brep.edges[seam];
        (e.half_edge, brep.half_edges[e.half_edge].twin)
    };
    let far_vertex = brep.half_edges[brep.half_edges[he].next].origin;
    let pos = brep.vertices[far_vertex].position;
    brep.vertices[far_vertex].position = Point3::new(pos.x, -1e-12, pos.z);
    brep.half_edges[he].pcurve = None;
    brep.half_edges[twin].pcurve = None;

    brep.recompute_face_pcurves(fh).unwrap();

    for heh in [he, twin] {
        let (s, e) = pcurve_span(&brep, heh);
        let du = (e.x - s.x).abs();
        assert!(
            du < 1e-9,
            "seam-hover axial edge must be an iso line in u, spans {du:.3e}"
        );
    }
}
