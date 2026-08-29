//! Interior spoke seam edge-triple authoring from committed cap
//! surfaces (PK-240).
//!
//! The N-quad poke writes straight-Line placeholder curves on the
//! interior spoke edges; with straight-spoke fans the committed cap
//! sides coincide with those lines within validation tolerance, so
//! the placeholders survived to serving unnoticed. The balanced-fan
//! retry curves the spokes (measured side bow 1e-3..4e-3 at
//! 1:20-1:50 radius asymmetry), and the surviving Line placeholder
//! then reads as `CurveSurfaceDeviation` against both cap surfaces.
//!
//! This pass derives each interior seam's edge triple — the 3-D
//! curve plus both half-edge pcurves — from the COMMITTED cap
//! geometry: the shared patch side is a degree-9 Bezier row of the
//! binonic grid (ring G0 across the shared side is gated at 1e-10
//! before commit, so either patch's row is the seam), and the
//! pcurves are exact UV iso lines. Correct by construction, no
//! sampling or fitting.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::pcurve::{Pcurve, PcurveGeometry};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{NurbsCurve, Point3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// One patch side: 10 boundary CPs plus its UV iso endpoints,
/// oriented along increasing row/column index.
struct PatchSide {
    cps: [Point3; 10],
    uv_start: (f64, f64),
    uv_end: (f64, f64),
}

/// Author every interior seam shared by two of `sub_faces` from the
/// committed binonic cap surfaces. Returns the number of seams
/// authored; seams whose faces are not binonic caps or whose
/// endpoints match no patch side are left unchanged.
pub(super) fn author_interior_spoke_seam_curves_from_committed_caps(
    brep: &mut BRep,
    sub_faces: &[FaceHandle],
) -> usize {
    let set: std::collections::BTreeSet<FaceHandle> = sub_faces.iter().copied().collect();
    let tol = brep.tolerance.linear().max(1.0e-12) * 10.0;
    let edges: Vec<EdgeHandle> = brep.edges.keys().collect();
    let mut authored = 0;
    for eh in edges {
        let Some(edge) = brep.edges.get(eh) else {
            continue;
        };
        let he = edge.half_edge;
        let Some(he_data) = brep.half_edges.get(he) else {
            continue;
        };
        let twin = he_data.twin;
        let Some(twin_data) = brep.half_edges.get(twin) else {
            continue;
        };
        let (Some(loop_a), Some(loop_b)) = (
            brep.loops.get(he_data.loop_),
            brep.loops.get(twin_data.loop_),
        ) else {
            continue;
        };
        let (face_a, face_b) = (loop_a.face, loop_b.face);
        if face_a == face_b || !set.contains(&face_a) || !set.contains(&face_b) {
            continue;
        }
        let (Some(origin), Some(dest)) = (
            brep.vertices.get(he_data.origin).map(|v| v.position),
            brep.vertices.get(twin_data.origin).map(|v| v.position),
        ) else {
            continue;
        };
        let Some(side_a) = matching_binonic_side(brep, face_a, origin, dest, tol) else {
            continue;
        };
        let Some(side_b) = matching_binonic_side(brep, face_b, origin, dest, tol) else {
            continue;
        };
        // The seam curve, oriented from the edge's canonical
        // half-edge origin to its destination (both sides carry the
        // same geometry to the ring G0 gate; side A is the author).
        let curve = NurbsCurve {
            degree: 9,
            knots: clamped_deg9_knots(),
            control_points: side_a.cps.to_vec(),
            weights: vec![1.0; 10],
        };
        if let Some(edge_mut) = brep.edges.get_mut(eh) {
            edge_mut.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(curve),
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        // Exact iso pcurves in each face's chart, parameterized by
        // the edge curve's t (both half-edges share the edge-curve
        // parameterization per the pcurve contract).
        if let Some(he_mut) = brep.half_edges.get_mut(he) {
            he_mut.pcurve = Some(Pcurve {
                geometry: PcurveGeometry::Linear {
                    start: side_a.uv_start,
                    end: side_a.uv_end,
                },
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        if let Some(twin_mut) = brep.half_edges.get_mut(twin) {
            twin_mut.pcurve = Some(Pcurve {
                geometry: PcurveGeometry::Linear {
                    start: side_b.uv_start,
                    end: side_b.uv_end,
                },
                t_start: 0.0,
                t_end: 1.0,
            });
        }
        authored += 1;
    }
    authored
}

/// The side of `face`'s committed binonic patch whose end control
/// points match `(origin, dest)` — oriented origin → dest.
fn matching_binonic_side(
    brep: &BRep,
    face: FaceHandle,
    origin: Point3,
    dest: Point3,
    tol: f64,
) -> Option<PatchSide> {
    let surface = brep.faces.get(face)?.surface.as_ref()?;
    let SurfaceGeometry::Gregory(patch) = &surface.geometry else {
        return None;
    };
    let GregoryPatch::RectangularBiNonic(binonic) = patch else {
        return None;
    };
    let grid = &binonic.grid;
    let sides = [
        PatchSide {
            cps: std::array::from_fn(|j| grid[0][j]),
            uv_start: (0.0, 0.0),
            uv_end: (0.0, 1.0),
        },
        PatchSide {
            cps: std::array::from_fn(|j| grid[9][j]),
            uv_start: (1.0, 0.0),
            uv_end: (1.0, 1.0),
        },
        PatchSide {
            cps: std::array::from_fn(|i| grid[i][0]),
            uv_start: (0.0, 0.0),
            uv_end: (1.0, 0.0),
        },
        PatchSide {
            cps: std::array::from_fn(|i| grid[i][9]),
            uv_start: (0.0, 1.0),
            uv_end: (1.0, 1.0),
        },
    ];
    for side in sides {
        let (s, e) = (side.cps[0], side.cps[9]);
        if s.distance_to(&origin) <= tol && e.distance_to(&dest) <= tol {
            return Some(side);
        }
        if s.distance_to(&dest) <= tol && e.distance_to(&origin) <= tol {
            let mut cps = side.cps;
            cps.reverse();
            return Some(PatchSide {
                cps,
                uv_start: side.uv_end,
                uv_end: side.uv_start,
            });
        }
    }
    None
}

fn clamped_deg9_knots() -> Vec<f64> {
    let mut knots = vec![0.0; 10];
    knots.extend(std::iter::repeat_n(1.0, 10));
    knots
}
