//! Triple-point refinement of boundary-edge × surface junctions
//! (PK-480).
//!
//! The exact-junction inventories (`pierce_crossing_params`,
//! `far_junction`) solve each junction as a 1D crossing of a
//! boundary EDGE CURVE with the other pair surface. That solve takes
//! the FITTED edge curve as position authority, so the junction
//! slides ALONG the edge by `edge fit slop / sin(crossing angle)`.
//! At healthy angles the slide stays inside the dead-end weld band
//! (the serving m6-nut class); at near-tangent crossings it explodes:
//! a tap entry cone's coaxial cone∩cylinder circles cross the helical
//! thread-crest edges at the ~3.5° lead angle, a ~1e-4 crest fit sag
//! became a 6.2e-3 junction slide, the cone×flank pair solved the
//! SAME triple point 6.2e-3 away, the crest edge was split at both
//! copies, the cone chart dead-ended at each, and the interior-graze
//! cascade stranded the whole chamfer window (honest MergeFailed on
//! the deep-plunged angled/perpendicular tap).
//!
//! The true junction is the common point of THREE surfaces — the
//! boundary edge's two adjacent face surfaces and the other pair
//! surface — which are mutually transversal at a genuine triple
//! point. A Newton on the three signed surface distances converges to
//! it from the 1D seed without any 1/sin amplification: the result is
//! accurate to the surface fits themselves (~1e-4 class, inside the
//! weld band). Every consumer downstream (window bounds, piece-end
//! retrims and pins, boundary-edge splits) inherits the refined
//! point, so the face pairs meeting at the junction agree by
//! construction.
//!
//! Honesty bounds: the refinement is kept only when the Newton
//! converges (residual to all three surfaces under the acceptance
//! gate) and the move from the seed stays inside the amplified-slide
//! cap; a seam edge (both sides one face), a degenerate normal
//! triple (near-parallel surfaces — the Jacobian goes singular), or a
//! diverging iteration keeps the 1D seed unchanged.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle};

/// Cap on how far the refinement may move a junction off its 1D
/// seed, as a multiple of the linear tolerance: covers the measured
/// near-tangent slide class (6.2e-3 at 1e-6 linear) with margin while
/// staying well under distinct-junction spacing (millimetres).
const MOVE_CAP_LINEAR: f64 = 20_000.0;

/// Acceptance residual for convergence, as a multiple of the linear
/// tolerance: the three surfaces genuinely meet at the point, so a
/// converged Newton sits orders of magnitude tighter; this gate only
/// rejects walks that stalled somewhere that is not a triple point.
const RESIDUAL_GATE_LINEAR: f64 = 20.0;

/// Refine an edge×surface junction onto the true triple point of the
/// edge's two adjacent face surfaces and `other_surf`. Returns the
/// seed unchanged whenever the configuration is out of scope or the
/// iteration fails a bound (see the module docs).
pub(crate) fn refine_to_triple_point(
    brep_x: &BRep,
    eh: EdgeHandle,
    other_surf: &SurfaceGeometry,
    seed: Point3,
    tolerance: &Tolerance,
) -> Point3 {
    let Some(edge) = brep_x.edges.get(eh) else {
        return seed;
    };
    let Some(he) = brep_x.half_edges.get(edge.half_edge) else {
        return seed;
    };
    let mut faces = [None, None];
    for (slot, h) in [(0, edge.half_edge), (1, he.twin)] {
        faces[slot] = brep_x
            .half_edges
            .get(h)
            .and_then(|hd| brep_x.loops.get(hd.loop_))
            .map(|l| l.face);
    }
    let (Some(fa), Some(fb)) = (faces[0], faces[1]) else {
        return seed;
    };
    // A seam edge bounds one face twice: only two distinct surfaces
    // meet there and the triple is degenerate by construction.
    if fa == fb {
        return seed;
    }
    let surf_of = |f| {
        brep_x
            .faces
            .get(f)
            .and_then(|face: &popkern_core::topology::Face| face.surface.as_ref())
            .map(|s| &s.geometry)
    };
    let (Some(s1), Some(s2)) = (surf_of(fa), surf_of(fb)) else {
        return seed;
    };
    let surfs: [&SurfaceGeometry; 3] = [s1, s2, other_surf];
    let lin = tolerance.linear();
    let residual_gate = lin * RESIDUAL_GATE_LINEAR;
    // The correction this Newton exists for is the along-edge slide
    // of the 1D solve: magnitude ≈ (seed residual to the missed
    // surface) / sin(crossing angle). Cap the move ADAPTIVELY at 32×
    // the seed's worst surface residual (crossings ≥ ~1.8°), bounded
    // by the tolerance-scaled ceiling — a pure tolerance multiple
    // explodes on coarse-tolerance corpora (linear = 1e-3 gives a
    // 20-unit cap on the multicurve cubes).
    let seed_residual = {
        let mut worst = 0.0_f64;
        for s in &surfs {
            let Some((u, v)) = s.closest_point(&seed) else {
                return seed;
            };
            worst = worst.max(s.point_at(u, v).distance_to(&seed));
        }
        worst
    };
    let move_cap = (seed_residual * 32.0).min(lin * MOVE_CAP_LINEAR);
    if move_cap <= 0.0 {
        return seed;
    }
    let jt = popkern_core::trace::trace_enabled("junction");
    let mut p = seed;
    let mut converged = false;
    for _ in 0..16 {
        let mut normals = [Vec3::new(0.0, 0.0, 0.0); 3];
        let mut dist = [0.0_f64; 3];
        for (i, s) in surfs.iter().enumerate() {
            let Some((u, v)) = s.closest_point(&p) else {
                if jt {
                    eprintln!(
                        "    [triple-junction] FAIL closest_point surf {i} at ({:.4},{:.4},{:.4})",
                        p.x, p.y, p.z
                    );
                }
                return seed;
            };
            let foot = s.point_at(u, v);
            let Some(n) = s.normal_at(u, v).normalized() else {
                return seed;
            };
            normals[i] = n;
            dist[i] = (p - foot).dot(&n);
        }
        if dist.iter().all(|d| d.abs() < residual_gate * 1e-3) {
            converged = true;
            break;
        }
        // Solve normals · delta = -dist (3×3 by Cramer).
        let det = normals[0].dot(&normals[1].cross(&normals[2]));
        if det.abs() < 1e-6 {
            if jt {
                eprintln!(
                    "    [triple-junction] FAIL degenerate det {det:.2e} at ({:.4},{:.4},{:.4})",
                    p.x, p.y, p.z
                );
            }
            return seed;
        }
        let rhs = Vec3::new(-dist[0], -dist[1], -dist[2]);
        // delta = N⁻¹ · rhs for the row matrix N of normals: the
        // inverse's columns are the scaled cross products of the
        // opposite normal pairs (adjugate over the row triple).
        let r0 = normals[1].cross(&normals[2]) * (1.0 / det);
        let r1 = normals[2].cross(&normals[0]) * (1.0 / det);
        let r2 = normals[0].cross(&normals[1]) * (1.0 / det);
        let delta = r0 * rhs.x + r1 * rhs.y + r2 * rhs.z;
        p = p + delta;
        if p.distance_to(&seed) > move_cap {
            if jt {
                eprintln!(
                    "    [triple-junction] FAIL move cap from ({:.4},{:.4},{:.4})",
                    seed.x, seed.y, seed.z
                );
            }
            return seed;
        }
        if delta.length() < 1e-13 {
            converged = true;
            break;
        }
    }
    if !converged {
        // A final residual check: accept a slow but arrived walk.
        converged = surfs.iter().all(|s| {
            s.closest_point(&p)
                .map(|(u, v)| s.point_at(u, v).distance_to(&p) < residual_gate)
                .unwrap_or(false)
        });
        if !converged && jt {
            eprintln!(
                "    [triple-junction] FAIL residual after iterations from ({:.4},{:.4},{:.4})",
                seed.x, seed.y, seed.z
            );
        }
    }
    if converged && p.distance_to(&seed) <= move_cap {
        if popkern_core::trace::trace_enabled("junction") {
            eprintln!(
                "    [triple-junction] ({:.6},{:.6},{:.6}) -> ({:.6},{:.6},{:.6}) moved {:.2e}",
                seed.x,
                seed.y,
                seed.z,
                p.x,
                p.y,
                p.z,
                p.distance_to(&seed)
            );
        }
        p
    } else {
        seed
    }
}

#[cfg(test)]
#[path = "clip_triple_junction_tests.rs"]
mod clip_triple_junction_tests;
