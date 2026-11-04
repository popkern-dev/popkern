//! The per-face data the edge-face scan reads for every edge crossing
//! a face: the solver's parameter window, the trim box a hit must land
//! in, and how periodicity relaxes it. Each depends on the face alone,
//! and on a face with many holes each walks every vertex: a cap of the
//! 64-slot disc recomputed them for each tool edge crossing it
//! (PK-1032), so the scan computes them once per face.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::intersect::bounded_parameter_range;
use crate::pcurve_box::UvBox;

/// What the edge-face scan reads of one face.
#[derive(Debug, Clone, Copy)]
pub(crate) struct FaceWindow {
    /// The solver's u window.
    pub u: (f64, f64),
    /// The solver's v window.
    pub v: (f64, f64),
    /// The box an accepted hit's (u, v) must lie in; `None` gives no
    /// verdict.
    pub trim_box: Option<UvBox>,
    /// The period a hit's u is tested modulo, `None` when u is
    /// aperiodic or the face carries its own seam.
    pub u_period: Option<f64>,
    /// The face carries its own seam: its u is not filtered.
    pub full_turn: bool,
}

/// The window of `face`, whose surface is `surface`.
pub(crate) fn face_window(
    brep: &BRep,
    face: FaceHandle,
    surface: &SurfaceGeometry,
    tolerance: &Tolerance,
) -> FaceWindow {
    let ((mut u_min, mut u_max), (v_min, v_max)) = bounded_parameter_range(brep, face, surface);
    // CLOSED-WALL CAVEAT (same family as the trim-box caveat
    // below): a periodic wall bounded only by rims has 1–2
    // loop vertices, ALL on the seam, so the vertex-derived u
    // window collapses to ±padding around the seam azimuth
    // and prunes legitimate far-side roots (measured: the
    // cross-drilled boss — each seam-line × wall pair lost
    // the root opposite the seam, the pierce curves never
    // split at the missing nodes, and the graze filter
    // cascaded ALL record pieces off both walls). A closed
    // wall spans the full period: widen u to one full turn.
    let distinct_positions = {
        let ps = brep.face_outer_loop_positions(face);
        let mut distinct: Vec<popkern_core::geometry::Point3> = Vec::new();
        for p in ps {
            if distinct
                .iter()
                .all(|q| q.distance_to(&p) > tolerance.linear())
            {
                distinct.push(p);
            }
        }
        distinct.len()
    };
    {
        use popkern_core::geometry::SurfaceEval as _;
        if let (Some(period), true) = (surface.parameter_periodicity().0, distinct_positions < 3)
            && u_max - u_min < period
        {
            let mid = 0.5 * (u_min + u_max);
            u_min = mid - 0.5 * period;
            u_max = mid + 0.5 * period;
        }
    }

    // Acceptance is trim-honest: the padded range above gives the
    // SOLVER headroom, but a hit landing in the padding is not on
    // the face (measured: a seat tool's wall plane crossed a deck
    // outline edge 4 mm outside the wall's bounds; the phantom hit
    // split the cellular graph edge there and the polluted loop
    // broke the NEXT boolean on the body). Bound acceptance by the
    // face's exact vertex-derived parameter box plus tolerance.
    // CLOSED-RIM CAVEAT: a face bounded by closed edges (a
    // cylinder cap's full-circle rim) has 1–2 loop vertices,
    // so its vertex-derived box collapses and would reject
    // every legitimate hit (measured: star-socket corner
    // hits on a cylinder cap all filtered, 427 open edges).
    // Fewer than 3 vertices → the pcurve window instead
    // (PK-921), or no verdict without pcurves.
    // Distinct-position guard (same reason as the widened u
    // window above): a closed wall's loop visits its 1–2 seam
    // vertices twice, so the RAW position count lies about
    // the trim information available.
    let face_positions = crate::intersect_edge_face_coincident::trim_box_positions(brep, face);
    let trim_box = if distinct_positions >= 3 {
        Some(crate::intersect_edge_face_coincident::vertex_trim_box(
            brep,
            face,
            &face_positions,
            surface,
        ))
    } else {
        crate::intersect_edge_face_coincident::pcurve_trim_box(brep, face)
    };
    // A face that carries its own SEAM (an edge whose two
    // half-edges both bound it) spans the whole u period, but
    // its seam vertices project to u = 0 only, never 2π: the
    // vertex box stops at the last vertex before the seam and
    // rejects the strip beyond it. Measured (PK-746): a fin
    // straddling a housing's seam put a root corner at
    // u = 6.2427; a boss rod's seam edge pierced the housing
    // at u = 6.2785, the hit was rejected, its ring carried
    // no node at the rod's seam, the rod's chart walk fell
    // back to 3D-Newell, and the union served without the
    // boss wall. Such a face's u is not filtered.
    let u_period = {
        use popkern_core::geometry::SurfaceEval as _;
        surface.parameter_periodicity().0
    }
    .filter(|_| !crate::intersect_edge_face_coincident::face_has_own_seam(brep, face));
    let full_turn = u_period.is_none() && {
        use popkern_core::geometry::SurfaceEval as _;
        surface.parameter_periodicity().0.is_some()
    };
    FaceWindow {
        u: (u_min, u_max),
        v: (v_min, v_max),
        trim_box,
        u_period,
        full_turn,
    }
}
