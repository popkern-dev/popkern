//! Folded-seam unfolding and chain-continuity period lifting for the
//! boundary-first ring builder (PK-457), plus the degenerate-sliver
//! drop its pole crossing requires. Split from `boundary_first.rs` /
//! `boundary_first_edge_samples.rs` (500-line cap).

use popkern_core::topology::BRep;

use crate::mesh::MeshVertex;

/// Chain-lift periods for a boundary ring (PK-457): the geometric or
/// analytic closure periods of the loop's face chart, but only for
/// MIXED loops — the same scope as the PK-454 closed-edge entry rule
/// (all-closed ring-tube families carry pcurves authored as one
/// consistent family and must stay untouched).
pub(crate) fn ring_chain_periods(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    loop_has_open: bool,
) -> (Option<f64>, Option<f64>) {
    if !loop_has_open {
        return (None, None);
    }
    brep.loops
        .get(loop_handle)
        .and_then(|l| brep.faces.get(l.face))
        .and_then(|f| f.surface.as_ref())
        .map(|sof| {
            use popkern_core::geometry::SurfaceEval;
            let analytic = sof.geometry.parameter_periodicity();
            let chart =
                popkern_core::geometry::chart_periods::geometric_chart_periods(&sof.geometry);
            (analytic.0.or(chart.0), analytic.1.or(chart.1))
        })
        .unwrap_or((None, None))
}

/// Lift one refined boundary sample to the period branch nearest its
/// predecessor (PK-457). Pcurve authoring keeps seam pairs FOLDED
/// onto one border (PK-441 — recomputation must not take chart
/// periods), so a mixed loop's refined uv chain can jump a whole
/// period between edges; smooth chains round to zero shifts and are
/// untouched.
pub(crate) fn lift_to_chain(
    (mut u, mut v): (f64, f64),
    prev: Option<(f64, f64)>,
    periods: (Option<f64>, Option<f64>),
) -> (f64, f64) {
    if let Some((pu, pv)) = prev {
        if let Some(p) = periods.0 {
            u += ((pu - u) / p).round() * p;
        }
        if let Some(p) = periods.1 {
            v += ((pv - v) / p).round() * p;
        }
    }
    (u, v)
}

/// Fold shift for a seam leg at emission time (PK-457): pcurve
/// authoring keeps seam pairs FOLDED onto one border (PK-441), so
/// when a ring re-enters an edge it has already walked (the second
/// leg of a seam pair) at the SAME border, the leg belongs one period
/// away — but only when the junction into it is chart-DEGENERATE (a
/// pole row: the period jump is geometrically free there). Left
/// folded, the ring doubles back as a slit, the dissolvers amputate
/// it, and the Delaunay coverage check requests splits forever (the
/// encroachment runaway measured on every revolved-NURBS dimple cap:
/// split counts doubling to the 512/edge cap). Deciding DURING the
/// walk lets every later direction/continuity choice (the PK-454
/// closed-edge entry rule, the chain lift) see the unfolded chart.
pub(crate) fn fold_shift_for_leg(
    surface: &popkern_core::geometry::SurfaceGeometry,
    first_leg_start: Option<(f64, f64)>,
    leg_start: (f64, f64),
    prev_point: Option<(f64, f64)>,
    periods: (Option<f64>, Option<f64>),
) -> (f64, f64) {
    use popkern_core::geometry::SurfaceEval;
    let (Some(first), Some(_prev)) = (first_leg_start, prev_point) else {
        return (0.0, 0.0);
    };
    for (axis, period) in [(0usize, periods.0), (1usize, periods.1)] {
        let Some(p) = period else { continue };
        let pick = |q: (f64, f64)| if axis == 0 { q.0 } else { q.1 };
        if (pick(leg_start) - pick(first)).abs() > 1.0e-3 * p {
            continue; // already on the other border
        }
        // Degenerate junction: the leg's own first sample is the
        // shared vertex the period jump crosses (a revolution cap's
        // pole); stepping the periodic axis there must not move in 3D
        // (the pole-row test the chart walk uses). A seam pair
        // meeting at a REGULAR vertex is a genuine slit and stays
        // folded for the dissolvers.
        let p0 = surface.point_at(leg_start.0, leg_start.1);
        let (du, dv) = if axis == 0 {
            (1.0e-4 * p, 0.0)
        } else {
            (0.0, 1.0e-4 * p)
        };
        let p1 = surface.point_at(leg_start.0 + du, leg_start.1 + dv);
        let scale = 1.0 + p0.distance_to(&popkern_core::geometry::Point3::new(0.0, 0.0, 0.0));
        if p0.distance_to(&p1) < 1.0e-11 * scale {
            return if axis == 0 { (p, 0.0) } else { (0.0, p) };
        }
    }
    (0.0, 0.0)
}

/// Drop triangles with a zero-length edge (PK-455/PK-457): a ring
/// carrying BOTH chart images of a degenerate (pole) vertex satisfies
/// Delaunay coverage with a triangle whose pole edge has zero 3D
/// length — a (P, P) half-edge no neighbour can cancel, the same
/// sliver class the uniform-grid path drops. Coverage is checked
/// BEFORE this runs; the slivers are dropped from the output only.
pub(crate) fn drop_degenerate_slivers(triangles: &mut Vec<[usize; 3]>, vertices: &[MeshVertex]) {
    triangles.retain(|t| {
        let (pa, pb, pc) = (
            vertices[t[0]].position,
            vertices[t[1]].position,
            vertices[t[2]].position,
        );
        let (ab, bc, ca) = (
            pa.distance_to(&pb),
            pb.distance_to(&pc),
            pc.distance_to(&pa),
        );
        ab.min(bc).min(ca) >= 1.0e-9 * ab.max(bc).max(ca)
    });
}

/// Canonicalize near-duplicate CONSECUTIVE ring positions (arc
/// joints: the two arcs meeting at an anchor evaluate the shared
/// point ~1e-17 apart, and the position-weld keys are BIT patterns).
/// Both entries snap to the pair's lexicographically smaller bit
/// pattern — a direction-independent choice, so the two faces sharing
/// the edges pick the SAME canonical point even though they traverse
/// in opposite orders (a keep-first dedupe diverged and left an
/// unpaired sliver at a strip-anchored rim, PK-344).
pub(crate) fn weld_consecutive_positions(pos: &mut [popkern_core::geometry::Point3]) {
    if pos.is_empty() {
        return;
    }
    let scale = pos
        .iter()
        .zip(pos.iter().cycle().skip(1))
        .map(|(a, b)| a.distance_to(b))
        .fold(0.0_f64, f64::max)
        .max(1e-12);
    let weld = scale * 1e-9;
    let n = pos.len();
    for i in 0..n {
        let j = (i + 1) % n;
        let (a, b) = (pos[i], pos[j]);
        let bits =
            |p: &popkern_core::geometry::Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
        if bits(&a) != bits(&b) && a.distance_to(&b) <= weld {
            let canonical = if bits(&a) <= bits(&b) { a } else { b };
            pos[i] = canonical;
            pos[j] = canonical;
        }
    }
}
