//! Fragment interior sampling for cellular classification: UV
//! scanline sample point selection (period-fold aware) and the
//! planar material-outward normal used for coorientation.
//!
//! Split from `cellular/boolean.rs` to keep it under the structural
//! line cap.

use popkern_core::topology::BRep;

use crate::merge::graph::IntersectionGraph;

/// Sample a point from a fragment's outer loop for classification.
/// Interior sample point for a fragment, aware of inner loops and
/// closed-rim outers. Works in the UV chart of the fragment's face:
/// builds pcurve-sampled polygons for the outer and inner loops,
/// scans across the outer polygon's centroid row, and picks the
/// midpoint of the widest interval covered by the outer but not by
/// any inner. Falls back to `None` when pcurves are missing.
pub(super) fn fragment_interior_sample(
    graph: &IntersectionGraph,
    fragment: &crate::merge::walk::Fragment,
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    orientation: bool,
    force_interior: bool,
) -> Option<popkern_core::geometry::Point3> {
    use popkern_core::geometry::surface::SurfaceEval;

    // `force_interior`: fragments of coincident-overlap faces must
    // sample strictly inside their region — the edge-midpoint
    // shortcut lands ON the other body's boundary corner, where the
    // offset-pair probes straddle an edge instead of a face and the
    // coorientation verdict inverts (measured: a flush relief's
    // flank wall probed at the tab's front corner read anti-oriented
    // and the mouth membrane was kept).
    let needs_care = force_interior
        || !fragment.inner_loops.is_empty()
        || fragment.outer_loop.iter().any(|&h| {
            let he = &graph.half_edges[h];
            he.origin == graph.half_edges[he.twin].origin
        });
    if !needs_care {
        return None;
    }

    // The per-half-edge UV run sampling (pcurve windows, projection
    // fallback, period-branch alignment) lives in
    // `fragment_sampling::loop_uv_runs`, shared with the cut-offset
    // sampler.
    let sample_loop = |hes: &[usize]| -> Option<Vec<(f64, f64)>> {
        let runs = super::fragment_sampling::loop_uv_runs(graph, &fragment.face, hes, surface)?;
        let poly: Vec<(f64, f64)> = runs.into_iter().flat_map(|(_, r)| r).collect();
        (poly.len() >= 3).then_some(poly)
    };
    // The chart walk's own loops come first: their runs sit on the
    // branches the walk embedded, pole rows included. Re-sampled, a
    // seam split by a hole straddling it projects both copies onto
    // one meridian (no pcurve survives an analytic split), the outer
    // polygon collapses onto the hole outline, and the kept sphere
    // face sampled inside its own hole (PK-937: a ball drilled across
    // its pole axis left the cellular path on every orientation).
    let (outer, inners) = match fragment
        .chart
        .as_ref()
        .filter(|c| c.loops.len() == 1 + fragment.inner_loops.len())
    {
        Some(chart) => {
            let mut loops = chart.loops.iter().map(|l| l.points().to_vec());
            let outer = loops.next()?;
            (outer, loops.collect::<Vec<_>>())
        }
        None => {
            let outer = sample_loop(&fragment.outer_loop)?;
            let inners: Vec<Vec<(f64, f64)>> = fragment
                .inner_loops
                .iter()
                .filter_map(|l| sample_loop(l))
                .collect();
            (outer, inners)
        }
    };
    if inners.len() != fragment.inner_loops.len() || outer.len() < 3 {
        return None;
    }

    // Scan row: outer centroid v.
    let cv = outer.iter().map(|p| p.1).sum::<f64>() / outer.len() as f64;
    let crossings_at = |poly: &[(f64, f64)], v: f64| -> Vec<f64> {
        let mut xs = Vec::new();
        let n = poly.len();
        for i in 0..n {
            let (u0, v0) = poly[i];
            let (u1, v1) = poly[(i + 1) % n];
            if (v0 > v) != (v1 > v) {
                xs.push(u0 + (v - v0) / (v1 - v0) * (u1 - u0));
            }
        }
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        xs
    };
    // DOUBLY periodic charts (tori) decide material exactly, by parity
    // from a reference known to lie in the region — a slit of the
    // loops, else the chart walk's handedness (`fragment_torus_sample`):
    // ray parity
    // has no escape there, and the nearest-boundary vote below is a
    // local guess. Fragments with neither reference, and loop sets that
    // bound no region, keep the vote, recorded under `face-selection`.
    if let (Some(pu), Some(pv)) = surface.parameter_periodicity() {
        match super::fragment_torus_sample::exact_torus_sample(
            graph,
            fragment,
            surface,
            cv,
            (pu, pv),
        ) {
            Ok((u, v)) => return Some(surface.point_at(u, v)),
            Err(reason) => {
                if popkern_core::trace::trace_enabled("face-selection") {
                    eprintln!("      [selpoly-torus] exact decision declined: {reason}");
                }
            }
        }
    }
    let outer_x = crossings_at(&outer, cv);
    if outer_x.len() < 2 {
        return None;
    }
    // PERIODIC WINDING GUARD: a fragment that winds around a periodic
    // chart yields a boundary polygon whose u-coordinates span a full
    // period or more, and its runs may sit on DIFFERENT period
    // branches (authored pcurves keep their stored branch, projected
    // runs chain from the walk position). Raw even-odd pairing of the
    // scan crossings then pairs across the branch gap and the sample
    // lands inside a cut window (measured: prism-through-cylinder
    // union — the wall fragment's window sides crossed the scan row
    // at u=0.253 and u=12.314, the paired midpoint 2π sat in the
    // boss footprint, and the whole wall classified Inside and was
    // dropped). Fold the crossings into one period cell — arcs
    // between adjacent fold-cell crossings partition the period
    // circle regardless of branch — and choose the arc whose midpoint
    // lies on the fragment's MATERIAL side, decided by a local 3D
    // test against the nearest boundary sample (interior = left of
    // the walk direction under the face's oriented normal). Width
    // alone cannot decide: on a corner-breached bore wall the OUTSIDE
    // region's widest fold arc ran through the kept band and the
    // whole shank classified Inside (measured: tilted_hole_plate).
    let u_period = {
        use popkern_core::geometry::surface::SurfaceEval as _;
        surface.parameter_periodicity().0
    };
    let (u_lo, u_hi) = outer.iter().fold((f64::MAX, f64::MIN), |(lo, hi), q| {
        (lo.min(q.0), hi.max(q.0))
    });
    if let Some(p) = u_period.filter(|p| u_hi - u_lo >= p * (1.0 - 1e-9)) {
        let anchor = outer_x[0];
        let fold = |u: f64| -> f64 { u - p * ((u - anchor) / p).floor() };
        // Open-v charts (cylinders, cones, closed-profile
        // extrusions): decide every fold arc by exact even-odd
        // containment over all period lifts of the loops, and cancel
        // coincident crossings in pairs (a seam edge's u = 0 and
        // u = period copies). See `crate::periodic_chart_parity`.
        if surface.parameter_periodicity().1.is_none() {
            use crate::periodic_chart_parity::{
                cancel_coincident_pairs, contains_on_u_periodic_chart_with_poles,
            };
            let pole_row = |v: f64| crate::merge::walk::chart_poles::row_collapses(surface, 0.0, v);
            let raw: Vec<f64> = outer_x
                .iter()
                .copied()
                .chain(inners.iter().flat_map(|i| crossings_at(i, cv)))
                .map(fold)
                .collect();
            let xs = cancel_coincident_pairs(raw, p * 1e-9);
            let loops: Vec<&[(f64, f64)]> = std::iter::once(outer.as_slice())
                .chain(inners.iter().map(|i| i.as_slice()))
                .collect();
            let arcs: Vec<(f64, f64)> = if xs.is_empty() {
                vec![(anchor, anchor + p)]
            } else {
                (0..xs.len())
                    .map(|i| {
                        let b = xs.get(i + 1).copied().unwrap_or(xs[0] + p);
                        (xs[i], b)
                    })
                    .collect()
            };
            let mut best: Option<(f64, f64)> = None;
            for (a, b) in arcs {
                let mid = 0.5 * (a + b);
                let material =
                    contains_on_u_periodic_chart_with_poles(&loops, mid, cv, p, &pole_row);
                if popkern_core::trace::trace_enabled("face-selection") {
                    eprintln!(
                        "      [selpoly-arc] [{a:.3},{b:.3}] mid={mid:.3} contains={material}"
                    );
                }
                if material && best.map(|(bw, _)| b - a > bw).unwrap_or(true) {
                    best = Some((b - a, mid));
                }
            }
            let (_, mid_u) = best?;
            if popkern_core::trace::trace_enabled("face-selection") {
                eprintln!(
                    "      [selpoly] periodic-parity scan cv={cv:.3} xs={xs:?} -> mid_u={mid_u:.3}"
                );
            }
            return Some(surface.point_at(mid_u, cv));
        }
        let inner_folded: Vec<Vec<f64>> = inners
            .iter()
            .map(|i| crossings_at(i, cv).into_iter().map(fold).collect())
            .collect();
        // Partition the period circle by OUTER and INNER crossings
        // alike: an arc between adjacent crossings is then uniformly
        // material and uniformly in/out of every inner loop, so the
        // midpoint tests below are sound. With outer crossings only,
        // an arc spanning most of the wall was rejected whole because
        // its midpoint happened to land inside an inner-loop hole,
        // and the sole surviving arc was the seam-straddling cut
        // window (measured: cross-drilled cylinder — the kept wall
        // sampled at the entry-hole centre, classified Inside, and
        // the cellular body shipped five half-claimed edges).
        let mut xs: Vec<f64> = outer_x
            .iter()
            .copied()
            .chain(inner_folded.iter().flatten().copied())
            .map(fold)
            .collect();
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        xs.dedup_by(|b, a| (*b - *a).abs() < p * 1e-9);
        // 3D boundary samples with walk-direction tangents for the
        // material-side test (the polygon is in walk order).
        let boundary_3d: Vec<(popkern_core::geometry::Point3, popkern_core::geometry::Vec3)> = {
            let n = outer.len();
            (0..n)
                .map(|i| {
                    let (u0, v0) = outer[i];
                    let (u1, v1) = outer[(i + 1) % n];
                    let q = surface.point_at(u0, v0);
                    let q1 = surface.point_at(u1, v1);
                    (q, q1 - q)
                })
                .collect()
        };
        let material_side = |mid_u: f64| -> Option<bool> {
            let pt = surface.point_at(mid_u, cv);
            let (mut best_d, mut best_i) = (f64::MAX, usize::MAX);
            for (i, (q, t)) in boundary_3d.iter().enumerate() {
                // Skip degenerate segments (branch-glue joints map to
                // one 3D point) — their tangent carries no direction.
                if t.length_squared() < 1e-24 {
                    continue;
                }
                let d = pt.distance_to(q);
                if d < best_d {
                    best_d = d;
                    best_i = i;
                }
            }
            let (q, t) = boundary_3d.get(best_i)?;
            let (qu, qv) = outer[best_i];
            let mut n = surface.normal_at(qu, qv);
            if !orientation {
                n = -n;
            }
            let interior_dir = n.cross(t);
            if interior_dir.length_squared() < 1e-30 {
                return None;
            }
            // FALLBACK ONLY (PK-685): reached on doubly periodic charts
            // where the exact decision (`fragment_torus_sample`) had no
            // reference. Runs are geometrically oriented, so region-LEFT
            // of travel under the oriented normal is the rule for both
            // operands — measured on oriented runs against exact parity
            // (u-periodic open-v fold arcs, bool + examples + fastener
            // suites): left agrees on 95% of A arcs and 84% of B arcs.
            // The retired operand rule ("B right") compensated for runs
            // left in bookkeeping order and agreed on 16% of B arcs
            // once they were oriented.
            Some((pt - *q).dot(&interior_dir) > 0.0)
        };
        if popkern_core::trace::trace_enabled("face-selection") {
            eprintln!("      [selpoly-fold] orientation={orientation}");
        }
        let mut best: Option<(f64, f64)> = None;
        let n = xs.len();
        // (width, mid, inside_inner, side) per fold arc.
        let arcs: Vec<(f64, f64, bool, Option<bool>)> = (0..n)
            .map(|i| {
                let a = xs[i];
                let b = if i + 1 < n { xs[i + 1] } else { xs[0] + p };
                let mid = 0.5 * (a + b);
                let inside_inner = inner_folded
                    .iter()
                    .any(|xs| xs.iter().filter(|&&x| x < mid).count() % 2 == 1);
                let side = material_side(mid);
                if popkern_core::trace::trace_enabled("face-selection") {
                    eprintln!(
                        "      [selpoly-arc] [{a:.3},{b:.3}] mid={mid:.3} inner={inside_inner} side={side:?}"
                    );
                }
                (b - a, mid, inside_inner, side)
            })
            .collect();
        // Every crossing — outer or inner boundary — toggles material,
        // so on an even crossing count the arcs ALTERNATE: material is
        // one parity class. The nearest-boundary side test is a local
        // vote that can misread an arc whose nearest boundary sample
        // belongs to a far strand (measured: a two-start thread's
        // runout band read T,T,T,F across four alternating arcs, the
        // widest "material" arc sat in a cut window, the kept shank
        // band classified Inside the tool and the body opened 665
        // mesh edges). Decide the class by majority over all arcs;
        // per-arc votes remain the rule when the count is odd
        // (tangential crossings merged by the dedup) or tied.
        let material_class = if n >= 2 && n.is_multiple_of(2) {
            let votes = |class: usize| {
                arcs.iter()
                    .enumerate()
                    .filter(|(i, arc)| arc.3 == Some(i % 2 == class))
                    .count()
            };
            let (even, odd) = (votes(0), votes(1));
            (even != odd).then_some(if even > odd { 0 } else { 1 })
        } else {
            None
        };
        for (i, &(width, mid, inside_inner, side)) in arcs.iter().enumerate() {
            let material = match material_class {
                Some(class) => i % 2 == class,
                None => side == Some(true),
            };
            if !inside_inner && material && best.map(|(bw, _)| width > bw).unwrap_or(true) {
                best = Some((width, mid));
            }
        }
        let (_, mid_u) = best?;
        if popkern_core::trace::trace_enabled("face-selection") {
            eprintln!(
                "      [selpoly] periodic-fold scan cv={cv:.3} xs={xs:?} -> mid_u={mid_u:.3}"
            );
        }
        return Some(surface.point_at(mid_u, cv));
    }
    // Candidate intervals inside the outer polygon (even-odd pairs).
    let mut best: Option<(f64, f64)> = None; // (width, mid_u)
    for pair in outer_x.chunks(2) {
        let [a, b] = pair else { continue };
        // Subtract inner coverage: split [a,b] at inner crossings.
        let mut cuts = vec![*a, *b];
        for inner in &inners {
            for x in crossings_at(inner, cv) {
                if x > *a && x < *b {
                    cuts.push(x);
                }
            }
        }
        cuts.sort_by(|p, q| p.partial_cmp(q).unwrap_or(std::cmp::Ordering::Equal));
        for w in cuts.windows(2) {
            let mid = 0.5 * (w[0] + w[1]);
            // Point-in-polygon (even-odd) for the mid against inners.
            let inside_inner = inners.iter().any(|inner| {
                let xs = crossings_at(inner, cv);
                xs.iter().filter(|&&x| x < mid).count() % 2 == 1
            });
            if !inside_inner {
                let width = w[1] - w[0];
                if best.map(|(bw, _)| width > bw).unwrap_or(true) {
                    best = Some((width, mid));
                }
            }
        }
    }
    let (_, mid_u) = best?;
    if popkern_core::trace::trace_enabled("face-selection") {
        let inner_x: Vec<Vec<f64>> = inners.iter().map(|i| crossings_at(i, cv)).collect();
        eprintln!(
            "      [selpoly] scan cv={cv:.3} outer_x={outer_x:?} inner_x={inner_x:?} -> mid_u={mid_u:.3}"
        );
    }
    Some(surface.point_at(mid_u, cv))
}

/// Material-true OUTWARD normal for a PLANAR face: the stored
/// outer-loop Newell vector (which points INTO the material in this
/// kernel's winding convention), negated. Returns `None` for
/// non-planar faces — coorientation classification then falls back to
/// the sign-arbitrary chart normal with orientation unknown.
pub(super) fn planar_outward_normal(
    brep: &BRep,
    face_handle: popkern_core::topology::FaceHandle,
) -> Option<popkern_core::geometry::Vec3> {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let face = brep.faces.get(face_handle)?;
    let surface = face.surface.as_ref()?;
    if !matches!(surface.geometry, SurfaceGeometry::Plane { .. }) {
        return None;
    }
    let pts = brep.face_outer_loop_positions(face_handle);
    if pts.len() < 3 {
        return None;
    }
    let mut n = popkern_core::geometry::Vec3::new(0.0, 0.0, 0.0);
    for i in 0..pts.len() {
        let p = pts[i];
        let q = pts[(i + 1) % pts.len()];
        n.x += (p.y - q.y) * (p.z + q.z);
        n.y += (p.z - q.z) * (p.x + q.x);
        n.z += (p.x - q.x) * (p.y + q.y);
    }
    (n * -1.0).normalized()
}
