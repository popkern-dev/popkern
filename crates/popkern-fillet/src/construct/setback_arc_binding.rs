//! Strip-face binding and projection repair for the setback arcs
//! (extracted from `setback_arc_crossings` for the 500-line file
//! budget; same nested-module pattern).

use super::setback_arc_crossings::arc_inside_trim_fraction;
use super::*;

/// Face carrying the strip piece the setback arc lives on.
///
/// A strip can be split into several faces by pre-carve operations
/// (trim-crossing anchors, far-corner surgery), and the apex-end
/// piece may carry a RE-AUTHORED sub-range chart that no longer
/// matches the sweep surface's parameterization (measured on
/// hexagon6_variable: fillet 0's arc sits 1.7e-4 from apex piece
/// 22v1 while the parameterization-exact test only sees the remote
/// piece 8v1, whose trim the arc never crosses). Bind by the arc
/// itself: the face whose surface is nearest the arc's interior
/// samples, within the cubic-vs-exact-row sagitta band (~2e-4
/// measured; 1e-3 accepts the class with margin). Fall back to the
/// parameterization-exact probe (apex-nearest piece) when no face
/// hugs the arc.
pub(super) fn find_fillet_face(
    brep: &BRep,
    fillet: &FilletSurfaceResult,
    arc_cubic: &[Point3; 4],
    apex: Point3,
) -> Result<FaceHandle, String> {
    // Band: cubic-vs-exact-row sagitta (~2e-4) plus the measured
    // extractor wobble on variable-radius corners (corner points sit
    // up to ~3e-3 off the strips). Same-band ties (wall planes near
    // tangency) are resolved by trim containment then distance.
    const ARC_ON_SURFACE_BAND: f64 = 5.0e-3;
    let arc = bezier_cubic_curve(arc_cubic);
    let arc_eval: &dyn CurveEval = &popkern_core::geometry::CurveGeometry::Nurbs(arc);
    let samples: Vec<Point3> = [0.25, 0.5, 0.75]
        .iter()
        .map(|&t| arc_eval.point_at(t))
        .collect();
    let mut best_by_arc: Option<(FaceHandle, f64, f64)> = None;
    for (fh, face) in brep.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        // Cap-family skip (PK-301): fillet strips are authored as
        // NURBS sweeps or analytic surfaces — the Gregory and
        // polyhedral-spline variants exist only as corner CAP
        // geometry and can never be the strip this binder looks for.
        // Skipping them keeps this scan O(strips) instead of paying a
        // Newton projection against every accumulated cap piece
        // (measured: +14 ms per corner on the 12-edge G2 cube, the
        // dominant per-corner growth term).
        if matches!(
            surface.geometry,
            popkern_core::geometry::SurfaceGeometry::Gregory(_)
                | popkern_core::geometry::SurfaceGeometry::PolyhedralSpline(_)
        ) {
            continue;
        }
        let mut worst = 0.0_f64;
        for sample in &samples {
            let d = popkern_core::geometry::uv_projection::project_point_to_uv(
                &surface.geometry,
                sample,
            )
            .map(|(u, v)| {
                (&surface.geometry as &dyn SurfaceEval)
                    .point_at(u, v)
                    .distance_to(sample)
            })
            .unwrap_or(f64::INFINITY);
            worst = worst.max(d);
        }
        if worst >= ARC_ON_SURFACE_BAND {
            continue;
        }
        // Several faces can share the underlying surface (split strip
        // pieces) — the trim that actually CONTAINS the arc wins;
        // containment fraction first, surface distance as tie-break.
        let contained = arc_inside_trim_fraction(brep, fh, arc_cubic).unwrap_or(0.0);
        if best_by_arc.is_none_or(|(_, bc, bd)| contained > bc || (contained == bc && worst < bd)) {
            best_by_arc = Some((fh, contained, worst));
        }
    }
    if let Some((fh, contained, worst)) = best_by_arc
        && contained > 0.0
    {
        if popkern_core::trace::trace_enabled("corner-carve") {
            eprintln!(
                "    [fillet-face] arc-bound {fh:?} contained={contained:.2} surface_d={worst:.3e}"
            );
        }
        return Ok(fh);
    }

    let probes = [(0.5, 0.5), (0.25, 0.75), (0.75, 0.25)];
    let targets: Vec<Point3> = probes
        .iter()
        .map(|&(u, v)| fillet.surface.point_at(u, v))
        .collect();
    let mut best: Option<(FaceHandle, f64)> = None;
    for (fh, face) in brep.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let matches = probes.iter().zip(&targets).all(|(&(u, v), target)| {
            (&surface.geometry as &dyn SurfaceEval)
                .point_at(u, v)
                .distance_to(target)
                < 1.0e-9
        });
        if !matches {
            continue;
        }
        let d = face_centroid(brep, fh)
            .map(|c| c.distance_to(&apex))
            .unwrap_or(f64::INFINITY);
        if best.is_none_or(|(_, bd)| d < bd) {
            best = Some((fh, d));
        }
    }
    best.map(|(fh, _)| fh)
        .ok_or_else(|| "no face carries the fillet sweep surface".to_string())
}

/// Repair point-inversion samples that diverged (Gregory
/// `closest_point` Newton can escape the patch domain near seam
/// ends — measured uv reaching v = 1.45 with 1.4e-2..2.5e-2 3-D
/// error on the hexagon6_variable fan while neighbours project at
/// 1e-4 class). A sample is an outlier when its surface-evaluation
/// residual exceeds 10x the median residual (plus an absolute
/// floor); outliers are replaced by linear UV interpolation between
/// the nearest sound neighbours.
pub(super) fn repair_diverged_projections(
    cps: &mut [Point3],
    samples3d: &[Point3],
    surface: &popkern_core::geometry::SurfaceGeometry,
) {
    use popkern_core::geometry::surface::SurfaceEval;
    let n = cps.len();
    if n < 3 || samples3d.len() != n {
        return;
    }
    let residual = |uv: &Point3, p: &Point3| surface.point_at(uv.x, uv.y).distance_to(p);
    let residuals: Vec<f64> = cps
        .iter()
        .zip(samples3d)
        .map(|(uv, p)| residual(uv, p))
        .collect();
    let mut sorted = residuals.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let median = sorted[n / 2];
    let band = (10.0 * median).max(1e-9);
    let good: Vec<bool> = residuals.iter().map(|&r| r <= band).collect();
    if good.iter().all(|&g| g) || good.iter().filter(|&&g| g).count() < 2 {
        return;
    }
    for i in 0..n {
        if good[i] {
            continue;
        }
        let prev = (0..i).rev().find(|&j| good[j]);
        let next = ((i + 1)..n).find(|&j| good[j]);
        let repaired = match (prev, next) {
            (Some(a), Some(b)) => {
                let w = (i - a) as f64 / (b - a) as f64;
                Point3::new(
                    cps[a].x + w * (cps[b].x - cps[a].x),
                    cps[a].y + w * (cps[b].y - cps[a].y),
                    0.0,
                )
            }
            // Endpoint runs: extrapolate from the two nearest sound
            // samples.
            (Some(a), None) if a >= 1 && good[a - 1] => {
                let d = (i - a) as f64;
                Point3::new(
                    cps[a].x + d * (cps[a].x - cps[a - 1].x),
                    cps[a].y + d * (cps[a].y - cps[a - 1].y),
                    0.0,
                )
            }
            (None, Some(b)) if b + 1 < n && good[b + 1] => {
                let d = (b - i) as f64;
                Point3::new(
                    cps[b].x + d * (cps[b].x - cps[b + 1].x),
                    cps[b].y + d * (cps[b].y - cps[b + 1].y),
                    0.0,
                )
            }
            _ => continue,
        };
        // Accept only measured improvement.
        if residual(&repaired, &samples3d[i]) < residuals[i] {
            cps[i] = repaired;
        }
    }
}
