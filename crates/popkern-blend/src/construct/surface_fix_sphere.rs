//! The inscribed corner sphere for a cube-style vertex blend (split
//! from `surface_fix.rs` under the 500-line budget).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;

/// Try to construct a sphere surface for a vertex blend corner.
///
/// Finds the nearest vertex of the original solid. The sphere center
/// is at (vertex ± R) in the direction of each adjacent face's inward
/// normal, with radius R = average fillet radius.
pub(super) fn try_sphere_for_corner(
    positions: &[popkern_core::geometry::Point3],
    original: &BRep,
    fillet_surfaces: &[SurfaceGeometry],
) -> Option<popkern_core::geometry::surface::SurfaceOnFace> {
    use popkern_core::geometry::surface::SurfaceOnFace;

    // Fillet radii from torus minor / cylinder radii. For the
    // sphere to be semantically correct as a corner patch, all
    // incident fillet strips must share one radius — otherwise
    // the rolling ball isn't the same at each contact and a
    // Gregory patch (or analogous N-sided construction) is the
    // right answer. Bail out on mixed radii so `build_smooth_blend`
    // keeps whatever N-sided NURBS it produced.
    // A swept NURBS strip is not a constant circular section — it is
    // a variable-radius or conic-profile blend (PK-189), and the
    // rolling ball that would justify a corner sphere does not exist
    // for it. Promoting anyway shipped a Sphere cap whose boundary
    // arcs lay the conic shoulder gap off the strips (measured
    // 2.7e-3 CurveSurfaceDeviation on the rho=0.35 cube corner).
    if fillet_surfaces
        .iter()
        .any(|s| matches!(s, SurfaceGeometry::Nurbs(_)))
    {
        return None;
    }
    let radii: Vec<f64> = fillet_surfaces
        .iter()
        .filter_map(|s| match s {
            SurfaceGeometry::Torus { minor_radius, .. } => Some(*minor_radius),
            SurfaceGeometry::Cylinder { radius, .. } => Some(*radius),
            _ => None,
        })
        .collect();
    if radii.len() >= 2 {
        let (min_r, max_r) = radii.iter().fold((f64::INFINITY, 0.0_f64), |(lo, hi), &r| {
            (lo.min(r), hi.max(r))
        });
        if max_r - min_r > 1e-6 {
            return None;
        }
    }
    let avg_radius = radii.first().copied().unwrap_or(0.1);

    // Orthogonality guard (PK-221): the centre construction below
    // offsets PER GLOBAL AXIS from the corner, which is only
    // meaningful for a cube-style corner — mutually orthogonal,
    // axis-aligned fillet axes. A pyramid apex (non-orthogonal
    // supports) satisfied every other gate and shipped a sphere with
    // g1 ~1.0 at its rim. Decline so the constructed N-sided cap
    // stands.
    let axes: Vec<popkern_core::geometry::Vec3> = fillet_surfaces
        .iter()
        .filter_map(|s| match s {
            SurfaceGeometry::Torus { axis, .. } | SurfaceGeometry::Cylinder { axis, .. } => {
                axis.normalized()
            }
            _ => None,
        })
        .collect();
    const AXIS_ALIGNMENT_TOLERANCE: f64 = 1.0e-6;
    for axis in &axes {
        let max_component = axis.x.abs().max(axis.y.abs()).max(axis.z.abs());
        if max_component < 1.0 - AXIS_ALIGNMENT_TOLERANCE {
            return None;
        }
    }
    for (i, a) in axes.iter().enumerate() {
        for b in axes.iter().skip(i + 1) {
            if a.dot(b).abs() > AXIS_ALIGNMENT_TOLERANCE
                && a.dot(b).abs() < 1.0 - AXIS_ALIGNMENT_TOLERANCE
            {
                return None;
            }
        }
    }

    // Use the per-coordinate **median** position rather than the
    // mean centroid. The boundary of a multi-fillet vertex blend
    // can include outlier vertices that bleed in from adjacent
    // fillets' regions; the mean centroid is biased by them, but
    // the median sits inside the dominant cluster and is robust.
    let median_axis = |xs: &mut Vec<f64>| -> f64 {
        xs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        xs[xs.len() / 2]
    };
    let mut xs: Vec<f64> = positions.iter().map(|p| p.x).collect();
    let mut ys: Vec<f64> = positions.iter().map(|p| p.y).collect();
    let mut zs: Vec<f64> = positions.iter().map(|p| p.z).collect();
    let cx = median_axis(&mut xs);
    let cy = median_axis(&mut ys);
    let cz = median_axis(&mut zs);

    // Find nearest original vertex to the median position.
    let mut best_dist = f64::INFINITY;
    let mut best_pos = None;
    for (_, v) in &original.vertices {
        let d = ((v.position.x - cx).powi(2)
            + (v.position.y - cy).powi(2)
            + (v.position.z - cz).powi(2))
        .sqrt();
        if d < best_dist {
            best_dist = d;
            best_pos = Some(v.position);
        }
    }

    let corner = best_pos?;
    if best_dist > avg_radius * 3.0 {
        return None; // Too far from any corner.
    }

    // Sphere center: offset from corner toward the solid interior.
    // For each axis, move by R in the direction away from the nearest
    // face (toward the interior).
    //
    // **Direction-check via centroid, not median.** The median sits
    // inside the dominant cluster (robust to outliers) but on a
    // symmetric N-vertex cap (e.g. cube_uniform's 3 rolling-ball
    // footprints at distance r along each cube axis), the median
    // equals the corner position exactly — defeating the
    // `>`-comparison and forcing the wrong (`else`) branch, which
    // produced sphere centres OUTSIDE the solid (cube_uniform
    // shipped with centre = (-r, -r, -r) instead of (+r, +r, +r),
    // explaining the long-standing 5-15 cm cap-vs-cube drift).
    //
    // The centroid (mean) breaks the symmetric tie consistently:
    // for the same 3 footprints, centroid = (r/3, r/3, r/3),
    // strictly greater than the corner on every axis. We use the
    // median for `cx/cy/cz` (it's the right choice for *position*),
    // but the centroid for the *direction* of the offset.
    let mean_x = positions.iter().map(|p| p.x).sum::<f64>() / positions.len() as f64;
    let mean_y = positions.iter().map(|p| p.y).sum::<f64>() / positions.len() as f64;
    let mean_z = positions.iter().map(|p| p.z).sum::<f64>() / positions.len() as f64;
    let sx = if mean_x > corner.x {
        corner.x + avg_radius
    } else {
        corner.x - avg_radius
    };
    let sy = if mean_y > corner.y {
        corner.y + avg_radius
    } else {
        corner.y - avg_radius
    };
    let sz = if mean_z > corner.z {
        corner.z + avg_radius
    } else {
        corner.z - avg_radius
    };

    Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Sphere {
            center: popkern_core::geometry::Point3::new(sx, sy, sz),
            radius: avg_radius,
            axis: popkern_core::geometry::Vec3::Z,
        },
        orientation: true,
    })
}
