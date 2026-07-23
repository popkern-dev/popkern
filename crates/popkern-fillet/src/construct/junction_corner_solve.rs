//! Corner-cap geometry derivation from a valence-2 junction's
//! context (PK-141): shared support by structural identity, corner-
//! side trim ends, end stations, then the rolling-ball corner solve.

use popkern_core::geometry::Point3;

use super::EdgeGeometry;

/// Solve the corner-cap geometry for a true-corner junction: find
/// the SHARED support (present in both fillets' support pairs, by
/// structural identity), each fillet's own support and corner-side
/// trim endpoints, reconstruct the end stations, and run the
/// rolling-ball corner solve.
pub(super) fn solve_corner_cap_for_junction(
    edge_geometries: &[EdgeGeometry],
    a: usize,
    b: usize,
    v_pos: Point3,
    radius: f64,
) -> Option<super::corner_cap_geometry::CornerCapGeometry> {
    use popkern_core::geometry::curve::CurveEval;
    let eg_a = edge_geometries.get(a)?;
    let eg_b = edge_geometries.get(b)?;
    let sig = |g: &popkern_core::geometry::surface::SurfaceGeometry| format!("{g:?}");
    let pairs_a = [&eg_a.analysis.surface_1, &eg_a.analysis.surface_2];
    let pairs_b = [&eg_b.analysis.surface_1, &eg_b.analysis.surface_2];
    let (shared, own_a, own_b) = pairs_a.iter().find_map(|&sa| {
        pairs_b.iter().find_map(|&sb| {
            (sig(sa) == sig(sb)).then(|| {
                let own_a = if sig(pairs_a[0]) == sig(sa) {
                    pairs_a[1]
                } else {
                    pairs_a[0]
                };
                let own_b = if sig(pairs_b[0]) == sig(sb) {
                    pairs_b[1]
                } else {
                    pairs_b[0]
                };
                (sa, own_a, own_b)
            })
        })
    })?;
    // Corner-side trim endpoints per fillet: the endpoint of each
    // trim nearest the junction vertex.
    let near_end = |trim: &popkern_core::geometry::NurbsCurve| {
        let (t0, t1) = trim.parameter_range();
        let p0 = trim.point_at(t0);
        let p1 = trim.point_at(t1);
        if v_pos.distance_to(&p0) <= v_pos.distance_to(&p1) {
            p0
        } else {
            p1
        }
    };
    let tol = 1e-9;
    let end_a = super::corner_cap_geometry::end_station_center(
        &eg_a.analysis.surface_1,
        &eg_a.analysis.surface_2,
        near_end(&eg_a.trim_3d_1),
        near_end(&eg_a.trim_3d_2),
        radius,
        tol,
    )?;
    let end_b = super::corner_cap_geometry::end_station_center(
        &eg_b.analysis.surface_1,
        &eg_b.analysis.surface_2,
        near_end(&eg_b.trim_3d_1),
        near_end(&eg_b.trim_3d_2),
        radius,
        tol,
    )?;
    super::corner_cap_geometry::solve_corner_cap(
        shared,
        [own_a, own_b],
        [end_a, end_b],
        radius,
        tol,
    )
}
