//! Surface identity checks used by exact-boundary support-ribbon recovery.

use popkern_core::geometry::surface::SurfaceGeometry;

pub(super) fn geometrically_unique_surface(
    candidates: &[SurfaceGeometry],
    tolerance: f64,
) -> Option<SurfaceGeometry> {
    let first = candidates.first()?;
    candidates
        .iter()
        .all(|candidate| equivalent_surface(first, candidate, tolerance))
        .then(|| first.clone())
}

fn equivalent_surface(a: &SurfaceGeometry, b: &SurfaceGeometry, tolerance: f64) -> bool {
    match (a, b) {
        (
            SurfaceGeometry::Plane {
                origin: origin_a,
                normal: normal_a,
                ..
            },
            SurfaceGeometry::Plane {
                origin: origin_b,
                normal: normal_b,
                ..
            },
        ) => {
            let len_a = normal_a.length();
            let len_b = normal_b.length();
            if len_a <= 1.0e-12 || len_b <= 1.0e-12 {
                return false;
            }
            let parallel = (normal_a.dot(normal_b) / (len_a * len_b)).abs() >= 1.0 - 1.0e-10;
            let offset = ((*origin_b - *origin_a).dot(normal_a) / len_a).abs();
            parallel && offset <= tolerance
        }
        _ => false,
    }
}
