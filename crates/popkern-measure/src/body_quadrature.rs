//! Divergence-theorem volume by quadrature over the tessellation's
//! parameter triangles (PK-807).
//!
//! `V = ⅓ ∬ r · n dA`. A flat facet integrates `r · n` over the chord
//! surface, so a mesh's volume is off by the facet sagitta over the
//! whole area (measured against OpenCASCADE: a propeller −0.53%, a keyed
//! shaft −0.13%). Here each triangle of a face mesh integrates the
//! exact surface instead: over its (u, v) triangle,
//! `∬ r(u,v) · (∂u × ∂v) du dv` with a degree-5 seven-point rule. The
//! integrand carries the chart's Jacobian, so the triangle's signed
//! (u, v) area orients the term exactly as the mesh's own winding (the
//! mesher orients the whole mesh outward). What remains is the region
//! between each boundary chord in (u, v) and the true trim curve; a
//! plane takes its exact term from its edge curves.
//!
//! A triangle whose (u, v) corners do not reproduce its vertices (a
//! seam-branch or pole corner) keeps its flat facet term.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_mesh::TriMesh;

/// Dunavant degree-5 rule: barycentric points and weights (sum 1).
const RULE: [([f64; 3], f64); 7] = [
    ([1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0], 0.225),
    (
        [
            0.059_715_871_789_77,
            0.470_142_064_105_115,
            0.470_142_064_105_115,
        ],
        0.132_394_152_788_506,
    ),
    (
        [
            0.470_142_064_105_115,
            0.059_715_871_789_77,
            0.470_142_064_105_115,
        ],
        0.132_394_152_788_506,
    ),
    (
        [
            0.470_142_064_105_115,
            0.470_142_064_105_115,
            0.059_715_871_789_77,
        ],
        0.132_394_152_788_506,
    ),
    (
        [
            0.797_426_985_353_087,
            0.101_286_507_323_456,
            0.101_286_507_323_456,
        ],
        0.125_939_180_544_827,
    ),
    (
        [
            0.101_286_507_323_456,
            0.797_426_985_353_087,
            0.101_286_507_323_456,
        ],
        0.125_939_180_544_827,
    ),
    (
        [
            0.101_286_507_323_456,
            0.101_286_507_323_456,
            0.797_426_985_353_087,
        ],
        0.125_939_180_544_827,
    ),
];

fn as_vec(p: Point3) -> Vec3 {
    Vec3::new(p.x, p.y, p.z)
}

/// A corner whose one chart partial is under this fraction of the other
/// is a pole of the chart (see `surface_term`).
const POLE_RATIO: f64 = 1e-4;

/// Flat facet term `⅙ a · (b × c)`.
fn facet_term(a: Point3, b: Point3, c: Point3) -> f64 {
    as_vec(a).dot(&as_vec(b).cross(&as_vec(c))) / 6.0
}

/// `⅓ ∬ r · (∂u × ∂v)` over the oriented (u, v) triangle `uv`, or `None`
/// when its corners do not map onto `points` within `gate`.
fn surface_term(
    surface: &SurfaceGeometry,
    uv: [(f64, f64); 3],
    points: [Point3; 3],
    gate: f64,
) -> Option<f64> {
    for (q, p) in uv.iter().zip(points) {
        if surface.point_at(q.0, q.1).distance_to(&p) > gate {
            return None;
        }
    }
    // A corner at a chart pole (one partial vanishes: a sphere's pole,
    // a cone's apex) carries an arbitrary value of the collapsed
    // coordinate, so the (u, v) triangle covers only part of the strip
    // its facet spans (measured: a unit sphere read -3.3e-4). The pole
    // is really an edge: integrate the quad from the pole row under the
    // other two corners' collapsed coordinates.
    for k in 0..3 {
        let (u, v) = uv[k];
        let (b, c) = (uv[(k + 1) % 3], uv[(k + 2) % 3]);
        // Collapsed by ratio: a pole vertex may sit a few micro-radians
        // off the pole in parameter (measured 3e-6 on a primitive
        // sphere), where the vanishing partial is tiny but not zero.
        let (su, sv) = (
            surface.partial_u(u, v).length(),
            surface.partial_v(u, v).length(),
        );
        let collapsed_u = su <= POLE_RATIO * sv;
        let collapsed_v = sv <= POLE_RATIO * su;
        if collapsed_u == collapsed_v {
            continue;
        }
        let (pb, pc) = if collapsed_u {
            ((b.0, v), (c.0, v))
        } else {
            ((u, b.1), (u, c.1))
        };
        return Some(oriented_term(surface, [pb, b, c]) + oriented_term(surface, [pb, c, pc]));
    }
    Some(oriented_term(surface, uv))
}

/// `⅓ ∬ r · (∂u × ∂v)` over the oriented (u, v) triangle `uv`.
fn oriented_term(surface: &SurfaceGeometry, uv: [(f64, f64); 3]) -> f64 {
    let signed_area = 0.5
        * ((uv[1].0 - uv[0].0) * (uv[2].1 - uv[0].1) - (uv[2].0 - uv[0].0) * (uv[1].1 - uv[0].1));
    let mut sum = 0.0;
    for (bary, w) in RULE {
        let u = bary[0] * uv[0].0 + bary[1] * uv[1].0 + bary[2] * uv[2].0;
        let v = bary[0] * uv[0].1 + bary[1] * uv[1].1 + bary[2] * uv[2].1;
        let r = as_vec(surface.point_at(u, v));
        let n = surface.partial_u(u, v).cross(&surface.partial_v(u, v));
        sum += w * r.dot(&n);
    }
    sum * signed_area / 3.0
}

/// Signed volume of `mesh` (an outward-oriented tessellation of `brep`)
/// with every face's triangles integrated on the exact surface.
pub(crate) fn quadrature_volume(
    brep: &popkern_core::topology::BRep,
    mesh: &TriMesh,
    gate: f64,
) -> f64 {
    let mut volume = 0.0;
    for fm in &mesh.face_meshes {
        let surface = brep
            .faces
            .get(fm.face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry);
        if matches!(surface, Some(SurfaceGeometry::Plane { .. })) {
            // A plane's term is exact from its edge curves (circular
            // holes and arcs included, PK-771); its facets only
            // supply the mesh's outward sense.
            let facets: f64 = fm
                .triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(|i| fm.vertices[i].position);
                    facet_term(a, b, c)
                })
                .sum();
            let exact = crate::planar_area::planar_face_divergence(brep, fm.face);
            volume += if (exact - facets).abs() <= (exact + facets).abs() {
                exact
            } else {
                -exact
            };
            continue;
        }
        for t in &fm.triangles {
            let [a, b, c] = t.map(|i| &fm.vertices[i]);
            let points = [a.position, b.position, c.position];
            let flat = facet_term(points[0], points[1], points[2]);
            volume += surface
                .and_then(|s| {
                    // Corners on one branch of a periodic chart: a facet
                    // across the seam keeps u = 0.1 and u = 6.2 as stored,
                    // a (u, v) triangle spanning the whole period that
                    // integrated a band of the surface (measured: an arm
                    // part read 3.86x its volume). The anchor is the
                    // corner with the largest chart Jacobian: a pole
                    // corner's collapsed coordinate is arbitrary.
                    let (pu, pv) = s.parameter_periodicity();
                    let jacobian = |q: (f64, f64)| {
                        s.partial_u(q.0, q.1).cross(&s.partial_v(q.0, q.1)).length()
                    };
                    let corners = [a.uv, b.uv, c.uv];
                    let anchor = corners
                        .iter()
                        .copied()
                        .max_by(|x, y| jacobian(*x).total_cmp(&jacobian(*y)))
                        .unwrap_or(a.uv);
                    let near = |x: f64, to: f64, period: Option<f64>| match period {
                        Some(p) if p > 0.0 => x - p * ((x - to) / p).round(),
                        _ => x,
                    };
                    let [a, b, c] =
                        corners.map(|q| (near(q.0, anchor.0, pu), near(q.1, anchor.1, pv)));
                    let exact =
                        surface_term(s, [a, b, c], points, gate).filter(|x| x.is_finite())?;
                    // The (u, v) triangle must be this facet's patch: its
                    // centroid's image lies within the facet's diameter
                    // of the facet's centroid.
                    let diameter = points[0]
                        .distance_to(&points[1])
                        .max(points[1].distance_to(&points[2]))
                        .max(points[2].distance_to(&points[0]));
                    let image = s.point_at((a.0 + b.0 + c.0) / 3.0, (a.1 + b.1 + c.1) / 3.0);
                    let centroid = Point3::new(
                        (points[0].x + points[1].x + points[2].x) / 3.0,
                        (points[0].y + points[1].y + points[2].y) / 3.0,
                        (points[0].z + points[1].z + points[2].z) / 3.0,
                    );
                    if image.distance_to(&centroid) > diameter + gate {
                        return None;
                    }
                    // The (u, v) orientation must agree with the facet's
                    // (a mirrored chart triangle flips the sign).
                    let normal = (points[1] - points[0]).cross(&(points[2] - points[0]));
                    let centroid_uv = ((a.0 + b.0 + c.0) / 3.0, (a.1 + b.1 + c.1) / 3.0);
                    let chart = s
                        .partial_u(centroid_uv.0, centroid_uv.1)
                        .cross(&s.partial_v(centroid_uv.0, centroid_uv.1));
                    let signed_area = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
                    let agree = normal.dot(&chart) * signed_area;
                    (agree > 0.0).then_some(exact)
                })
                .unwrap_or(flat);
        }
    }
    volume
}
