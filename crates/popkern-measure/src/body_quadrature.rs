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
//! seam-branch or pole corner) keeps its flat facet term. A triangle
//! integrated on its surface also integrates the sliver between each of
//! its boundary chords and the face's trim (`body_quadrature_sliver`).
//!
//! The same machinery measures a face's AREA (`∬ |∂u × ∂v| du dv`):
//! summed facets fall short of a curved face by their sagitta's share
//! (measured: an elliptic cylinder 4.5e-6, a paraboloid 5.8e-6), the
//! exact integrand does not.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_mesh::TriMesh;

#[path = "body_quadrature_sliver.rs"]
mod sliver;

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

/// What a quadrature integrates over the surface.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Integrand {
    /// `⅓ r · (∂u × ∂v)`: the divergence-theorem volume term, signed
    /// by the (u, v) winding.
    Volume,
    /// `|∂u × ∂v|`: area, signed by the (u, v) winding (the caller
    /// multiplies by the face's winding sign to make it positive).
    Area,
}

/// A corner whose one chart partial is under this fraction of the other
/// is a pole of the chart (see `surface_term`).
const POLE_RATIO: f64 = 1e-4;

/// A corner whose one chart partial is under this fraction of the same
/// partial at the triangle's other corners is at a pole (see
/// `pole_corner`).
const POLE_RELATIVE: f64 = 1e-3;

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
    kind: Integrand,
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
    if let Some((k, collapsed_u)) = pole_corner(surface, uv) {
        let (u, v) = uv[k];
        let (b, c) = (uv[(k + 1) % 3], uv[(k + 2) % 3]);
        let (pb, pc) = if collapsed_u {
            ((b.0, v), (c.0, v))
        } else {
            ((u, b.1), (u, c.1))
        };
        return Some(
            oriented_term(surface, [pb, b, c], kind) + oriented_term(surface, [pb, c, pc], kind),
        );
    }
    Some(oriented_term(surface, uv, kind))
}

/// The corner of `uv` at a chart pole, and whether its `u` (else its
/// `v`) is the collapsed coordinate. Collapsed by ratio to the other
/// partial: a pole vertex may sit a few micro-radians off the pole in
/// parameter (measured 3e-6 on a primitive sphere), where the vanishing
/// partial is tiny but not zero. Or relative to the triangle's other
/// corners: a boolean's apex vertex lies within its tolerance of the
/// apex, not at it, on a kilometre-scale cone 0.016 up the slant where
/// the partial ratio is 1e-2 but the rim's `∂u` is 1e6 times the apex's
/// (PK-1133; read as a triangle, the apex fan integrated one row of the
/// cone and lost a third of it). Both are ratios, so neither reads a
/// small ball's rows near its pole (0.5 to their neighbours) or a long
/// sliver on a cylinder (1) as a pole.
fn pole_corner(surface: &SurfaceGeometry, uv: [(f64, f64); 3]) -> Option<(usize, bool)> {
    let partials = uv.map(|(u, v)| {
        (
            surface.partial_u(u, v).length(),
            surface.partial_v(u, v).length(),
        )
    });
    (0..3).find_map(|k| {
        let (su, sv) = partials[k];
        let others = [partials[(k + 1) % 3], partials[(k + 2) % 3]];
        let su_far = others[0].0.max(others[1].0);
        let sv_far = others[0].1.max(others[1].1);
        let collapsed_u = su <= POLE_RATIO * sv || su <= POLE_RELATIVE * su_far;
        let collapsed_v = sv <= POLE_RATIO * su || sv <= POLE_RELATIVE * sv_far;
        (collapsed_u != collapsed_v).then_some((k, collapsed_u))
    })
}

/// `kind`'s integrand over the oriented (u, v) triangle `uv`.
pub(super) fn oriented_term(
    surface: &SurfaceGeometry,
    uv: [(f64, f64); 3],
    kind: Integrand,
) -> f64 {
    let signed_area = 0.5
        * ((uv[1].0 - uv[0].0) * (uv[2].1 - uv[0].1) - (uv[2].0 - uv[0].0) * (uv[1].1 - uv[0].1));
    let mut sum = 0.0;
    for (bary, w) in RULE {
        let u = bary[0] * uv[0].0 + bary[1] * uv[1].0 + bary[2] * uv[2].0;
        let v = bary[0] * uv[0].1 + bary[1] * uv[1].1 + bary[2] * uv[2].1;
        let n = surface.partial_u(u, v).cross(&surface.partial_v(u, v));
        sum += w * match kind {
            Integrand::Volume => as_vec(surface.point_at(u, v)).dot(&n) / 3.0,
            Integrand::Area => n.length(),
        };
    }
    sum * signed_area
}

/// The (u, v) winding sign of a face mesh's outward-wound triangles:
/// the sign of the chart normal against the facets' outward normals,
/// area-weighted over the whole face (a facet's cross product is twice
/// its area). `None` when the face does not decide (zero or non-finite):
/// its facets stand.
fn face_winding(surface: &SurfaceGeometry, fm: &popkern_mesh::FaceMesh) -> Option<f64> {
    let mut vote = 0.0;
    for t in &fm.triangles {
        let [a, b, c] = t.map(|i| fm.vertices[i].uv);
        let [p, q, r] = t.map(|i| fm.vertices[i].position);
        let centroid = ((a.0 + b.0 + c.0) / 3.0, (a.1 + b.1 + c.1) / 3.0);
        let chart = surface
            .partial_u(centroid.0, centroid.1)
            .cross(&surface.partial_v(centroid.0, centroid.1));
        vote += (q - p).cross(&(r - p)).dot(&chart);
    }
    (vote.is_finite() && vote != 0.0).then(|| vote.signum())
}

/// Signed volume of `mesh` (an outward-oriented tessellation of `brep`)
/// with every face's triangles integrated on the exact surface.
pub(crate) fn quadrature_volume(
    brep: &popkern_core::topology::BRep,
    mesh: &TriMesh,
    gate: f64,
) -> f64 {
    quadrature_face_volumes(brep, mesh, gate)
        .into_iter()
        .map(|(_, v)| v)
        .sum()
}

/// Each face's divergence term `⅓ ∬ r · n dA`, `n` its folded normal
/// (the mesh's outward sense): the volume behind the face less the
/// volume in front, summed over a region's boundary.
pub(crate) fn quadrature_face_volumes(
    brep: &popkern_core::topology::BRep,
    mesh: &TriMesh,
    gate: f64,
) -> Vec<(popkern_core::topology::FaceHandle, f64)> {
    let mut out = Vec::with_capacity(mesh.face_meshes.len());
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
            let term = if (exact - facets).abs() <= (exact + facets).abs() {
                exact
            } else {
                -exact
            };
            out.push((fm.face, term));
            continue;
        }
        out.push((
            fm.face,
            face_quadrature(brep, fm, surface, gate, Integrand::Volume),
        ));
    }
    out
}

/// `kind` over one face's mesh triangles, each on the exact `surface`
/// with its boundary slivers, or its flat facet where the (u, v)
/// triangle is not the facet's patch.
fn face_quadrature(
    brep: &popkern_core::topology::BRep,
    fm: &popkern_mesh::FaceMesh,
    surface: Option<&SurfaceGeometry>,
    gate: f64,
    kind: Integrand,
) -> f64 {
    let mut out = 0.0;
    // Chords on the face's boundary: mesh edges of one triangle.
    let mut uses = std::collections::HashMap::<(usize, usize), u32>::new();
    for t in &fm.triangles {
        for k in 0..3 {
            let (i, j) = (t[k], t[(k + 1) % 3]);
            *uses.entry((i.min(j), i.max(j))).or_default() += 1;
        }
    }
    let trims = surface
        .map(|_| sliver::face_trims(brep, fm.face))
        .unwrap_or_default();
    let winding = surface.and_then(|s| face_winding(s, fm));
    for t in &fm.triangles {
        let boundary = [0, 1, 2].map(|k| {
            let (i, j) = (t[k], t[(k + 1) % 3]);
            uses.get(&(i.min(j), i.max(j))) == Some(&1)
        });
        let [a, b, c] = t.map(|i| &fm.vertices[i]);
        let points = [a.position, b.position, c.position];
        let flat = match kind {
            Integrand::Volume => facet_term(points[0], points[1], points[2]),
            Integrand::Area => facet_area(points[0], points[1], points[2]),
        };
        out += surface
            .and_then(|s| {
                // Corners on one branch of a periodic chart: a facet
                // across the seam keeps u = 0.1 and u = 6.2 as stored,
                // a (u, v) triangle spanning the whole period that
                // integrated a band of the surface (measured: an arm
                // part read 3.86x its volume). The anchor is the
                // corner with the largest chart Jacobian: a pole
                // corner's collapsed coordinate is arbitrary.
                let (pu, pv) = s.parameter_periodicity();
                let jacobian =
                    |q: (f64, f64)| s.partial_u(q.0, q.1).cross(&s.partial_v(q.0, q.1)).length();
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
                let [a, b, c] = corners.map(|q| (near(q.0, anchor.0, pu), near(q.1, anchor.1, pv)));
                let exact =
                    surface_term(s, [a, b, c], points, gate, kind).filter(|x| x.is_finite())?;
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
                // The (u, v) orientation must be the face's: a
                // mirrored chart triangle flips the sign. Decided
                // against the face's own winding, not the facet's
                // normal: a sliver facet chording 27° of a helix
                // tips its normal across the surface's, and 50477 of
                // an M6 flank's triangles fell back to facets whose
                // errors no longer cancelled (0.7% of a tapped
                // block, PK-901).
                //
                // A triangle collinear in (u, v) covers no surface: its
                // term is its slivers alone. The mesher closes a curved
                // rim with such triangles (three rim vertices on one
                // chart row), whose facets span the arc's sagitta off
                // the surface: taken flat, two spherical sectors'
                // common read 1.6e-3 low with its area exact (PK-1138).
                let signed_area = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
                let uv_sign = winding?;
                if signed_area * uv_sign <= 0.0 && !collinear_uv([a, b, c], signed_area) {
                    return None;
                }
                let corners = [a, b, c];
                // A chord from a pole corner runs in (u, v) from the pole
                // row under its other end: read at the pole's stored
                // coordinate, a turned cone's apex half a period from its
                // seam made the seam sliver a wedge of the whole fan (a
                // kilometre ice-cream cone read 4/3 of its area,
                // PK-1133).
                let pole = pole_corner(s, corners);
                let chord = |k: usize| {
                    let (mut p, mut q) = (corners[k], corners[(k + 1) % 3]);
                    match pole {
                        Some((i, true)) if i == k => p.0 = q.0,
                        Some((i, false)) if i == k => p.1 = q.1,
                        Some((i, true)) if i == (k + 1) % 3 => q.0 = p.0,
                        Some((i, false)) if i == (k + 1) % 3 => q.1 = p.1,
                        _ => {}
                    }
                    (p, q)
                };
                let slivers: f64 = (0..3)
                    .filter(|&k| boundary[k] && !trims.is_empty())
                    .filter_map(|k| {
                        let (p, q) = chord(k);
                        sliver::sliver_term(s, &trims, p, q, gate, kind)
                    })
                    .sum();
                Some(match kind {
                    Integrand::Volume => exact + slivers,
                    // The winding makes a face's (u, v) triangles positive.
                    Integrand::Area => (exact + slivers) * uv_sign,
                })
            })
            .unwrap_or(flat);
    }
    out
}

/// Whether the (u, v) triangle `uv` of doubled signed area
/// `signed_area` is collinear up to rounding: its area under 64 ulps of
/// its longest edge's square.
fn collinear_uv(uv: [(f64, f64); 3], signed_area: f64) -> bool {
    let edge = |p: (f64, f64), q: (f64, f64)| (q.0 - p.0).powi(2) + (q.1 - p.1).powi(2);
    let longest = edge(uv[0], uv[1])
        .max(edge(uv[1], uv[2]))
        .max(edge(uv[2], uv[0]));
    let scale = uv
        .iter()
        .fold(1.0_f64, |m, q| m.max(q.0.abs()).max(q.1.abs()));
    signed_area.abs() <= 64.0 * f64::EPSILON * scale * longest.sqrt()
}

/// Area of the flat facet `a b c`.
fn facet_area(a: Point3, b: Point3, c: Point3) -> f64 {
    0.5 * (b - a).cross(&(c - a)).length()
}

/// Area of a non-planar face from its measuring mesh: every triangle
/// integrated on the exact surface (see the module notes).
pub(crate) fn quadrature_face_area(
    brep: &popkern_core::topology::BRep,
    fm: &popkern_mesh::FaceMesh,
    gate: f64,
) -> f64 {
    let surface = brep
        .faces
        .get(fm.face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry);
    face_quadrature(brep, fm, surface, gate, Integrand::Area)
}

#[cfg(test)]
#[path = "body_quadrature_tests.rs"]
mod tests;
