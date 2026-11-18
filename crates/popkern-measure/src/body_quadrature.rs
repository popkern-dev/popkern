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

use popkern_core::geometry::chart::ChartDescriptor;
use popkern_core::geometry::chart::ops::aligned;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_mesh::TriMesh;

#[path = "body_quadrature_analytic.rs"]
mod analytic;
#[path = "body_quadrature_rule.rs"]
mod rule;
#[path = "body_quadrature_singular.rs"]
mod singular;
#[path = "body_quadrature_sliver.rs"]
mod sliver;
pub(super) use rule::oriented_term;
use singular::{FaceRows, face_rows, singular_corner};

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

/// Flat facet term `⅙ a · (b × c)`.
fn facet_term(a: Point3, b: Point3, c: Point3) -> f64 {
    as_vec(a).dot(&as_vec(b).cross(&as_vec(c))) / 6.0
}

/// `⅓ ∬ r · (∂u × ∂v)` over the oriented (u, v) triangle `uv`, or `None`
/// when its corners do not map onto `points` within `gate`.
fn surface_term(
    surface: &SurfaceGeometry,
    descriptor: &ChartDescriptor,
    uv: [(f64, f64); 3],
    points: [Point3; 3],
    (gate, vertices): (f64, &FaceRows),
    kind: Integrand,
) -> Option<f64> {
    // A corner at a chart pole (one partial vanishes: a sphere's pole,
    // a cone's apex) carries an arbitrary value of the collapsed
    // coordinate, so the (u, v) triangle covers only part of the strip
    // its facet spans (measured: a unit sphere read -3.3e-4). The pole
    // is really an edge: integrate the quad from the pole row (at its
    // own parameter) under the other two corners' collapsed coordinates.
    // That corner's stored parameters are not read, so they need not
    // reproduce it: a kilometre cone's apex vertex sat at v = 0.0018,
    // 2e-3 off its point, past the gate, and its whole fan fell to flat
    // facets (the volume 6.3e-6 low at a micron tolerance).
    let corner = singular_corner(descriptor, points, (vertices, gate));
    for (k, (q, p)) in uv.iter().zip(points).enumerate() {
        if corner.is_some_and(|c| c.k == k) {
            continue;
        }
        if surface.point_at(q.0, q.1).distance_to(&p) > gate {
            return None;
        }
    }
    if let Some(c) = corner {
        let (b, cc) = (uv[(c.k + 1) % 3], uv[(c.k + 2) % 3]);
        let (pb, pc) = if c.collapsed_u {
            ((b.0, c.at), (cc.0, c.at))
        } else {
            ((c.at, b.1), (c.at, cc.1))
        };
        let c = cc;
        return Some(
            oriented_term(surface, [pb, b, c], kind) + oriented_term(surface, [pb, c, pc], kind),
        );
    }
    Some(oriented_term(surface, uv, kind))
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
        // A sphere, cylinder or cone face the mesh misread takes its
        // exact term (`body_quadrature_analytic`).
        let exact = brep
            .faces
            .get(fm.face)
            .and_then(|f| f.surface.as_ref())
            .and_then(|sof| {
                let mesh_area = || face_quadrature(brep, fm, surface, gate, Integrand::Area);
                analytic::misread_face_term(brep, fm, sof, &mesh_area)
            });
        let term =
            exact.unwrap_or_else(|| face_quadrature(brep, fm, surface, gate, Integrand::Volume));
        out.push((fm.face, term));
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
    let mut uses = popkern_math::hash::FastHashMap::<(usize, usize), u32>::default();
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
    // The face's chart artifacts, from the chart service.
    let descriptor = surface.map(|s| ChartDescriptor::of(s, &brep.tolerance));
    // The face's vertices (where its loops meet a pole) and each pole's
    // resolution.
    let vertices = face_rows(brep, fm.face, descriptor.as_deref());
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
            .zip(descriptor.as_deref())
            .and_then(|(s, descriptor)| {
                // Corners on one branch of a periodic chart: a facet
                // across the seam keeps u = 0.1 and u = 6.2 as stored,
                // a (u, v) triangle spanning the whole period that
                // integrated a band of the surface (measured: an arm
                // part read 3.86x its volume). The anchor is the
                // corner with the largest chart Jacobian: a pole
                // corner's collapsed coordinate is arbitrary.
                let jacobian =
                    |q: (f64, f64)| s.partial_u(q.0, q.1).cross(&s.partial_v(q.0, q.1)).length();
                let corners = [a.uv, b.uv, c.uv];
                let anchor = corners
                    .iter()
                    .copied()
                    .max_by(|x, y| jacobian(*x).total_cmp(&jacobian(*y)))
                    .unwrap_or(a.uv);
                let [a, b, c] =
                    into_domain(descriptor, corners.map(|q| aligned(descriptor, q, anchor)));
                let exact = surface_term(s, descriptor, [a, b, c], points, (gate, &vertices), kind)
                    .filter(|x| x.is_finite())?;
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
                let corners = [a, b, c];
                let at_row = singular_corner(descriptor, points, (&vertices, gate));
                // A pole corner's stored collapsed coordinate is arbitrary:
                // the orientation is the integrated quad's, from the row
                // under the other two corners (a ball's own pole vertex
                // at the longitude it was made with read the triangle
                // mirrored, and its facet stood: two bodies one volume
                // read 2.3e-6 apart).
                let oriented = match at_row {
                    Some(r) => {
                        let (q, w) = (corners[(r.k + 1) % 3], corners[(r.k + 2) % 3]);
                        let (pq, pw) = if r.collapsed_u {
                            ((q.0, r.at), (w.0, r.at))
                        } else {
                            ((r.at, q.1), (r.at, w.1))
                        };
                        let twice = |x: (f64, f64), y: (f64, f64), z: (f64, f64)| {
                            (y.0 - x.0) * (z.1 - x.1) - (z.0 - x.0) * (y.1 - x.1)
                        };
                        twice(pq, q, w) + twice(pq, w, pw)
                    }
                    None => signed_area,
                };
                if oriented * uv_sign <= 0.0 && !collinear_uv([a, b, c], signed_area) {
                    return None;
                }
                // A chord from a pole corner runs in (u, v) from the pole
                // row under its other end: read at the pole's stored
                // coordinate, a turned cone's apex half a period from its
                // seam made the seam sliver a wedge of the whole fan (a
                // kilometre ice-cream cone read 4/3 of its area,
                // PK-1133).
                let chord = |k: usize| {
                    let (mut p, mut q) = (corners[k], corners[(k + 1) % 3]);
                    // The row corner at the row's own parameter, under its
                    // chord partner's collapsed coordinate.
                    match at_row {
                        Some(r) if r.k == k && r.collapsed_u => p = (q.0, r.at),
                        Some(r) if r.k == k => p = (r.at, q.1),
                        Some(r) if r.k == (k + 1) % 3 && r.collapsed_u => q = (p.0, r.at),
                        Some(r) if r.k == (k + 1) % 3 => q = (r.at, p.1),
                        _ => {}
                    }
                    (p, q)
                };
                let slivers: f64 = (0..3)
                    .filter(|&k| boundary[k] && !trims.is_empty())
                    .filter_map(|k| {
                        let (p, q) = chord(k);
                        sliver::sliver_term(s, descriptor, &trims, p, q, gate, kind)
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

/// A triangle aligned across a geometric closure's cut moved by whole
/// periods into the domain: a clamped closed NURBS evaluates past its
/// domain by extrapolation, and a triangle anchored on a corner at u = 0
/// put its neighbour at 0.992 on -0.008 (an atlas-meshed NURBS ball's
/// seam triangles fell back to their facets, 1.3e-5 of its volume low).
/// Periodic directions keep the anchor's branch.
fn into_domain(d: &ChartDescriptor, uv: [(f64, f64); 3]) -> [(f64, f64); 3] {
    use popkern_core::geometry::chart::{ChartPeriod, PeriodKind};
    let ((u0, u1), (v0, v1)) = d.domain;
    let shift = |period: Option<&ChartPeriod>, (lo, hi): (f64, f64), at: [f64; 3]| match period {
        Some(c) if c.kind == PeriodKind::GeometricClosure => {
            let (low, high) = (at[0].min(at[1]).min(at[2]), at[0].max(at[1]).max(at[2]));
            if low < lo {
                c.period * ((lo - low) / c.period).ceil()
            } else if high > hi {
                -c.period * ((high - hi) / c.period).ceil()
            } else {
                0.0
            }
        }
        _ => 0.0,
    };
    let du = shift(d.u.as_ref(), (u0, u1), uv.map(|q| q.0));
    let dv = shift(d.v.as_ref(), (v0, v1), uv.map(|q| q.1));
    uv.map(|(u, v)| (u + du, v + dv))
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
