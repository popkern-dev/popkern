//! Exact divergence terms of sphere, cylinder and cone faces whose
//! measuring mesh misread them (PK-1186).
//!
//! A face's term `⅓ ∬ r · n dA` needs the mesh only where `r · n` varies
//! over the surface in a way no boundary integral states. On a sphere
//! (centre `c`, radius `R`) and a cylinder (`c` on its axis) the point's
//! offset from `c` across the axis is `±R n`, so `r · n = c · n ± R`; on
//! a cone every point is the apex plus a generator normal to `n`, so
//! `r · n = apex · n`. The term is then `⅓ (c · A ± R a)`, with `a` the
//! face's area integrated over its stored pcurves and `A = ∬ n dA` its
//! vector area, half its loops' moments on the edge curves (Stokes).
//! The measuring mesh supplies only two signs, each the sign of a sum
//! whose every facet contributes its own area with that sign: whether
//! the folded normal the loops run about points out of the mesh, and
//! whether the outward normal points away from the axis.
//!
//! The term carries the trims' fit (a tenth of the tolerance: two
//! balls' rims moved it 1.5e-6 between two paths whose meshes agreed to
//! 2e-7), so it stands in only where the mesh misread the face: where
//! the mesh's area and the pcurves' disagree past the pcurves' band, the
//! perimeter times the linear tolerance. Beside a sphere's pole a ring
//! passing 1e-5 from it fell to the fallback grid and read the face's
//! area 1.93 over, the body's volume 4.8e-2 low; another seam turn read
//! it 1.04e-4 over.

use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::BRep;
use popkern_mesh::FaceMesh;

/// The exact term of a sphere, cylinder or cone face whose mesh misread
/// it (`mesh_area` the mesh's quadrature of its area), `None` where the
/// mesh read it (its term stands), for other surfaces, or where the
/// loops, area or mesh signs do not read.
pub(super) fn misread_face_term(
    brep: &BRep,
    fm: &FaceMesh,
    sof: &SurfaceOnFace,
    mesh_area: &dyn Fn() -> f64,
) -> Option<f64> {
    // The point whose offset is radial, the axis it is radial to (none:
    // a sphere's centre), and the radius (none: a cone).
    let (c, axis, radius) = match &sof.geometry {
        SurfaceGeometry::Sphere { center, radius, .. } => (*center, None, Some(*radius)),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => (*origin, Some(*axis), Some(*radius)),
        SurfaceGeometry::Cone { apex, .. } => (*apex, None, None),
        _ => return None,
    };
    let area = misread_area(brep, fm, mesh_area)?;
    let loops = crate::planar_area::face_vector_area(brep, fm.face, c)?;
    // Σ ± facet area: the folded normal against the facet's (outward),
    // the radial direction against the facet's.
    let (mut out_sense, mut away, mut facets) = (0.0, 0.0, 0.0);
    for t in &fm.triangles {
        let [p, q, r] = t.map(|i| fm.vertices[i].position);
        let n2 = (q - p).cross(&(r - p)) * 0.5;
        let facet = n2.length();
        let Some(n) = n2.normalized() else {
            continue;
        };
        facets += facet;
        let corner = fm.vertices[t[0]].uv;
        if let Some(folded) = sof.folded_normal(corner.0, corner.1) {
            out_sense += facet * folded.dot(&n).signum();
        }
        let centroid = Point3::new(
            (p.x + q.x + r.x) / 3.0,
            (p.y + q.y + r.y) / 3.0,
            (p.z + q.z + r.z) / 3.0,
        );
        let offset = centroid - c;
        let radial = match axis {
            Some(a) => offset - a * offset.dot(&a),
            None => offset,
        };
        if let Some(radial) = radial.normalized() {
            away += facet * radial.dot(&n).signum();
        }
    }
    // Each sign must hold over most of the mesh's area.
    let decided = |s: f64| (s.abs() > 0.5 * facets).then(|| s.signum());
    let outward: Vec3 = loops * decided(out_sense)?;
    let radial = match radius {
        Some(r) => decided(away)? * r * area,
        None => 0.0,
    };
    Some((c.x * outward.x + c.y * outward.y + c.z * outward.z + radial) / 3.0)
}

/// The face's area over its pcurves where the mesh's quadrature of it
/// disagrees past the pcurves' band (the perimeter times the linear
/// tolerance), `None` where the mesh read it.
fn misread_area(brep: &BRep, fm: &FaceMesh, mesh_area: &dyn Fn() -> f64) -> Option<f64> {
    let trims = crate::face::face_area(brep, fm.face, &brep.tolerance).ok()?;
    let mut uses = std::collections::HashMap::<(usize, usize), u32>::new();
    for t in &fm.triangles {
        for k in 0..3 {
            let (i, j) = (t[k], t[(k + 1) % 3]);
            *uses.entry((i.min(j), i.max(j))).or_default() += 1;
        }
    }
    let perimeter: f64 = uses
        .iter()
        .filter(|(_, n)| **n == 1)
        .map(|(&(i, j), _)| {
            fm.vertices[i]
                .position
                .distance_to(&fm.vertices[j].position)
        })
        .sum();
    ((mesh_area().abs() - trims).abs() > perimeter * brep.tolerance.linear()).then_some(trims)
}
