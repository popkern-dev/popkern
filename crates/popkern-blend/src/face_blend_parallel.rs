//! Line-contact (parallel-support) face–face blend path (PK-142):
//! antiparallel normals make the contact Jacobian singular and
//! tangency to both walls exists only at radius = gap/2, so the strip
//! is a CANONICAL CHOICE, not a solve — seed at wall A's outer-loop
//! centroid, march closed-form along the trim polygon's dominant
//! extent, and bow the half-circle sections toward the closed side of
//! the gap.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

use crate::contact::ContactSolution;
use crate::error::BlendError;
use crate::face_blend::{FaceBlendResult, FaceBlendSpec, FaceDomain, clamp_uv};

/// Closed-form march for (anti)parallel supports: centres on the
/// medial line, contacts as alternating closest-point feet; sections
/// are half-circles (line-contact family, radius = gap/2).
pub(crate) fn parallel_gap_blend(
    brep: &BRep,
    spec: &FaceBlendSpec,
    a: &FaceDomain,
    b: &FaceDomain,
    seed_uva: (f64, f64),
    seed_uvb: (f64, f64),
    tolerance: Tolerance,
) -> Result<FaceBlendResult, BlendError> {
    let radius = spec.radius;
    // Parallel supports admit a two-parameter family of tangent
    // lines, so the strip placement is a canonical CHOICE, not a
    // solve: seed at face A's outer-loop centroid (interior by
    // construction for the supported convex-ish trims) and march
    // along the trim polygon's dominant extent, so the strip spans
    // the wall's long dimension deterministically.
    let positions = brep.face_outer_loop_positions(spec.face_a);
    let centroid = loop_position_centroid(&positions);
    let seed_uva = a
        .surface
        .closest_point(&centroid)
        .map(|uv| clamp_uv(uv, a.uv_range))
        .unwrap_or(seed_uva);
    let seed_pa = a.surface.point_at(seed_uva.0, seed_uva.1);
    let seed_uvb = b
        .surface
        .closest_point(&seed_pa)
        .map(|uv| clamp_uv(uv, b.uv_range))
        .unwrap_or(seed_uvb);
    let na = a.surface.normal_at(seed_uva.0, seed_uva.1);
    let tangent = dominant_extent_direction(&positions, na)
        .or_else(|| {
            let du = a.surface.partial_u(seed_uva.0, seed_uva.1);
            (du - na * du.dot(&na)).normalized()
        })
        .ok_or_else(|| BlendError::InvalidParameter {
            parameter: "face_a".into(),
            value: "parallel gap".into(),
            reason: "degenerate support tangent".into(),
        })?;
    let sample = |uva: (f64, f64), uvb: (f64, f64)| -> ContactSolution {
        let pa = a.surface.point_at(uva.0, uva.1);
        let pb = b.surface.point_at(uvb.0, uvb.1);
        let centre = Point3::new(
            (pa.x + pb.x) * 0.5,
            (pa.y + pb.y) * 0.5,
            (pa.z + pb.z) * 0.5,
        );
        ContactSolution {
            uv0: uva,
            uv1: uvb,
            center: centre,
            contact0: pa,
            contact1: pb,
            iterations: 0,
            step_rejections: 0,
            residual: 0.0,
            condition: 1.0,
            jacobian: [[0.0; 4]; 4],
        }
    };
    let step = (radius * 0.5).max(tolerance.linear() * 100.0);
    let advance =
        |uva: (f64, f64), _uvb: (f64, f64), sign: f64| -> Option<((f64, f64), (f64, f64))> {
            let pa = a.surface.point_at(uva.0, uva.1);
            let target_a = pa + tangent * (step * sign);
            let next_a = clamp_uv(a.surface.closest_point(&target_a)?, a.uv_range);
            let pa2 = a.surface.point_at(next_a.0, next_a.1);
            let next_b = clamp_uv(b.surface.closest_point(&pa2)?, b.uv_range);
            Some((next_a, next_b))
        };
    let mut chain = vec![sample(seed_uva, seed_uvb)];
    for sign in [1.0, -1.0] {
        let (mut uva, mut uvb) = (seed_uva, seed_uvb);
        for _ in 0..400 {
            let Some((na2, nb2)) = advance(uva, uvb, sign) else {
                break;
            };
            let moved = a
                .surface
                .point_at(na2.0, na2.1)
                .distance_to(&a.surface.point_at(uva.0, uva.1));
            let sol = sample(na2, nb2);
            if sign > 0.0 {
                chain.push(sol);
            } else {
                chain.insert(0, sol);
            }
            uva = na2;
            uvb = nb2;
            if moved < step * 0.1 {
                break;
            }
        }
    }
    if chain.len() < 2 {
        return Err(BlendError::InvalidParameter {
            parameter: "radius".into(),
            value: format!("{radius}"),
            reason: "parallel-gap march produced fewer than two sections".into(),
        });
    }
    let resolver = |target: Point3, _warm: &ContactSolution| -> Option<ContactSolution> {
        let uva2 = clamp_uv(a.surface.closest_point(&target)?, a.uv_range);
        let pa2 = a.surface.point_at(uva2.0, uva2.1);
        let uvb2 = clamp_uv(b.surface.closest_point(&pa2)?, b.uv_range);
        Some(sample(uva2, uvb2))
    };
    let seed_pb = b.surface.point_at(seed_uvb.0, seed_uvb.1);
    let bulge = parallel_bulge_direction(brep, spec, seed_pa, seed_pb, tangent);
    crate::face_blend_attach::finish_face_blend(brep, spec, chain, &resolver, bulge, tolerance)
}

/// Bow side for the line-contact (parallel-support) family: the
/// section half-circles must bow toward the CLOSED side of the gap —
/// the side where a neighbor face of wall A lies inside the slab
/// between the two support planes (e.g. a slot floor). Returns `None`
/// when no neighbor disambiguates (a free-floating parallel pair),
/// which downstream degrades to the chord sheet and declines attach.
fn parallel_bulge_direction(
    brep: &BRep,
    spec: &FaceBlendSpec,
    seed_pa: Point3,
    seed_pb: Point3,
    tangent: Vec3,
) -> Option<Vec3> {
    let across = seed_pb - seed_pa;
    let gap = across.length();
    let na_dir = across.normalized()?;
    let wall = brep.faces.get(spec.face_a)?;
    let mut he = brep.loops[wall.outer_loop].half_edge;
    let start = he;
    for _ in 0..brep.half_edges.len() + 1 {
        let twin = brep.half_edges[he].twin;
        let neighbor = brep.loops[brep.half_edges[twin].loop_].face;
        if neighbor != spec.face_a && neighbor != spec.face_b {
            let c = loop_position_centroid(&brep.face_outer_loop_positions(neighbor));
            let along = (c - seed_pa).dot(&na_dir);
            if along > gap * 0.05 && along < gap * 0.95 {
                let d = c - seed_pa;
                let in_plane = d - na_dir * d.dot(&na_dir);
                let perp = in_plane - tangent * in_plane.dot(&tangent);
                if let Some(dir) = perp.normalized() {
                    return Some(dir);
                }
            }
        }
        he = brep.half_edges[he].next;
        if he == start {
            break;
        }
    }
    None
}

fn loop_position_centroid(points: &[Point3]) -> Point3 {
    let n = points.len().max(1) as f64;
    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
    for p in points {
        x += p.x;
        y += p.y;
        z += p.z;
    }
    Point3::new(x / n, y / n, z / n)
}

/// Principal in-plane axis of a boundary polygon: the eigenvector of
/// the planar covariance with the largest extent.
fn dominant_extent_direction(points: &[Point3], normal: Vec3) -> Option<Vec3> {
    let c = loop_position_centroid(points);
    let e1 = points
        .iter()
        .find_map(|p| ((*p - c) - normal * (*p - c).dot(&normal)).normalized())?;
    let e2 = normal.cross(&e1);
    let (mut saa, mut sab, mut sbb) = (0.0_f64, 0.0_f64, 0.0_f64);
    for p in points {
        let d = *p - c;
        let pa = d.dot(&e1);
        let pb = d.dot(&e2);
        saa += pa * pa;
        sab += pa * pb;
        sbb += pb * pb;
    }
    let theta = 0.5 * (2.0 * sab).patan2(saa - sbb);
    (e1 * theta.pcos() + e2 * theta.psin()).normalized()
}
