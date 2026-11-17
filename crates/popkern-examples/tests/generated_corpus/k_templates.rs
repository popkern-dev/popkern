//! Curve-tangency templates (configuration K, design §4): B tangent to A
//! along a ruling or a circle, placed from the operands' own geometry.

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Transform3, Vec3};

use super::operand::{Kind, Operand};
use super::rng::Rng;

/// A placed curve tangency: B's rigid placement and uniform scale,
/// whether the normals are aligned (B inside A along the curve), the
/// contact normal (A's outward), and samples of the contact curve, all
/// in A's frame.
pub struct Tangency {
    pub place: Transform3,
    pub b_scale: f64,
    pub aligned: bool,
    pub normal: Vec3,
    pub curve: Vec<Point3>,
}

fn surface(op: &Operand, pick: impl Fn(&SurfaceGeometry) -> bool) -> Option<&SurfaceGeometry> {
    op.brep
        .faces
        .values()
        .filter_map(|f| f.surface.as_ref().map(|s| &s.geometry))
        .find(|g| pick(g))
}

fn cylinder(op: &Operand) -> Option<(Point3, Vec3, f64)> {
    match surface(op, |g| matches!(g, SurfaceGeometry::Cylinder { .. }))? {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some((*origin, *axis, *radius)),
        _ => None,
    }
}

fn sphere(op: &Operand) -> Option<(Point3, f64)> {
    match surface(op, |g| matches!(g, SurfaceGeometry::Sphere { .. }))? {
        SurfaceGeometry::Sphere { center, radius, .. } => Some((*center, *radius)),
        _ => None,
    }
}

fn cone(op: &Operand) -> Option<(Point3, Vec3, f64)> {
    match surface(op, |g| matches!(g, SurfaceGeometry::Cone { .. }))? {
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => Some((*apex, *axis, *half_angle)),
        _ => None,
    }
}

fn torus(op: &Operand) -> Option<(f64, f64)> {
    match surface(op, |g| matches!(g, SurfaceGeometry::Torus { .. }))? {
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } => Some((*major_radius, *minor_radius)),
        _ => None,
    }
}

fn circle(c: Point3, r: f64, z: f64) -> Vec<Point3> {
    (0..16)
        .map(|k| {
            let t = k as f64 * std::f64::consts::TAU / 16.0;
            Point3::new(c.x + r * t.cos(), c.y + r * t.sin(), z)
        })
        .collect()
}

fn segment(a: Point3, b: Point3) -> Vec<Point3> {
    (0..=8)
        .map(|k| {
            let t = k as f64 / 8.0;
            Point3::new(
                a.x + (b.x - a.x) * t,
                a.y + (b.y - a.y) * t,
                a.z + (b.z - a.z) * t,
            )
        })
        .collect()
}

/// The template for `(a, b)` with B moved `delta` apart along the
/// contact normal (negative: overlapping). `None` when the pair has none.
pub fn place(a: &Operand, b: &Operand, delta: f64, rng: &mut Rng) -> Option<Tangency> {
    use Kind as K;
    let (lo, hi) = super::distance::mesh_box(&a.mesh);
    struct Bb {
        min: Point3,
        max: Point3,
    }
    let bb = Bb { min: lo, max: hi };
    let up = Vec3::new(0.0, 0.0, 1.0);
    let origin = Point3::new(0.0, 0.0, 0.0);
    let x_axis = Vec3::new(1.0, 0.0, 0.0);
    match (a.kind, b.kind) {
        // A cylinder, a cone or a frustum lying on the box's top along a
        // ruling, placed exactly (the mesh's lowest vertex would float a
        // sagitta above the plane when the ruling falls between stations).
        (K::P, K::Y) | (K::P, K::C) | (K::P, K::F) => {
            let top = bb.max.z;
            let len = 0.5 * (bb.max.x - bb.min.x).min(bb.max.y - bb.min.y);
            let b_scale = (len / b.size).min(1.0);
            let (lay, e0, e1) = lying(b, b_scale)?;
            let shift_y = rng.range(-0.1, 0.1) * len;
            let lift = Transform3::translation(0.0, shift_y, top - e0.z + delta);
            let place = lift.compose(&lay);
            let curve = segment(
                Point3::new(e0.x, e0.y + shift_y, top),
                Point3::new(e1.x, e1.y + shift_y, top),
            );
            Some(Tangency {
                place,
                b_scale,
                aligned: false,
                normal: up,
                curve,
            })
        }
        // Parallel cylinders touching along a line.
        (K::Y, K::Y) => {
            let (_, _, ra) = cylinder(a)?;
            let (_, _, rb) = cylinder(b)?;
            let gap = ra + rb + delta;
            let place = Transform3::translation(gap, 0.0, rng.range(-0.2, 0.2) * a.size);
            let curve = segment(
                Point3::new(ra, 0.0, bb.min.z),
                Point3::new(ra, 0.0, bb.max.z),
            );
            Some(Tangency {
                place,
                b_scale: 1.0,
                aligned: false,
                normal: x_axis,
                curve,
            })
        }
        // A torus lying flat on the box's top along its lowest circle.
        (K::P, K::T) => {
            let (big, small) = torus(b)?;
            let b_scale =
                (0.45 * (bb.max.x - bb.min.x).min(bb.max.y - bb.min.y) / (big + small)).min(1.0);
            let top = bb.max.z;
            let place = Transform3::translation(0.0, 0.0, top + small * b_scale + delta);
            let curve = circle(origin, big * b_scale, top);
            Some(Tangency {
                place,
                b_scale,
                aligned: false,
                normal: up,
                curve,
            })
        }
        // A ball seated in the torus's hole, along the tube's inner equator.
        (K::T, K::S) => {
            let (big, small) = torus(a)?;
            let (_, rb) = sphere(b)?;
            let rho = big - small - delta;
            let b_scale = rho / rb;
            let curve = circle(origin, big - small, 0.0);
            Some(Tangency {
                place: Transform3::identity(),
                b_scale,
                aligned: false,
                normal: x_axis * -1.0,
                curve,
            })
        }
        // A ball inside the cylinder at its radius, along the equator.
        (K::Y, K::S) => {
            let (o, _, r) = cylinder(a)?;
            let (_, rb) = sphere(b)?;
            let b_scale = (r - delta) / rb;
            let zc = 0.5 * (bb.min.z + bb.max.z)
                + rng.range(-0.1, 0.1) * (bb.max.z - bb.min.z - 2.0 * r).max(0.0);
            let place = Transform3::translation(o.x, o.y, zc);
            Some(Tangency {
                place,
                b_scale,
                aligned: true,
                normal: x_axis,
                curve: circle(Point3::new(o.x, o.y, 0.0), r, zc),
            })
        }
        // A ball inside the cone, along the circle where it is tangent.
        (K::C, K::S) => {
            let (apex, axis, ha) = cone(a)?;
            let (_, rb) = sphere(b)?;
            // Base radius at the bottom; a ball of a third of it.
            let rho = 0.3 * (bb.max.x - bb.min.x) * 0.5;
            let d = rho / ha.sin();
            // `axis` points from the apex into the cone.
            let c = apex + axis * d;
            let rho = rho - delta;
            let b_scale = rho / rb;
            let place = Transform3::translation(c.x, c.y, c.z);
            let ring_z = (apex + axis * (d * ha.cos() * ha.cos())).z;
            let ring_r = d * ha.cos() * ha.sin();
            let normal = Vec3::new(1.0, 0.0, 0.0);
            Some(Tangency {
                place,
                b_scale,
                aligned: true,
                normal,
                curve: circle(Point3::new(c.x, c.y, 0.0), ring_r, ring_z),
            })
        }
        _ => None,
    }
}

/// Lay a surface of revolution about z (cylinder, cone, frustum) on its
/// side: a rotation about y bringing the ruling of the plane y = 0 on one
/// side horizontal and lowest. Returns the rotation (B unscaled in, its
/// scale applied by the caller's placement) and the ruling's ends in the
/// laid, scaled frame. The rims are the planar caps (exact in the mesh)
/// and the side surface's own radius at their heights.
fn lying(b: &Operand, b_scale: f64) -> Option<(Transform3, Point3, Point3)> {
    let zs: Vec<f64> = b.mesh.vertices.iter().map(|v| v.position.z).collect();
    let (z0, z1) = zs
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), z| {
            (lo.min(*z), hi.max(*z))
        });
    let radius_at = |z: f64| -> Option<f64> {
        match b.kind {
            Kind::Y => cylinder(b).map(|(_, _, r)| r),
            _ => cone(b).map(|(apex, _, ha)| (apex.z - z).abs() * ha.tan()),
        }
    };
    let (r0, r1) = (radius_at(z0)?, radius_at(z1)?);
    let origin = Point3::new(0.0, 0.0, 0.0);
    for side in [1.0, -1.0] {
        let (p0, p1) = (
            Point3::new(side * r0, 0.0, z0),
            Point3::new(side * r1, 0.0, z1),
        );
        let d = p1 - p0;
        // Rotating about +y by θ maps (x, z) to (x cos θ + z sin θ,
        // −x sin θ + z cos θ): the ruling is horizontal where
        // −dx sin θ + dz cos θ = 0.
        for theta in [d.z.atan2(d.x), d.z.atan2(d.x) + std::f64::consts::PI] {
            let turn =
                Transform3::rotation_around_axis(origin, Vec3::new(0.0, 1.0, 0.0), theta).ok()?;
            let laid = turn.compose(&Transform3::scale_uniform(b_scale));
            let (q0, q1) = (laid.apply_point(&p0), laid.apply_point(&p1));
            let lowest_mesh = b
                .mesh
                .vertices
                .iter()
                .map(|v| laid.apply_point(&v.position).z)
                .fold(f64::INFINITY, f64::min);
            if (q0.z - q1.z).abs() <= 1e-9 * b.size && lowest_mesh >= q0.z - 4.0 * b.chord * b_scale
            {
                return Some((turn, q0, q1));
            }
        }
    }
    None
}
