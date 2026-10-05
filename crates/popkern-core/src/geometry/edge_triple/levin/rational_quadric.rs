use num_rational::BigRational;
use num_traits::{One, Zero};
use popkern_math::PortableMath;

use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::{Point3, Vec3};

use super::rational::{Mat4SymRat, Rat};

/// Convert an analytic surface to its exact rational quadric form.
/// Returns `None` for unsupported surfaces or non-finite parameters.
pub(crate) fn surface_to_quadric_rat(surf: &SurfaceGeometry) -> Option<Mat4SymRat> {
    match surf {
        SurfaceGeometry::Plane { origin, normal, .. } => Some(plane_quadric_rat(*origin, *normal)?),
        SurfaceGeometry::Sphere { center, radius } => Some(sphere_quadric_rat(*center, *radius)?),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => Some(cylinder_quadric_rat(*origin, *axis, *radius)?),
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => Some(cone_quadric_rat(*apex, *axis, *half_angle)?),
        _ => None,
    }
}

pub(super) fn mat4_zero() -> [[Rat; 4]; 4] {
    let z = Rat::zero();
    [
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z.clone()],
        [z.clone(), z.clone(), z.clone(), z],
    ]
}

fn rat_from_f64(x: f64) -> Option<Rat> {
    BigRational::from_float(x)
}

fn rat_one() -> Rat {
    Rat::one()
}

fn plane_quadric_rat(origin: Point3, normal: Vec3) -> Option<Mat4SymRat> {
    let nx = rat_from_f64(normal.x)?;
    let ny = rat_from_f64(normal.y)?;
    let nz = rat_from_f64(normal.z)?;
    let ox = rat_from_f64(origin.x)?;
    let oy = rat_from_f64(origin.y)?;
    let oz = rat_from_f64(origin.z)?;
    let np = &nx * &ox + &ny * &oy + &nz * &oz;
    let v = [nx, ny, nz, -np];
    let mut a = mat4_zero();
    for (i, row) in a.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            *slot = &v[i] * &v[j];
        }
    }
    Some(Mat4SymRat(a))
}

fn sphere_quadric_rat(center: Point3, radius: f64) -> Option<Mat4SymRat> {
    let cx = rat_from_f64(center.x)?;
    let cy = rat_from_f64(center.y)?;
    let cz = rat_from_f64(center.z)?;
    let r = rat_from_f64(radius)?;
    let cdot = &cx * &cx + &cy * &cy + &cz * &cz;
    let r_sq = &r * &r;
    let mut a = mat4_zero();
    a[0][0] = rat_one();
    a[1][1] = rat_one();
    a[2][2] = rat_one();
    a[0][3] = -&cx;
    a[3][0] = -cx;
    a[1][3] = -&cy;
    a[3][1] = -cy;
    a[2][3] = -&cz;
    a[3][2] = -cz;
    a[3][3] = cdot - r_sq;
    Some(Mat4SymRat(a))
}

fn cylinder_quadric_rat(origin: Point3, axis: Vec3, radius: f64) -> Option<Mat4SymRat> {
    let ax = rat_from_f64(axis.x)?;
    let ay = rat_from_f64(axis.y)?;
    let az = rat_from_f64(axis.z)?;
    let ox = rat_from_f64(origin.x)?;
    let oy = rat_from_f64(origin.y)?;
    let oz = rat_from_f64(origin.z)?;
    let r = rat_from_f64(radius)?;
    let m00 = rat_one() - &ax * &ax;
    let m01 = -(&ax * &ay);
    let m02 = -(&ax * &az);
    let m11 = rat_one() - &ay * &ay;
    let m12 = -(&ay * &az);
    let m22 = rat_one() - &az * &az;
    let m = [
        [m00.clone(), m01.clone(), m02.clone()],
        [m01.clone(), m11.clone(), m12.clone()],
        [m02.clone(), m12.clone(), m22.clone()],
    ];
    let o = [ox.clone(), oy.clone(), oz.clone()];
    let mo = mat3_vec3_rat(&m, &o);
    let omo = &o[0] * &mo[0] + &o[1] * &mo[1] + &o[2] * &mo[2];
    let r_sq = &r * &r;
    let mut a = mat4_zero();
    a[0][0] = m00;
    a[0][1] = m01.clone();
    a[1][0] = m01;
    a[0][2] = m02.clone();
    a[2][0] = m02;
    a[1][1] = m11;
    a[1][2] = m12.clone();
    a[2][1] = m12;
    a[2][2] = m22;
    a[0][3] = -mo[0].clone();
    a[3][0] = -mo[0].clone();
    a[1][3] = -mo[1].clone();
    a[3][1] = -mo[1].clone();
    a[2][3] = -mo[2].clone();
    a[3][2] = -mo[2].clone();
    a[3][3] = omo - r_sq;
    Some(Mat4SymRat(a))
}

fn cone_quadric_rat(apex: Point3, axis: Vec3, half_angle: f64) -> Option<Mat4SymRat> {
    let cos2 = rat_from_f64(half_angle.pcos().powi(2))?;
    let ax = rat_from_f64(axis.x)?;
    let ay = rat_from_f64(axis.y)?;
    let az = rat_from_f64(axis.z)?;
    let vx = rat_from_f64(apex.x)?;
    let vy = rat_from_f64(apex.y)?;
    let vz = rat_from_f64(apex.z)?;
    let q00 = &cos2 - &ax * &ax;
    let q01 = -(&ax * &ay);
    let q02 = -(&ax * &az);
    let q11 = &cos2 - &ay * &ay;
    let q12 = -(&ay * &az);
    let q22 = &cos2 - &az * &az;
    let q = [
        [q00.clone(), q01.clone(), q02.clone()],
        [q01.clone(), q11.clone(), q12.clone()],
        [q02.clone(), q12.clone(), q22.clone()],
    ];
    let v = [vx, vy, vz];
    let qv = mat3_vec3_rat(&q, &v);
    let vqv = &v[0] * &qv[0] + &v[1] * &qv[1] + &v[2] * &qv[2];
    let mut a = mat4_zero();
    a[0][0] = q00;
    a[0][1] = q01.clone();
    a[1][0] = q01;
    a[0][2] = q02.clone();
    a[2][0] = q02;
    a[1][1] = q11;
    a[1][2] = q12.clone();
    a[2][1] = q12;
    a[2][2] = q22;
    a[0][3] = -qv[0].clone();
    a[3][0] = -qv[0].clone();
    a[1][3] = -qv[1].clone();
    a[3][1] = -qv[1].clone();
    a[2][3] = -qv[2].clone();
    a[3][2] = -qv[2].clone();
    a[3][3] = vqv;
    Some(Mat4SymRat(a))
}

fn mat3_vec3_rat(m: &[[Rat; 3]; 3], v: &[Rat; 3]) -> [Rat; 3] {
    [
        &m[0][0] * &v[0] + &m[0][1] * &v[1] + &m[0][2] * &v[2],
        &m[1][0] * &v[0] + &m[1][1] * &v[1] + &m[1][2] * &v[2],
        &m[2][0] * &v[0] + &m[2][1] * &v[1] + &m[2][2] * &v[2],
    ]
}
