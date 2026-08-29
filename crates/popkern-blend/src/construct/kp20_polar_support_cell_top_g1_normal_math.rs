//! Linearized normal-equation helpers for the top-G1 solve.

use std::collections::HashMap;

use popkern_core::error::CoreError;
use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::{
    RectangularGregoryBiNonic, bernstein_binonic_all, bernstein_binonic_deriv_all,
};
use popkern_core::topology::{BRep, FaceHandle};

use super::worst_g1_normal_target::{ControlKey, ControlVariable};
use super::{rectangular_binonic_on_face, set_rectangular_binonic_on_face};

#[derive(Debug, Clone, Copy)]
pub(super) struct NormalSample {
    pub(super) normal: Vec3,
    raw_normal_magnitude: f64,
    tangent_u: Vec3,
    tangent_v: Vec3,
    u_basis: [f64; 10],
    v_basis: [f64; 10],
    u_deriv_basis: [f64; 10],
    v_deriv_basis: [f64; 10],
}

pub(super) fn normal_sample(
    patch: &RectangularGregoryBiNonic,
    uv: (f64, f64),
) -> Option<NormalSample> {
    let u_basis = bernstein_binonic_all(uv.0);
    let v_basis = bernstein_binonic_all(uv.1);
    let u_deriv_basis = bernstein_binonic_deriv_all(uv.0);
    let v_deriv_basis = bernstein_binonic_deriv_all(uv.1);
    let tangent_u = patch.partial_u(uv.0, uv.1);
    let tangent_v = patch.partial_v(uv.0, uv.1);
    let raw_normal = tangent_u.cross(&tangent_v);
    let raw_normal_magnitude = raw_normal.dot(&raw_normal).sqrt();
    Some(NormalSample {
        normal: raw_normal.normalized()?,
        raw_normal_magnitude,
        tangent_u,
        tangent_v,
        u_basis,
        v_basis,
        u_deriv_basis,
        v_deriv_basis,
    })
}

pub(super) fn normal_sensitivity(sample: NormalSample, variable: ControlVariable) -> Vec3 {
    let direction = axis_direction(variable.axis);
    let delta_u =
        direction * (sample.u_deriv_basis[variable.key.u] * sample.v_basis[variable.key.v]);
    let delta_v =
        direction * (sample.u_basis[variable.key.u] * sample.v_deriv_basis[variable.key.v]);
    let raw_delta = delta_u.cross(&sample.tangent_v) + sample.tangent_u.cross(&delta_v);
    (raw_delta - sample.normal * raw_delta.dot(&sample.normal))
        * (1.0 / sample.raw_normal_magnitude)
}

pub(super) fn limited_control_deltas(
    variables: &[ControlVariable],
    solution: &[f64],
    limit: f64,
) -> (HashMap<ControlKey, Vec3>, usize, f64) {
    let mut deltas: HashMap<ControlKey, Vec3> = HashMap::new();
    for (variable, value) in variables.iter().zip(solution) {
        deltas
            .entry(variable.key)
            .and_modify(|delta| *delta = *delta + axis_direction(variable.axis) * *value)
            .or_insert(axis_direction(variable.axis) * *value);
    }
    let mut limited_count = 0usize;
    let mut max_delta = 0.0_f64;
    for delta in deltas.values_mut() {
        let magnitude = delta.dot(delta).sqrt();
        if magnitude > limit {
            *delta = *delta * (limit / magnitude);
            limited_count += 1;
            max_delta = max_delta.max(limit);
        } else {
            max_delta = max_delta.max(magnitude);
        }
    }
    (deltas, limited_count, max_delta)
}

pub(super) fn apply_control_deltas(
    brep: &mut BRep,
    deltas: &HashMap<ControlKey, Vec3>,
    damping: f64,
) -> Result<(), CoreError> {
    let mut patches: HashMap<FaceHandle, RectangularGregoryBiNonic> = HashMap::new();
    for (&key, delta) in deltas {
        if let std::collections::hash_map::Entry::Vacant(entry) = patches.entry(key.face) {
            entry.insert(rectangular_binonic_on_face(brep, key.face)?);
        }
        let Some(patch) = patches.get_mut(&key.face) else {
            return Err(CoreError::PreconditionFailed(
                "kp20 top-G1 normal solve control patch missing".into(),
            ));
        };
        patch.grid[key.u][key.v] = patch.grid[key.u][key.v] + *delta * damping;
    }
    for (face, patch) in patches {
        set_rectangular_binonic_on_face(brep, face, patch)?;
    }
    Ok(())
}

pub(super) fn add_normal_equation(
    matrix: &mut [Vec<f64>],
    rhs: &mut [f64],
    coefficients: &[f64],
    target: f64,
) {
    for i in 0..coefficients.len() {
        if coefficients[i].abs() <= 1.0e-18 {
            continue;
        }
        rhs[i] += coefficients[i] * target;
        for j in 0..coefficients.len() {
            if coefficients[j].abs() > 1.0e-18 {
                matrix[i][j] += coefficients[i] * coefficients[j];
            }
        }
    }
}

pub(super) fn component(vector: Vec3, axis: usize) -> f64 {
    match axis {
        0 => vector.x,
        1 => vector.y,
        _ => vector.z,
    }
}

fn axis_direction(axis: usize) -> Vec3 {
    match axis {
        0 => Vec3::new(1.0, 0.0, 0.0),
        1 => Vec3::new(0.0, 1.0, 0.0),
        _ => Vec3::new(0.0, 0.0, 1.0),
    }
}
