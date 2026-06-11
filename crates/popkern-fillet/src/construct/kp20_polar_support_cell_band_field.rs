//! Transverse ribbon-field authoring for constructive bands.
//!
//! Builds, per sub-chain, the degree-8 field `R(s)` that spans —
//! together with the chain tangent — every member row's shared
//! tangent plane. Real node directions are interpolated exactly;
//! chained transverse legs are authoritative and constrain the
//! field to contain them.

use std::collections::HashMap;

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::topology::{FaceHandle, VertexHandle};

/// Sample count for the per-sub-chain ribbon-field fit.
const FIELD_FIT_SAMPLES: usize = 17;

use super::cell_reauthor::{AuthoredSubChain, second_row_controls, transverse_row_side};
use super::chain_fit::{ChainRowKey, boundary_controls};
use super::chain_ribbon_field::{eval_ribbon_field, fit_ribbon_field, row_derivative_coefficients};

/// Author the sub-chain's transverse ribbon field: direction
/// targets `W(s) = N(s) × T̂(s)` (the in-plane transverse of the
/// interpolated node-normal field), interpolated exactly at the
/// nodes so the corner solves close.
pub(super) fn author_ribbon_field(
    sub: &AuthoredSubChain,
    node_normals: &HashMap<VertexHandle, Vec3>,
    patches: &HashMap<FaceHandle, RectangularGregoryBiNonic>,
    chained_rows: &std::collections::HashSet<ChainRowKey>,
) -> Option<[Vec3; 9]> {
    let mut normal_field: Vec<(f64, Vec3)> = Vec::with_capacity(sub.nodes.len());
    for &(vertex, param) in &sub.nodes {
        if let Some(&normal) = node_normals.get(&vertex) {
            normal_field.push((param, normal));
        }
    }
    normal_field.sort_by(|a, b| a.0.total_cmp(&b.0));
    normal_field.dedup_by(|a, b| (a.0 - b.0).abs() <= 1.0e-9);
    if normal_field.is_empty() {
        return None;
    }
    for index in 1..normal_field.len() {
        if normal_field[index - 1].1.dot(&normal_field[index].1) < 0.0 {
            normal_field[index].1 = normal_field[index].1 * -1.0;
        }
    }
    let derivative = row_derivative_coefficients(&sub.chain_cps);
    let transverse_at = |t: f64| -> Option<Vec3> {
        let tangent = eval_ribbon_field(&derivative, t).normalized()?;
        let normal = interpolate_normal(&normal_field, t);
        (normal.cross(&tangent)).normalized()
    };
    let mut samples = Vec::with_capacity(FIELD_FIT_SAMPLES);
    let mut previous: Option<Vec3> = None;
    for index in 0..FIELD_FIT_SAMPLES {
        let t = index as f64 / (FIELD_FIT_SAMPLES - 1) as f64;
        let Some(mut direction) = transverse_at(t) else {
            continue;
        };
        if let Some(previous) = previous
            && previous.dot(&direction) < 0.0
        {
            direction = direction * -1.0;
        }
        previous = Some(direction);
        samples.push((t, direction));
    }
    let mut node_constraints = Vec::with_capacity(normal_field.len());
    for &(param, _) in &normal_field {
        let Some(mut direction) = transverse_at(param) else {
            continue;
        };
        // Align with the fitted sample stream's orientation.
        if let Some((_, nearest)) = samples
            .iter()
            .min_by(|a, b| (a.0 - param).abs().total_cmp(&(b.0 - param).abs()))
            && nearest.dot(&direction) < 0.0
        {
            direction = direction * -1.0;
        }
        node_constraints.push((param, direction));
    }
    // Chained transverse legs are authoritative: where a member
    // row's corner leg belongs to a chain-authored transverse
    // row, the field must CONTAIN that leg (its tangent-
    // perpendicular component) so the band derivation closes
    // without moving chain-owned controls. Soft-heavy rows: no
    // budget limit (the hard-budget version silently dropped
    // legs on node-rich sub-chains, leaving frozen endpoints out
    // of span through the corner zones).
    let mut soft_directions: Vec<(f64, Vec3)> = Vec::new();
    for &(key, (start, end)) in &sub.rows {
        for (corner_k, param) in [(0usize, start), (9, end)] {
            if !chained_rows.contains(&ChainRowKey {
                face: key.face,
                side: transverse_row_side(key.side, corner_k),
            }) {
                continue;
            }
            if node_constraints
                .iter()
                .any(|&(existing, _)| (existing - param).abs() <= 1.0e-6)
            {
                continue;
            }
            let Some(patch) = patches.get(&key.face) else {
                continue;
            };
            let row0 = boundary_controls(patch, key.side);
            let row1 = second_row_controls(patch, key.side);
            let leg = row1[corner_k] - row0[corner_k];
            let Some(unit_tangent) = eval_ribbon_field(&derivative, param).normalized() else {
                continue;
            };
            let Some(direction) = (leg - unit_tangent * leg.dot(&unit_tangent)).normalized() else {
                continue;
            };
            soft_directions.push((param, direction));
        }
    }
    fit_ribbon_field(&samples, &node_constraints, &soft_directions).ok()
}

/// Piecewise-linear interpolation of the sign-aligned node
/// normals, renormalized.
fn interpolate_normal(field: &[(f64, Vec3)], param: f64) -> Vec3 {
    let first = field[0];
    let last = field[field.len() - 1];
    if param <= first.0 {
        return first.1;
    }
    if param >= last.0 {
        return last.1;
    }
    for window in field.windows(2) {
        let (p0, n0) = window[0];
        let (p1, n1) = window[1];
        if param <= p1 + 1.0e-9 {
            let span = p1 - p0;
            if span <= 1.0e-9 {
                return n0;
            }
            let fraction = (param - p0) / span;
            let blended = n0 * (1.0 - fraction) + n1 * fraction;
            return blended.normalized().unwrap_or(n0);
        }
    }
    last.1
}
