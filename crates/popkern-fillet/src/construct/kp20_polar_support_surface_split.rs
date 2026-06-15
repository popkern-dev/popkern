//! Support/pole ownership helpers for KP20 support-surface candidates.

use popkern_core::geometry::gregory::{
    Kp20Bi6Patch, Kp20MacroPatch, Kp20PieceId, bernstein_bi6_all, bernstein_bi6_deriv_all,
    kp20_apply_g2_strip, kp20_b_coeffs_near_o, kp20_propagate_c2_upper_from_lower,
};
use popkern_core::geometry::{Point3, Vec3};

pub(super) fn support_owner_piece(piece: Kp20PieceId) -> Kp20PieceId {
    match piece {
        Kp20PieceId::NearOriginLeft | Kp20PieceId::NearMidpointLeft => {
            Kp20PieceId::NearMidpointLeft
        }
        Kp20PieceId::NearOriginRight | Kp20PieceId::NearMidpointRight => {
            Kp20PieceId::NearMidpointRight
        }
    }
}

pub(super) fn support_side_index(piece: Kp20PieceId) -> usize {
    match piece {
        Kp20PieceId::NearMidpointLeft => 0,
        Kp20PieceId::NearMidpointRight => 1,
        _ => unreachable!("support-pole split trials target support-owned pieces"),
    }
}

pub(super) fn propagate_support_pole_interfaces(
    macro_patches: &mut [Kp20MacroPatch],
    changed_patches: &mut [[bool; 4]],
    support_changed: Vec<[bool; 2]>,
) {
    let n = macro_patches.len();
    let b0 = kp20_b_coeffs_near_o((core::f64::consts::TAU / n as f64).cos())[0];
    for (sector, changed) in support_changed.into_iter().enumerate() {
        if changed[0] {
            let lower = *macro_patches[sector].piece(Kp20PieceId::NearMidpointLeft);
            kp20_propagate_c2_upper_from_lower(
                &lower,
                macro_patches[sector].piece_mut(Kp20PieceId::NearOriginLeft),
            );
            changed_patches[sector][Kp20PieceId::NearOriginLeft.index()] = true;
        }
        if changed[1] {
            let lower = *macro_patches[sector].piece(Kp20PieceId::NearMidpointRight);
            kp20_propagate_c2_upper_from_lower(
                &lower,
                macro_patches[sector].piece_mut(Kp20PieceId::NearOriginRight),
            );
            changed_patches[sector][Kp20PieceId::NearOriginRight.index()] = true;
        }
        if changed[0] || changed[1] {
            let [left, right, _, _] = &mut macro_patches[sector].pieces;
            kp20_apply_g2_strip(left, right, b0);
        }
    }
}

pub(super) fn evaluate_bi6_patch(piece: &Kp20Bi6Patch, u: f64, v: f64) -> Point3 {
    let bu = bernstein_bi6_all(u.clamp(0.0, 1.0));
    let bv = bernstein_bi6_all(v.clamp(0.0, 1.0));
    let mut out = Point3::new(0.0, 0.0, 0.0);
    for (i, &bu_i) in bu.iter().enumerate() {
        for (j, &bv_j) in bv.iter().enumerate() {
            let w = bu_i * bv_j;
            out.x += w * piece.cps[i][j].x;
            out.y += w * piece.cps[i][j].y;
            out.z += w * piece.cps[i][j].z;
        }
    }
    out
}

pub(super) fn evaluate_bi6_normal(piece: &Kp20Bi6Patch, u: f64, v: f64) -> Option<Vec3> {
    let u = u.clamp(0.0, 1.0);
    let v = v.clamp(0.0, 1.0);
    let bu = bernstein_bi6_all(u);
    let bv = bernstein_bi6_all(v);
    let du = bernstein_bi6_deriv_all(u);
    let dv = bernstein_bi6_deriv_all(v);
    let mut tangent_u = Vec3::new(0.0, 0.0, 0.0);
    let mut tangent_v = Vec3::new(0.0, 0.0, 0.0);
    for (i, (&bu_i, &du_i)) in bu.iter().zip(du.iter()).enumerate() {
        for (j, (&bv_j, &dv_j)) in bv.iter().zip(dv.iter()).enumerate() {
            let p = piece.cps[i][j];
            tangent_u = tangent_u + (p - Point3::new(0.0, 0.0, 0.0)) * (du_i * bv_j);
            tangent_v = tangent_v + (p - Point3::new(0.0, 0.0, 0.0)) * (bu_i * dv_j);
        }
    }
    tangent_u.cross(&tangent_v).normalized()
}
