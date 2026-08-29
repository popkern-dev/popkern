//! Endpoint-local row fitting for KP20 support-cell boundary clusters.

use popkern_core::error::CoreError;
use popkern_core::geometry::Point3;
use popkern_core::topology::BRep;

use super::super::{
    BoundarySide, FullBoundaryTrim, rectangular_binonic_on_face, row_fit,
    set_rectangular_binonic_on_face,
};
use super::{PatchCornerKey, average_point, lerp_scalar, trim_endpoint_pairs};

pub(super) fn fit_boundary_bands(
    brep: &mut BRep,
    cluster: &[PatchCornerKey],
    trims: &[FullBoundaryTrim],
) -> Result<(), CoreError> {
    for &trim in trims {
        for (corner_a, corner_b) in trim_endpoint_pairs(trim) {
            if !cluster.contains(&corner_a) || !cluster.contains(&corner_b) {
                continue;
            }
            let patch_a = rectangular_binonic_on_face(brep, trim.face_a)?;
            let patch_b = rectangular_binonic_on_face(brep, trim.face_b)?;
            let samples = endpoint_band_samples(&patch_a, &patch_b, trim, corner_a);
            let mut candidate_a = patch_a.clone();
            let mut candidate_b = patch_b.clone();
            row_fit::fit_boundary_samples(&mut candidate_a, trim.side_a, &samples.a);
            row_fit::fit_boundary_samples(&mut candidate_b, trim.side_b, &samples.b);
            set_rectangular_binonic_on_face(brep, trim.face_a, candidate_a)?;
            set_rectangular_binonic_on_face(brep, trim.face_b, candidate_b)?;
        }
    }
    Ok(())
}

struct EndpointBandSamples {
    a: [row_fit::BoundaryRowFitSample; 4],
    b: [row_fit::BoundaryRowFitSample; 4],
}

fn endpoint_band_samples(
    patch_a: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    patch_b: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    trim: FullBoundaryTrim,
    corner_a: PatchCornerKey,
) -> EndpointBandSamples {
    let endpoint_fraction = endpoint_parameter_fraction(corner_a, trim.side_a);
    let direction = if endpoint_fraction > 0.5 { -1.0 } else { 1.0 };
    let fractions: [f64; 4] = [
        endpoint_fraction,
        endpoint_fraction + direction / 9.0,
        endpoint_fraction + direction * 2.0 / 9.0,
        endpoint_fraction + direction / 3.0,
    ];
    let a = std::array::from_fn(|index| {
        let fraction = fractions[index].clamp(0.0, 1.0);
        let param_a = lerp_scalar(trim.start_a, trim.end_a, fraction);
        let param_b = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
        row_fit::BoundaryRowFitSample {
            parameter: param_a,
            target: endpoint_band_target(patch_a, patch_b, trim, param_a, param_b),
        }
    });
    let b = std::array::from_fn(|index| {
        let fraction = fractions[index].clamp(0.0, 1.0);
        let param_a = lerp_scalar(trim.start_a, trim.end_a, fraction);
        let param_b = lerp_scalar(trim.start_b, trim.end_b, 1.0 - fraction);
        row_fit::BoundaryRowFitSample {
            parameter: param_b,
            target: endpoint_band_target(patch_a, patch_b, trim, param_a, param_b),
        }
    });
    EndpointBandSamples { a, b }
}

fn endpoint_parameter_fraction(corner: PatchCornerKey, side: BoundarySide) -> f64 {
    let at_max = match side {
        BoundarySide::UMin | BoundarySide::UMax => corner.v_max,
        BoundarySide::VMin | BoundarySide::VMax => corner.u_max,
    };
    if at_max { 1.0 } else { 0.0 }
}

fn endpoint_band_target(
    patch_a: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    patch_b: &popkern_core::geometry::gregory::RectangularGregoryBiNonic,
    trim: FullBoundaryTrim,
    param_a: f64,
    param_b: f64,
) -> Point3 {
    average_point(
        row_fit::eval_boundary(patch_a, trim.side_a, param_a),
        row_fit::eval_boundary(patch_b, trim.side_b, param_b),
    )
}
