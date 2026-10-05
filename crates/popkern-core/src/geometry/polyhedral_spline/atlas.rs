use super::types::{PatchEdge, PolyhedralPatch, PolyhedralSpline, PolyhedralSplineError};
use crate::geometry::{Point3, Vec3};
use popkern_math::PortableMath;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeamOrientation {
    Same,
    Reversed,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolyhedralSeam {
    pub patch_a: usize,
    pub edge_a: PatchEdge,
    pub patch_b: usize,
    pub edge_b: PatchEdge,
    pub orientation: SeamOrientation,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PolyhedralPatchAtlas {
    pub patch_count: usize,
    pub seams: Vec<PolyhedralSeam>,
    pub boundary_edges: Vec<(usize, PatchEdge)>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct SeamContinuityReport {
    pub seam: PolyhedralSeam,
    pub samples: usize,
    pub max_position_gap: f64,
    pub max_boundary_tangent_angle: f64,
    pub max_normal_angle: f64,
    pub c0_tolerance: f64,
    pub angle_tolerance: f64,
    pub accepted: bool,
}

pub fn build_patch_atlas(
    spline: &PolyhedralSpline,
    samples: usize,
    tolerance: f64,
) -> PolyhedralPatchAtlas {
    let samples = samples.max(2);
    let mut seams = Vec::new();
    let mut used = vec![[false; 4]; spline.patches.len()];
    for a in 0..spline.patches.len() {
        for b in (a + 1)..spline.patches.len() {
            for edge_a in all_edges() {
                for edge_b in all_edges() {
                    if let Some(orientation) = matching_orientation(
                        &spline.patches[a],
                        edge_a,
                        &spline.patches[b],
                        edge_b,
                        samples,
                        tolerance,
                    ) {
                        seams.push(PolyhedralSeam {
                            patch_a: a,
                            edge_a,
                            patch_b: b,
                            edge_b,
                            orientation,
                        });
                        used[a][edge_index(edge_a)] = true;
                        used[b][edge_index(edge_b)] = true;
                    }
                }
            }
        }
    }

    let mut boundary_edges = Vec::new();
    for (patch, edge_used) in used.iter().enumerate() {
        for edge in all_edges() {
            if !edge_used[edge_index(edge)] {
                boundary_edges.push((patch, edge));
            }
        }
    }
    PolyhedralPatchAtlas {
        patch_count: spline.patches.len(),
        seams,
        boundary_edges,
    }
}

pub fn validate_spline_atlas(
    spline: &PolyhedralSpline,
    samples: usize,
    c0_tolerance: f64,
    angle_tolerance: f64,
) -> Result<Vec<SeamContinuityReport>, PolyhedralSplineError> {
    let atlas = build_patch_atlas(spline, samples, c0_tolerance);
    let mut reports = Vec::with_capacity(atlas.seams.len());
    for seam in atlas.seams {
        let report = validate_seam(
            &spline.patches[seam.patch_a],
            &spline.patches[seam.patch_b],
            seam,
            samples,
            c0_tolerance,
            angle_tolerance,
        );
        if !report.accepted {
            return Err(PolyhedralSplineError::SeamValidation {
                family: spline.family,
                reason: format!(
                    "patch {} {:?} to patch {} {:?}: gap={}, tangent_angle={}, normal_angle={}",
                    report.seam.patch_a,
                    report.seam.edge_a,
                    report.seam.patch_b,
                    report.seam.edge_b,
                    report.max_position_gap,
                    report.max_boundary_tangent_angle,
                    report.max_normal_angle
                ),
            });
        }
        reports.push(report);
    }
    Ok(reports)
}

fn validate_seam(
    a: &PolyhedralPatch,
    b: &PolyhedralPatch,
    seam: PolyhedralSeam,
    samples: usize,
    c0_tolerance: f64,
    angle_tolerance: f64,
) -> SeamContinuityReport {
    let samples = samples.max(2);
    let mut max_position_gap = 0.0_f64;
    let mut max_boundary_tangent_angle = 0.0_f64;
    let mut max_normal_angle = 0.0_f64;
    let singular_edge =
        edge_span(a, seam.edge_a) <= c0_tolerance && edge_span(b, seam.edge_b) <= c0_tolerance;
    for i in 0..samples {
        let t = i as f64 / (samples - 1) as f64;
        let s = orient_t(t, seam.orientation);
        let pa = eval_edge(a, seam.edge_a, t);
        let pb = eval_edge(b, seam.edge_b, s);
        max_position_gap = max_position_gap.max(pa.distance_to(&pb));
        if !singular_edge {
            max_boundary_tangent_angle = max_boundary_tangent_angle.max(angle_between(
                edge_tangent(a, seam.edge_a, t),
                edge_tangent(b, seam.edge_b, s) * orient_sign(seam.orientation),
            ));
            max_normal_angle = max_normal_angle.max(angle_between(
                tangent_plane_normal(a, seam.edge_a, t),
                tangent_plane_normal(b, seam.edge_b, s),
            ));
        }
    }
    SeamContinuityReport {
        seam,
        samples,
        max_position_gap,
        max_boundary_tangent_angle,
        max_normal_angle,
        c0_tolerance,
        angle_tolerance,
        accepted: max_position_gap <= c0_tolerance
            && max_boundary_tangent_angle <= angle_tolerance
            && max_normal_angle <= angle_tolerance,
    }
}

fn edge_span(patch: &PolyhedralPatch, edge: PatchEdge) -> f64 {
    eval_edge(patch, edge, 0.0).distance_to(&eval_edge(patch, edge, 1.0))
}

fn matching_orientation(
    a: &PolyhedralPatch,
    edge_a: PatchEdge,
    b: &PolyhedralPatch,
    edge_b: PatchEdge,
    samples: usize,
    tolerance: f64,
) -> Option<SeamOrientation> {
    for orientation in [SeamOrientation::Same, SeamOrientation::Reversed] {
        let mut max_gap = 0.0_f64;
        for i in 0..samples {
            let t = i as f64 / (samples - 1) as f64;
            let pa = eval_edge(a, edge_a, t);
            let pb = eval_edge(b, edge_b, orient_t(t, orientation));
            max_gap = max_gap.max(pa.distance_to(&pb));
        }
        if max_gap <= tolerance {
            return Some(orientation);
        }
    }
    None
}

fn eval_edge(patch: &PolyhedralPatch, edge: PatchEdge, t: f64) -> Point3 {
    match edge {
        PatchEdge::UMin => patch.evaluate(0.0, t),
        PatchEdge::UMax => patch.evaluate(1.0, t),
        PatchEdge::VMin => patch.evaluate(t, 0.0),
        PatchEdge::VMax => patch.evaluate(t, 1.0),
    }
}

fn edge_tangent(patch: &PolyhedralPatch, edge: PatchEdge, t: f64) -> Vec3 {
    match edge {
        PatchEdge::UMin | PatchEdge::UMax => patch.partial_v(edge_u(edge), t),
        PatchEdge::VMin | PatchEdge::VMax => patch.partial_u(t, edge_v(edge)),
    }
}

fn tangent_plane_normal(patch: &PolyhedralPatch, edge: PatchEdge, t: f64) -> Vec3 {
    let (u, v) = match edge {
        PatchEdge::UMin => (0.0, t),
        PatchEdge::UMax => (1.0, t),
        PatchEdge::VMin => (t, 0.0),
        PatchEdge::VMax => (t, 1.0),
    };
    patch.partial_u(u, v).cross(&patch.partial_v(u, v))
}

fn angle_between(a: Vec3, b: Vec3) -> f64 {
    let an = a.length();
    let bn = b.length();
    if an <= 1e-14 || bn <= 1e-14 {
        return 0.0;
    }
    let cos = (a.dot(&b) / (an * bn)).clamp(-1.0, 1.0).abs();
    cos.pacos()
}

fn edge_u(edge: PatchEdge) -> f64 {
    match edge {
        PatchEdge::UMin => 0.0,
        PatchEdge::UMax => 1.0,
        PatchEdge::VMin | PatchEdge::VMax => unreachable!(),
    }
}

fn edge_v(edge: PatchEdge) -> f64 {
    match edge {
        PatchEdge::VMin => 0.0,
        PatchEdge::VMax => 1.0,
        PatchEdge::UMin | PatchEdge::UMax => unreachable!(),
    }
}

fn orient_t(t: f64, orientation: SeamOrientation) -> f64 {
    match orientation {
        SeamOrientation::Same => t,
        SeamOrientation::Reversed => 1.0 - t,
    }
}

fn orient_sign(orientation: SeamOrientation) -> f64 {
    match orientation {
        SeamOrientation::Same => 1.0,
        SeamOrientation::Reversed => -1.0,
    }
}

fn all_edges() -> [PatchEdge; 4] {
    [
        PatchEdge::UMin,
        PatchEdge::UMax,
        PatchEdge::VMin,
        PatchEdge::VMax,
    ]
}

fn edge_index(edge: PatchEdge) -> usize {
    match edge {
        PatchEdge::UMin => 0,
        PatchEdge::UMax => 1,
        PatchEdge::VMin => 2,
        PatchEdge::VMax => 3,
    }
}
