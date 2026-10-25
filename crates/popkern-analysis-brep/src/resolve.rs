//! Resolving a setup's face references against a model.
//!
//! Selectors match face surface geometry (exact for planes, cylinders
//! and spheres), not triangles, except `NearestTo`, which measures the
//! distance from its point to each face's triangles. Faces are visited
//! in `FaceId` order and results are sorted, so resolution is
//! deterministic. A reference that matches a different number of faces
//! than it expects is a typed decline, never a silent pick.

use std::collections::BTreeSet;

use popkern_analysis::report::{ResolvedFace, Warning};
use popkern_analysis::{
    AnalysisError, AnalysisModel, Expect, FaceDef, FaceId, FaceRef, Selector, Surface, Tolerance,
};
use popkern_math::PortableMath;

use crate::geometry::{dot, norm, point_line_distance, point_triangle_distance, sub, unit};

/// Faces each definition resolved to, plus warnings, in setup order.
#[derive(Debug, Clone, PartialEq)]
pub struct Resolution {
    /// One entry per face definition.
    pub faces: Vec<ResolvedFace>,
    /// Warnings raised while resolving (raw-id references).
    pub warnings: Vec<Warning>,
}

/// Explicit tolerances one selector evaluates with.
#[derive(Debug, Clone, Copy)]
struct Tol {
    cos_angle: f64,
    distance: f64,
}

impl Tol {
    fn new(default: &Tolerance, angle_deg: Option<f64>, distance: Option<f64>) -> Tol {
        let deg = angle_deg.unwrap_or(default.angle_deg);
        Tol {
            cos_angle: (deg * std::f64::consts::PI / 180.0).pcos(),
            distance: distance.unwrap_or(default.distance),
        }
    }
}

fn surface_matches(selector: &Selector, surface: &Surface, default: &Tolerance) -> bool {
    match (selector, surface) {
        (
            Selector::Plane {
                normal,
                through,
                angle_tol,
                dist_tol,
            },
            Surface::Plane { normal: n, offset },
        ) => {
            let tol = Tol::new(default, *angle_tol, *dist_tol);
            unit(*normal).is_some_and(|want| dot(want, *n) >= tol.cos_angle)
                && (dot(*n, *through) - offset).abs() <= tol.distance
        }
        (
            Selector::Cylinder {
                axis,
                radius,
                through,
                angle_tol,
                dist_tol,
            },
            Surface::Cylinder {
                origin,
                axis: a,
                radius: r,
            },
        ) => {
            let tol = Tol::new(default, *angle_tol, *dist_tol);
            let parallel = unit(*axis).is_some_and(|want| dot(want, *a).abs() >= tol.cos_angle);
            let on_axis =
                through.is_none_or(|p| point_line_distance(p, *origin, *a) <= tol.distance);
            parallel && on_axis && (r - radius).abs() <= tol.distance
        }
        (
            Selector::Sphere {
                centre,
                radius,
                dist_tol,
            },
            Surface::Sphere {
                centre: c,
                radius: r,
            },
        ) => {
            let tol = Tol::new(default, None, *dist_tol);
            norm(sub(*centre, *c)) <= tol.distance && (r - radius).abs() <= tol.distance
        }
        _ => false,
    }
}

/// Distance from `point` to the nearest triangle of each face, indexed
/// by `FaceId`.
fn face_distances(model: &AnalysisModel, point: [f64; 3]) -> Vec<f64> {
    let mut best = vec![f64::INFINITY; model.faces.len()];
    for solid in &model.solids {
        for (corners, face) in solid.triangles.iter().zip(&solid.triangle_face) {
            let [a, b, c] = corners.map(|v| solid.vertices[v as usize]);
            let d = point_triangle_distance(point, a, b, c);
            if let Some(slot) = best.get_mut(face.0 as usize)
                && d < *slot
            {
                *slot = d;
            }
        }
    }
    best
}

fn select(model: &AnalysisModel, selector: &Selector, default: &Tolerance) -> BTreeSet<FaceId> {
    match selector {
        Selector::NearestTo(point) => {
            let d = face_distances(model, *point);
            let min = d.iter().copied().fold(f64::INFINITY, f64::min);
            if !min.is_finite() {
                return BTreeSet::new();
            }
            // Faces within the distance tolerance of the nearest tie.
            model
                .faces
                .iter()
                .filter(|f| d[f.id.0 as usize] <= min + default.distance)
                .map(|f| f.id)
                .collect()
        }
        Selector::AllOf(parts) => {
            let mut sets = parts.iter().map(|p| select(model, p, default));
            let first = sets.next().unwrap_or_default();
            sets.fold(first, |acc, s| acc.intersection(&s).copied().collect())
        }
        Selector::AnyOf(parts) => parts
            .iter()
            .flat_map(|p| select(model, p, default))
            .collect(),
        _ => model
            .faces
            .iter()
            .filter(|f| surface_matches(selector, &f.surface, default))
            .map(|f| f.id)
            .collect(),
    }
}

fn expected(expect: Expect, matched: usize) -> bool {
    match expect {
        Expect::Exactly(n) => matched == n as usize,
        Expect::AtLeast(n) => matched >= n as usize,
    }
}

/// Resolves every face definition against `model`.
pub fn resolve_faces(
    model: &AnalysisModel,
    faces: &[FaceDef],
    default: &Tolerance,
) -> Result<Resolution, AnalysisError> {
    let mut resolution = Resolution {
        faces: Vec::new(),
        warnings: Vec::new(),
    };
    for def in faces {
        let matched: BTreeSet<FaceId> = match &def.reference {
            FaceRef::Name(name) => model
                .faces
                .iter()
                .filter(|f| f.name.as_deref() == Some(name.as_str()))
                .map(|f| f.id)
                .collect(),
            FaceRef::Id(id) => {
                resolution.warnings.push(Warning::FaceIdReference {
                    name: def.name.clone(),
                });
                model
                    .faces
                    .get(id.0 as usize)
                    .map(|f| f.id)
                    .into_iter()
                    .collect()
            }
            FaceRef::Select(selector) => select(model, selector, default),
        };
        let candidates: Vec<FaceId> = matched.into_iter().collect();
        if candidates.is_empty() {
            return Err(AnalysisError::SelectorUnresolved {
                name: def.name.clone(),
            });
        }
        if !expected(def.expect, candidates.len()) {
            return Err(AnalysisError::SelectorAmbiguous {
                name: def.name.clone(),
                matched: candidates.len(),
                candidates,
            });
        }
        resolution.faces.push(ResolvedFace {
            name: def.name.clone(),
            faces: candidates,
        });
    }
    Ok(resolution)
}
