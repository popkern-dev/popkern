//! Resolution check of stress raisers (map §1.6 step 5b): each concave
//! fillet, and each concave cylinder without a role (bores, holes, slot
//! ends: PH-056), in a stressed region is compared with the cell size.

use serde::Serialize;

use crate::model::{AnalysisModel, FaceId, FaceInfo, FaceRole, Surface};
use crate::structural::StructuralResult;

use super::{AssessError, AssessOptions};

/// Worst measured under-read of a fillet peak by r/h (PH-027: σ₁ on the
/// cell-resolution surface samples of the C8/C9 stepped bars, worst of
/// tension and bending over 4 grid shifts). Linear between rows; below
/// the first row the first value is a lower bound.
pub const UNDER_READ: [(f64, f64); 4] = [(1.0, 0.464), (2.0, 0.237), (3.0, 0.144), (4.0, 0.067)];

/// How well a fillet is resolved.
#[derive(Debug, Clone, Copy, PartialEq, Serialize)]
pub enum FilletStatus {
    /// r/h at or above `fillet_resolved_cells`.
    Resolved,
    /// Peak under-resolved: it may read low by about `expected_low_by`
    /// (a fraction, from [`UNDER_READ`]).
    UnderResolved {
        /// Expected relative under-read.
        expected_low_by: f64,
    },
    /// r/h below `fillet_sharp_cells`: treated as a sharp corner.
    Sharp,
}

/// One concave fillet or bore in a stressed region.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FilletCheck {
    /// The fillet or bore face.
    pub face: FaceId,
    /// Its radius, mm.
    pub radius: f64,
    /// r / h.
    pub cells: f64,
    /// Largest von Mises on the face, MPa.
    pub peak_von_mises: f64,
    /// The verdict.
    pub status: FilletStatus,
}

/// Worst measured under-read of a bore's peak by r/h, the worst of
/// three sources at each row: C6, a plate with a central hole in
/// tension over 4 grid shifts against its converged peak (PH-058: 1 to
/// 8 cells; 0.552, 0.460, 0.258, 0.175, 0.121, 0.083); C10 tension, a
/// semicircular notch (PH-027: 3 to 6 cells); and bored_block's 0.18
/// bore in tension against CalculiX at one shift (PH-050: 1 to 4
/// cells; 0.613 at 1 cell sets that row). A hole in tension under-reads
/// more than a fillet at the same r/h.
pub const BORE_UNDER_READ: [(f64, f64); 6] = [
    (1.0, 0.613),
    (2.0, 0.460),
    (3.0, 0.258),
    (4.0, 0.175),
    (5.0, 0.121),
    (6.0, 0.083),
];

/// Expected under-read of a fillet at `cells` from [`UNDER_READ`].
pub fn expected_under_read(cells: f64) -> f64 {
    under_read(&UNDER_READ, cells)
}

/// `table` at `cells`: linear between rows, the end values beyond them.
fn under_read(table: &[(f64, f64)], cells: f64) -> f64 {
    let Some(&first) = table.first() else {
        return 0.0;
    };
    let mut prev = first;
    if cells <= prev.0 {
        return prev.1;
    }
    for &(c, v) in &table[1..] {
        if cells <= c {
            let t = (cells - prev.0) / (c - prev.0);
            return prev.1 + t * (v - prev.1);
        }
        prev = (c, v);
    }
    prev.1
}

/// Whether the blend `surface` is concave as seen from the material at
/// `p` with outward normal `n` (the axis lies outside the solid).
pub(super) fn concave_at(surface: &Surface, p: [f64; 3], n: [f64; 3]) -> Option<bool> {
    let (centre, axis, ring) = match *surface {
        Surface::Cylinder { origin, axis, .. } => (origin, axis, None),
        Surface::Torus {
            centre,
            axis,
            major_radius,
            ..
        } => (centre, axis, Some(major_radius)),
        _ => return None,
    };
    let v = [p[0] - centre[0], p[1] - centre[1], p[2] - centre[2]];
    let a = v[0] * axis[0] + v[1] * axis[1] + v[2] * axis[2];
    let mut r = [v[0] - a * axis[0], v[1] - a * axis[1], v[2] - a * axis[2]];
    if let Some(major) = ring {
        // From the tube's centre circle instead of the axis.
        let len = (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]).sqrt();
        if len == 0.0 {
            return None;
        }
        let c = r.map(|x| major * x / len);
        r = [v[0] - c[0], v[1] - c[1], v[2] - c[2]];
    }
    Some(r[0] * n[0] + r[1] * n[1] + r[2] * n[2] < 0.0)
}

/// Checks every concave fillet whose peak von Mises is at least
/// `stressed_fraction` of its solid's.
pub fn fillet_checks(
    model: &AnalysisModel,
    result: &StructuralResult,
    options: &AssessOptions,
) -> Result<Vec<FilletCheck>, AssessError> {
    let radius = |f: &FaceInfo| match f.role {
        Some(FaceRole::Fillet { radius }) => Some(radius),
        _ => None,
    };
    let rule = Rule {
        radius: &radius,
        resolved_cells: options.fillet_resolved_cells,
        table: &UNDER_READ,
    };
    raiser_checks(model, result, options, &rule)
}

/// Checks every concave cylinder without a role (a bore, hole or slot
/// end) whose peak von Mises is at least `stressed_fraction` of its
/// solid's, against [`BORE_UNDER_READ`].
pub fn bore_checks(
    model: &AnalysisModel,
    result: &StructuralResult,
    options: &AssessOptions,
) -> Result<Vec<FilletCheck>, AssessError> {
    let radius = |f: &FaceInfo| match (f.role, &f.surface) {
        (None, Surface::Cylinder { radius, .. }) => Some(*radius),
        _ => None,
    };
    let rule = Rule {
        radius: &radius,
        resolved_cells: options.bore_resolved_cells,
        table: &BORE_UNDER_READ,
    };
    raiser_checks(model, result, options, &rule)
}

/// Per face of `mesh` with a `radius`, the concavity vote of its
/// triangles (concave minus convex).
fn concavity_votes(
    model: &AnalysisModel,
    mesh: &crate::model::SolidMesh,
    radius: &dyn Fn(&FaceInfo) -> Option<f64>,
) -> Result<std::collections::BTreeMap<FaceId, i64>, AssessError> {
    let mut votes: std::collections::BTreeMap<FaceId, i64> = Default::default();
    for (t, &f) in mesh.triangles.iter().zip(&mesh.triangle_face) {
        let Some(info) = model.faces.get(f.0 as usize) else {
            return Err(AssessError::Mismatch("face id"));
        };
        if radius(info).is_none() {
            continue;
        }
        let [a, b, c] = t.map(|i| mesh.vertices[i as usize]);
        let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
        let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
        let n = [
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        ];
        let centroid = [0, 1, 2].map(|k| (a[k] + b[k] + c[k]) / 3.0);
        let entry = votes.entry(f).or_insert(0);
        match concave_at(&info.surface, centroid, n) {
            Some(true) => *entry += 1,
            Some(false) => *entry -= 1,
            None => {}
        }
    }
    Ok(votes)
}

/// Every concave fillet and every concave cylinder without a role (the
/// faces [`fillet_checks`] and [`bore_checks`] cover), with its radius,
/// in face order: the raisers of any field, regardless of its peaks.
pub fn concave_raisers(model: &AnalysisModel) -> Result<Vec<(FaceId, f64)>, AssessError> {
    let radius = |f: &FaceInfo| match (f.role, &f.surface) {
        (Some(FaceRole::Fillet { radius }), _) => Some(radius),
        (None, Surface::Cylinder { radius, .. }) => Some(*radius),
        _ => None,
    };
    let mut out = Vec::new();
    for mesh in &model.solids {
        for (face, vote) in concavity_votes(model, mesh, &radius)? {
            if let Some(r) = model.faces.get(face.0 as usize).and_then(radius)
                && vote > 0
            {
                out.push((face, r));
            }
        }
    }
    out.sort_by_key(|(f, _)| *f);
    Ok(out)
}

/// Which faces a check covers and how it rates them.
struct Rule<'a> {
    /// The raiser radius of a covered face.
    radius: &'a dyn Fn(&FaceInfo) -> Option<f64>,
    /// r/h from which the peak is resolved.
    resolved_cells: f64,
    /// Expected under-read by r/h.
    table: &'a [(f64, f64)],
}

fn raiser_checks(
    model: &AnalysisModel,
    result: &StructuralResult,
    options: &AssessOptions,
    rule: &Rule<'_>,
) -> Result<Vec<FilletCheck>, AssessError> {
    let h = result.report.cell_size.ok_or(AssessError::NoCellSize)?;
    if result.stresses.len() != model.solids.len() {
        return Err(AssessError::Mismatch("solid count"));
    }
    let radius_of = |f: FaceId| model.faces.get(f.0 as usize).and_then(|i| (rule.radius)(i));
    let mut out = Vec::new();
    for (mesh, stresses) in model.solids.iter().zip(&result.stresses) {
        let solid_peak = stresses
            .surface
            .iter()
            .map(|s| s.stress.von_mises)
            .fold(0.0, f64::max);
        // Per face: peak von Mises and the concavity vote.
        let votes = concavity_votes(model, mesh, rule.radius)?;
        let mut faces: std::collections::BTreeMap<FaceId, (f64, i64)> =
            votes.iter().map(|(f, v)| (*f, (0.0, *v))).collect();
        for s in &stresses.surface {
            if let Some(entry) = faces.get_mut(&s.face) {
                entry.0 = entry.0.max(s.stress.von_mises);
            } else if radius_of(s.face).is_some() {
                faces.insert(s.face, (s.stress.von_mises, 0));
            }
        }
        for (face, (peak, vote)) in faces {
            let Some(radius) = radius_of(face) else {
                continue;
            };
            if vote <= 0 || peak < options.stressed_fraction * solid_peak {
                continue;
            }
            let cells = radius / h;
            let status = if cells < options.fillet_sharp_cells {
                FilletStatus::Sharp
            } else if cells < rule.resolved_cells {
                FilletStatus::UnderResolved {
                    expected_low_by: under_read(rule.table, cells),
                }
            } else {
                FilletStatus::Resolved
            };
            out.push(FilletCheck {
                face,
                radius,
                cells,
                peak_von_mises: peak,
                status,
            });
        }
    }
    Ok(out)
}
