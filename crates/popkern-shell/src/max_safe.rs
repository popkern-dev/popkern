//! Pre-flight maximum-safe-offset query (PK-328 row 481, design §4):
//! the largest offset magnitude for which no face vanishes, no edge
//! collapses, and no self-intersection trimming is required — the
//! slider upper bound a GUI reads before running the operation.
//!
//! Three bounds, min-combined, no surface offsetting performed:
//!
//! 1. **Surface vanish** — per face, the distance at which the offset
//!    surface degenerates. Analytic classes are exact and mirror
//!    `offset_analytic_classified`'s existence rules (cylinder/sphere
//!    radius through zero, torus minor through zero or spindle
//!    threshold, cone apex overtaking the face's axial window);
//!    free-form faces sample the offset-map fold quadratic
//!    `((Sᵤ + m·d·nᵤ) × (Sᵥ + m·d·nᵥ))·n = 0` over the trimmed UV
//!    domain — no shape-operator sign conventions to trust.
//! 2. **Edge collapse** — the linear vertex-motion model of the chord
//!    gate ([`crate::collapse_check`]) run pre-flight: per-vertex
//!    velocity from the least-squares solution of `nᵢ·v = mᵢ` over
//!    adjacent-face normals (exact for plane triple points,
//!    first-order otherwise); an edge whose chord shrinks at rate
//!    `ℓ̇ < 0` collapses at `ℓ₀/(−ℓ̇)`.
//! 3. **Thin-wall clearance** — non-adjacent face pairs approach at a
//!    closing speed `s = (v_a − v_b)·û` per unit distance; the wall
//!    closes at `gap/s`. Gaps are measured against TRIMMED patches
//!    ([`crate::body_distance`]), so faraway surface extensions cannot
//!    fake a collision.
//!
//! The reported maximum is `0.99 ×` the tightest bound (the
//! fillet-offset convention: the reported value must itself serve);
//! the acceptance gate is slider consistency on the corpus — serve at
//! the reported max, decline or trim at `1.05 ×` it.

use popkern_math::PortableMath;
use std::collections::HashMap;
use std::fmt;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::body_distance::FaceTrimEntry;
use crate::classify::classify_faces;
use crate::error::ShellError;
use crate::max_safe_motion::{thin_wall_candidates, vertex_velocities};
use crate::types::ShellDirection;
use popkern_blend::offset::analytic::radial_normal_sign;

/// Fraction of the tightest exact bound reported as the safe maximum
/// (the fillet-offset convention: the reported value must itself
/// serve).
const SAFE_FRACTION: f64 = 0.99;

/// UV grid resolution per axis for the free-form fold sampling.
const FOLD_GRID: usize = 24;

/// The mechanism that limits the offset distance, with the binding
/// entity as a typed handle.
#[derive(Debug, Clone, Copy, PartialEq)]
#[non_exhaustive]
pub enum OffsetConstraint {
    /// Nothing binds: every face can offset arbitrarily far in the
    /// requested direction (e.g. a convex plane-bounded solid growing
    /// outward).
    Unbounded,
    /// The face's offset surface degenerates (a radius driven through
    /// zero, a torus spindling, a cone apex overtaking the face, or a
    /// free-form offset map folding).
    SurfaceVanish {
        /// The binding face in the input body.
        face: FaceHandle,
    },
    /// The edge's chord shrinks to zero length under the linear
    /// vertex-motion model.
    EdgeCollapse {
        /// The binding edge in the input body.
        edge: EdgeHandle,
    },
    /// Two non-adjacent faces approach each other and their wall
    /// closes.
    ThinWall {
        /// First face of the closing pair.
        face_a: FaceHandle,
        /// Second face of the closing pair.
        face_b: FaceHandle,
    },
}

impl fmt::Display for OffsetConstraint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            OffsetConstraint::Unbounded => write!(f, "no geometric limit"),
            OffsetConstraint::SurfaceVanish { face } => {
                write!(
                    f,
                    "the surface of face {face:?} degenerates at this distance"
                )
            }
            OffsetConstraint::EdgeCollapse { edge } => {
                write!(f, "edge {edge:?} collapses to zero length")
            }
            OffsetConstraint::ThinWall { face_a, face_b } => {
                write!(
                    f,
                    "the wall between face {face_a:?} and face {face_b:?} closes"
                )
            }
        }
    }
}

/// Result of a [`max_safe_offset`] query: the slider upper bound and
/// the mechanism/entity that binds it.
#[derive(Debug, Clone, Copy)]
pub struct SafeOffsetReport {
    /// Largest safe offset magnitude ([`f64::INFINITY`] when
    /// unbounded). Already scaled by the safe fraction — the full
    /// operation serves at this value.
    pub max_distance: f64,
    /// What limits the offset.
    pub binding: OffsetConstraint,
}

impl fmt::Display for SafeOffsetReport {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.max_distance.is_finite() {
            write!(
                f,
                "maximum safe offset {:.6}: {}",
                self.max_distance, self.binding
            )
        } else {
            write!(f, "offset is unlimited: {}", self.binding)
        }
    }
}

/// Upgrade a mechanism decline to [`ShellError::OffsetTooLarge`]
/// carrying the pre-flight bound (PK-328): when the requested distance
/// provably exceeds the safe maximum, "offset 1.0 exceeds the maximum
/// 0.894 (the surface of face … degenerates)" is the GUI-readable
/// verdict — not the trim or join internals (a spindle-torus join
/// decline past the bound IS an offset-too-large condition). Declines
/// within the safe bound (genuine mechanism gaps) pass through
/// unchanged.
pub(crate) fn upgrade_trim_decline(
    brep: &BRep,
    direction: ShellDirection,
    requested: f64,
    tolerance: Tolerance,
    policy: crate::types::TrimPolicy,
    err: ShellError,
) -> ShellError {
    // Under TrimPolicy::Decline the caller asked for the raw
    // self-intersection reports — never rewrite those.
    if matches!(policy, crate::types::TrimPolicy::Decline) {
        return err;
    }
    if matches!(
        err,
        ShellError::SelfIntersecting { .. } | ShellError::JoinFailed { .. }
    ) && let Ok(report) = max_safe_offset(brep, direction, tolerance)
        && report.max_distance.is_finite()
        && requested.abs() > report.max_distance
    {
        return ShellError::OffsetTooLarge {
            requested: requested.abs(),
            max_safe: report.max_distance,
            binding: report.binding.to_string(),
        };
    }
    err
}

/// Compute the largest safe offset magnitude for `direction`
/// (design §4; interactive — no surface offsetting is performed).
///
/// `Midplane` moves each skin by half the distance, so its report is
/// twice the tighter of the two per-side bounds.
///
/// # Errors
///
/// [`ShellError::NoFaces`] on an empty body;
/// [`ShellError::MaterialSenseUnresolved`] when a component's material
/// side cannot be measured (offset direction would be meaningless).
pub fn max_safe_offset(
    brep: &BRep,
    direction: ShellDirection,
    tolerance: Tolerance,
) -> Result<SafeOffsetReport, ShellError> {
    if matches!(direction, ShellDirection::Midplane) {
        let inward = max_safe_offset(brep, ShellDirection::Inward, tolerance)?;
        let outward = max_safe_offset(brep, ShellDirection::Outward, tolerance)?;
        let side = if inward.max_distance <= outward.max_distance {
            inward
        } else {
            outward
        };
        return Ok(SafeOffsetReport {
            max_distance: side.max_distance * 2.0,
            binding: side.binding,
        });
    }

    // Unit-rate classification: only the per-face geometric sign is
    // consumed (motion = m·n per unit distance).
    let classification = classify_faces(brep, 1.0, &[], direction, &[])?;
    let motions: HashMap<FaceHandle, f64> = classification
        .offset_faces
        .iter()
        .map(|&(fh, _, inward)| (fh, if inward { -1.0 } else { 1.0 }))
        .collect();

    let mut candidates: Vec<(f64, OffsetConstraint)> = Vec::new();

    // Bound 1: surface vanish, per face in sorted order.
    let mut face_handles: Vec<FaceHandle> = motions.keys().copied().collect();
    face_handles.sort();
    for &fh in &face_handles {
        let m = motions[&fh];
        if let Some(d) = surface_vanish_bound(brep, fh, m) {
            candidates.push((d, OffsetConstraint::SurfaceVanish { face: fh }));
        }
    }

    // Bound 2: edge collapse from linear vertex motion, with a
    // curved-aware secant refinement of the near-binding candidates
    // (the linear rate is first-order at curved corners).
    let adjacency = crate::max_safe_motion::vertex_face_adjacency(brep);
    let velocities = vertex_velocities(brep, &motions);
    let mut edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    edge_handles.sort();
    let lin_tol = tolerance.linear();
    let mut edge_candidates: Vec<(
        f64,
        EdgeHandle,
        popkern_core::topology::VertexHandle,
        popkern_core::topology::VertexHandle,
    )> = Vec::new();
    for eh in edge_handles {
        let edge = &brep.edges[eh];
        // Free rims (open sheets) may carry twin-less half-edges;
        // defensive access — an edge without both endpoints simply
        // contributes no chord bound.
        let Some(he) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(he.twin) else {
            continue;
        };
        let (Some(v0d), Some(v1d)) = (brep.vertices.get(he.origin), brep.vertices.get(twin.origin))
        else {
            continue;
        };
        let p0 = v0d.position;
        let p1 = v1d.position;
        let chord = p1 - p0;
        let len0 = chord.length();
        if len0 <= lin_tol {
            continue; // closed edge: no chord direction
        }
        let dir = chord * (1.0 / len0);
        let (Some(v0), Some(v1)) = (velocities.get(&he.origin), velocities.get(&twin.origin))
        else {
            continue;
        };
        let rate = (*v1 - *v0).dot(&dir);
        if rate < -1.0e-9 {
            edge_candidates.push((len0 / (-rate), eh, he.origin, twin.origin));
        }
    }
    let min_edge_linear = edge_candidates
        .iter()
        .map(|(d, ..)| *d)
        .fold(f64::INFINITY, f64::min);
    for (d, eh, v0, v1) in edge_candidates {
        let refined = if d <= min_edge_linear * 1.3 {
            crate::max_safe_motion::refine_edge_collapse_bound(
                brep, &adjacency, &motions, v0, v1, d,
            )
            .unwrap_or(d)
        } else {
            d
        };
        candidates.push((refined, OffsetConstraint::EdgeCollapse { edge: eh }));
    }

    // Bound 3: thin-wall clearance between non-adjacent faces.
    thin_wall_candidates(brep, &motions, lin_tol, &mut candidates);

    let best = candidates
        .into_iter()
        .filter(|(d, _)| d.is_finite() && *d > 0.0)
        .fold(None::<(f64, OffsetConstraint)>, |best, cand| match best {
            Some((bd, _)) if bd <= cand.0 => best,
            _ => Some(cand),
        });
    Ok(match best {
        Some((d, binding)) => SafeOffsetReport {
            max_distance: d * SAFE_FRACTION,
            binding,
        },
        None => SafeOffsetReport {
            max_distance: f64::INFINITY,
            binding: OffsetConstraint::Unbounded,
        },
    })
}

/// Distance at which the face's offset surface degenerates in the
/// motion direction (`m·n` per unit distance). `None` = unbounded.
/// Analytic rules mirror `offset_analytic_classified`.
fn surface_vanish_bound(brep: &BRep, fh: FaceHandle, m: f64) -> Option<f64> {
    let surface = &brep.faces.get(fh)?.surface.as_ref()?.geometry;
    match surface {
        SurfaceGeometry::Plane { .. } => None,
        SurfaceGeometry::Cylinder { radius, .. } | SurfaceGeometry::Sphere { radius, .. } => {
            let radial = radial_normal_sign(surface)?;
            (m * radial < 0.0).then_some(*radius)
        }
        SurfaceGeometry::Torus {
            major_radius,
            minor_radius,
            ..
        } => {
            let radial = radial_normal_sign(surface)?;
            if m * radial < 0.0 {
                Some(*minor_radius) // tube radius through zero
            } else {
                Some(major_radius - minor_radius) // spindle threshold: trimming begins
            }
        }
        SurfaceGeometry::Cone {
            apex,
            axis,
            half_angle,
        } => {
            // The offset cone's apex translates −m·d/sin α along the
            // unit axis; for m < 0 it advances INTO the face's axial
            // window and the face fully degenerates when the apex
            // passes the window's far edge (`z_far · sin α` — the
            // exact bound the cone-through-apex gate names).
            if m >= 0.0 {
                return None;
            }
            let sin_a = half_angle.psin();
            let axis_n = axis.normalized()?;
            if sin_a <= 0.0 {
                return None;
            }
            let z_far = face_boundary_positions(brep, fh)
                .into_iter()
                .map(|p| (p - *apex).dot(&axis_n))
                .fold(f64::NEG_INFINITY, f64::max);
            (z_far > 0.0).then_some(z_far * sin_a)
        }
        SurfaceGeometry::Nurbs(_)
        | SurfaceGeometry::Gregory(_)
        | SurfaceGeometry::PolyhedralSpline(_) => free_form_fold_bound(brep, fh, surface, m),
    }
}

/// Outer-loop vertex positions of a face.
fn face_boundary_positions(brep: &BRep, fh: FaceHandle) -> Vec<Point3> {
    let Some(face) = brep.faces.get(fh) else {
        return Vec::new();
    };
    let Some(loop_data) = brep.loops.get(face.outer_loop) else {
        return Vec::new();
    };
    let start = loop_data.half_edge;
    let mut out = Vec::new();
    let mut cursor = start;
    let bound = brep.half_edges.len();
    let mut steps = 0usize;
    while let Some(he) = brep.half_edges.get(cursor) {
        if let Some(v) = brep.vertices.get(he.origin) {
            out.push(v.position);
        }
        cursor = he.next;
        steps += 1;
        if cursor == start || steps > bound {
            break;
        }
    }
    out
}

/// Smallest positive fold distance of the free-form offset map over
/// the face's trimmed UV domain: at each in-trim grid sample the area
/// form `((Sᵤ + m·d·nᵤ) × (Sᵥ + m·d·nᵥ))·n` is quadratic in `d`; its
/// smallest positive root is the exact per-sample fold distance.
fn free_form_fold_bound(
    brep: &BRep,
    fh: FaceHandle,
    surface: &SurfaceGeometry,
    m: f64,
) -> Option<f64> {
    let entry = FaceTrimEntry::build(brep, fh)?;
    let (u0, v0, u1, v1) = entry.uv_bbox()?;
    let (hu, hv) = (
        (u1 - u0).abs().max(1e-9) * 1e-4,
        (v1 - v0).abs().max(1e-9) * 1e-4,
    );
    let mut best: Option<f64> = None;
    for i in 0..=FOLD_GRID {
        for j in 0..=FOLD_GRID {
            let u = u0 + (u1 - u0) * i as f64 / FOLD_GRID as f64;
            let v = v0 + (v1 - v0) * j as f64 / FOLD_GRID as f64;
            if !entry.contains(u, v) {
                continue;
            }
            let su = surface.partial_u(u, v);
            let sv = surface.partial_v(u, v);
            let n = surface.normal_at(u, v);
            if n.length() < 1e-12 {
                continue;
            }
            let nu = (surface.normal_at(u + hu, v) - surface.normal_at(u - hu, v)) * (0.5 / hu);
            let nv = (surface.normal_at(u, v + hv) - surface.normal_at(u, v - hv)) * (0.5 / hv);
            let c0 = su.cross(&sv).dot(&n);
            let c1 = m * (su.cross(&nv).dot(&n) + nu.cross(&sv).dot(&n));
            let c2 = m * m * nu.cross(&nv).dot(&n);
            if let Some(d) = smallest_positive_root(c0, c1, c2)
                && best.is_none_or(|b| d < b)
            {
                best = Some(d);
            }
        }
    }
    best
}

/// Smallest positive root of `c0 + c1·d + c2·d² = 0` (with `c0 > 0`
/// the un-offset area form; degenerate leading terms fall back to the
/// linear root).
fn smallest_positive_root(c0: f64, c1: f64, c2: f64) -> Option<f64> {
    if c2.abs() < 1e-14 * (c0.abs() + c1.abs()).max(1.0) {
        if c1.abs() < 1e-14 {
            return None;
        }
        let r = -c0 / c1;
        return (r > 0.0).then_some(r);
    }
    let disc = c1 * c1 - 4.0 * c2 * c0;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    let (r1, r2) = ((-c1 - sq) / (2.0 * c2), (-c1 + sq) / (2.0 * c2));
    [r1, r2]
        .into_iter()
        .filter(|r| *r > 0.0)
        .fold(None, |best: Option<f64>, r| {
            Some(best.map_or(r, |b| b.min(r)))
        })
}

#[cfg(test)]
#[path = "max_safe_tests.rs"]
mod tests;
