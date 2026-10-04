//! Interference and clearance between part occurrences (design
//! `docs/architecture/assemblies.md` §8.1, PK-630).
//!
//! For every pair of part occurrences whose world bounds come within
//! the requested clearance:
//!
//! 1. **Mesh distance with a certified lower bound.** Each definition
//!    body is tessellated once and every triangle's deviation from its
//!    true face is measured ([`crate::mesh_bvh::deviation`]); a BVH
//!    branch-and-bound finds the closest triangle pair and a lower
//!    bound `min (d(ta, tb) − eps(ta) − eps(tb))` on the true distance.
//! 2. **Refinement onto the true faces.** From the closest mesh points,
//!    alternating projection onto the two faces' surfaces; when both
//!    projected points lie on their (trimmed) faces — within their
//!    face's mesh deviation of its triangles — the pair is a feasible
//!    pair of surface points, so its distance is an UPPER bound and the
//!    reported distance. Otherwise the mesh distance is reported and
//!    `refined` says so.
//! 3. **Interference.** Crossing triangles, or one body's vertex inside
//!    the other (containment), make the pair a candidate; the
//!    intersection solid is built by the boolean and its volume
//!    measured. A positive volume is `Interfering`; an empty one is
//!    `Touching`. When the boolean declines, the pair is `Interfering`
//!    with no volume and the decline text — never an estimate.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::uv_projection::project_point_to_uv;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};
use serde::Serialize;

use crate::assembly::{Assembly, DefinitionId};
use crate::error::AssemblyError;
use crate::mesh_bvh::{LocalMesh, WorldMesh};

/// What to check.
#[derive(Debug, Clone, Copy)]
pub struct InterferenceOptions {
    /// Report pairs whose distance may be below this (0 reports only
    /// touching and interfering pairs).
    pub clearance: f64,
    /// Distances at or below this are contact; `None` takes the
    /// largest build tolerance among the parts.
    pub tolerance: Option<f64>,
    /// Tessellation chord height as a fraction of each part's size.
    pub chord_fraction: f64,
}

impl Default for InterferenceOptions {
    fn default() -> Self {
        InterferenceOptions {
            clearance: 0.0,
            tolerance: None,
            chord_fraction: 1e-3,
        }
    }
}

/// How a pair of occurrences relate.
#[derive(Debug, Clone, PartialEq, Serialize)]
#[serde(tag = "status", rename_all = "snake_case")]
pub enum PairStatus {
    /// Apart: `distance` (upper bound when `refined`) and a certified
    /// `lower_bound` above the contact tolerance.
    Clear {
        /// The reported distance.
        distance: f64,
        /// Certified lower bound on the true distance.
        lower_bound: f64,
        /// Whether `distance` was refined onto the true faces.
        refined: bool,
    },
    /// In contact (distance within the tolerance), no overlap volume.
    Touching {
        /// The reported distance.
        distance: f64,
    },
    /// Overlapping solids.
    Interfering {
        /// Volume of the overlap, when the boolean built it.
        volume: Option<f64>,
        /// Why the boolean declined, when it did.
        decline: Option<String>,
    },
}

/// One reported pair.
#[derive(Debug, Clone, Serialize)]
pub struct PairReport {
    /// First occurrence path.
    pub a: String,
    /// Second occurrence path.
    pub b: String,
    /// The relation.
    #[serde(flatten)]
    pub status: PairStatus,
    /// The closest points (world), when apart or touching.
    pub closest: Option<[[f64; 3]; 2]>,
    /// The faces those points lie on, as entity paths.
    pub faces: Option<[String; 2]>,
}

/// Every reported pair, worst first (interfering, touching, then by
/// distance).
#[derive(Debug, Clone, Serialize)]
pub struct InterferenceReport {
    /// The pairs within the clearance.
    pub pairs: Vec<PairReport>,
    /// Part occurrences checked.
    pub occurrences: usize,
    /// Pairs whose bounds came within the clearance (each measured).
    pub candidate_pairs: usize,
    /// The contact tolerance used.
    pub tolerance: f64,
    /// The clearance asked for.
    pub clearance: f64,
}

struct Occurrence {
    path: String,
    definition: DefinitionId,
    bodies: Vec<(BRep, WorldMesh, usize)>,
}

impl Assembly {
    /// Check every pair of part occurrences (visible or not).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Script`] / [`AssemblyError::Geometry`] when a
    /// part cannot be built, placed or tessellated.
    pub fn interference(
        &mut self,
        options: &InterferenceOptions,
    ) -> Result<InterferenceReport, AssemblyError> {
        let tolerance = options.tolerance.unwrap_or_else(|| {
            self.part_models_mut()
                .into_iter()
                .map(|(_, m)| m.tolerance().linear())
                .fold(0.0, f64::max)
        });
        let parts = self.definition_parts()?;
        // One tessellation per definition body, measured once.
        let mut local: HashMap<(DefinitionId, usize), LocalMesh> = HashMap::new();
        for (id, bodies) in &parts {
            for (i, part) in bodies.iter().enumerate() {
                let brep = part.brep();
                let size = bounds_diagonal(brep).max(1e-12);
                let chord = popkern_mesh::TessellationOptions {
                    chord_height: options.chord_fraction * size,
                    ..popkern_mesh::TessellationOptions::default()
                };
                let mesh = popkern_mesh::tessellate(brep, &chord).map_err(|e| {
                    AssemblyError::Geometry {
                        instance: self
                            .definition(*id)
                            .map(|d| d.name.clone())
                            .unwrap_or_default(),
                        detail: format!("tessellation for the interference check: {e}"),
                    }
                })?;
                local.insert((*id, i), LocalMesh::new(brep, &mesh));
            }
        }
        let mut occurrences: Vec<Occurrence> = Vec::new();
        for body in self.flatten()? {
            let path = body.path.to_string();
            let Some(mesh) = local.get(&(body.definition, body.body)) else {
                continue;
            };
            let world = WorldMesh::new(mesh, &body.world);
            match occurrences.last_mut() {
                Some(o) if o.path == path => o.bodies.push((body.brep, world, body.body)),
                _ => occurrences.push(Occurrence {
                    path,
                    definition: body.definition,
                    bodies: vec![(body.brep, world, body.body)],
                }),
            }
        }
        let mut pairs = Vec::new();
        let mut candidates = 0;
        for i in 0..occurrences.len() {
            for j in i + 1..occurrences.len() {
                let (a, b) = (&occurrences[i], &occurrences[j]);
                let near = a.bodies.iter().any(|(_, ma, _)| {
                    b.bodies.iter().any(|(_, mb, _)| {
                        ma.bounds.distance(&mb.bounds) - ma.max_eps() - mb.max_eps()
                            <= options.clearance + tolerance
                    })
                });
                if !near {
                    continue;
                }
                candidates += 1;
                if let Some(report) = self.check_pair(a, b, &local, tolerance, options.clearance) {
                    pairs.push(report);
                }
            }
        }
        pairs.sort_by(|x, y| rank(&x.status).total_cmp(&rank(&y.status)));
        Ok(InterferenceReport {
            pairs,
            occurrences: occurrences.len(),
            candidate_pairs: candidates,
            tolerance,
            clearance: options.clearance,
        })
    }

    fn check_pair(
        &self,
        a: &Occurrence,
        b: &Occurrence,
        local: &HashMap<(DefinitionId, usize), LocalMesh>,
        tolerance: f64,
        clearance: f64,
    ) -> Option<PairReport> {
        // Overlap first: crossing meshes or containment, then the
        // boolean decides interfering versus touching.
        for (brep_a, mesh_a, _) in &a.bodies {
            for (brep_b, mesh_b, _) in &b.bodies {
                let crossing = mesh_a.crosses(mesh_b)
                    || brep_a
                        .vertices
                        .values()
                        .next()
                        .is_some_and(|v| mesh_b.contains(&v.position))
                    || brep_b
                        .vertices
                        .values()
                        .next()
                        .is_some_and(|v| mesh_a.contains(&v.position));
                if !crossing {
                    continue;
                }
                let status = match popkern_bool::intersection_ref(
                    brep_a,
                    brep_b,
                    Tolerance::default().with_linear(tolerance),
                    None,
                ) {
                    Ok(overlap) if overlap.faces.is_empty() => {
                        PairStatus::Touching { distance: 0.0 }
                    }
                    Ok(overlap) => {
                        match popkern_measure::body_volume(&overlap, &Tolerance::default()) {
                            Ok(volume) if volume.abs() > tolerance.powi(3) => {
                                PairStatus::Interfering {
                                    volume: Some(volume.abs()),
                                    decline: None,
                                }
                            }
                            Ok(_) => PairStatus::Touching { distance: 0.0 },
                            Err(e) => PairStatus::Interfering {
                                volume: None,
                                decline: Some(format!("overlap volume: {e}")),
                            },
                        }
                    }
                    Err(e) => PairStatus::Interfering {
                        volume: None,
                        decline: Some(e.to_string()),
                    },
                };
                return Some(PairReport {
                    a: a.path.clone(),
                    b: b.path.clone(),
                    status,
                    closest: None,
                    faces: None,
                });
            }
        }
        // Apart: the closest body pair.
        let mut best: Option<(f64, bool, Point3, Point3, [String; 2])> = None;
        let mut lower = f64::INFINITY;
        for (brep_a, mesh_a, body_a) in &a.bodies {
            for (brep_b, mesh_b, body_b) in &b.bodies {
                let Some(c) = mesh_a.closest(mesh_b) else {
                    continue;
                };
                lower = lower.min(c.lower_bound);
                let la = local.get(&(a.definition, *body_a))?;
                let lb = local.get(&(b.definition, *body_b))?;
                let fa = la.faces[mesh_a.local[c.tri_here]];
                let fb = lb.faces[mesh_b.local[c.tri_there]];
                // Mesh points on triangles that lie ON their faces (zero
                // measured deviation: planar facets) are already feasible
                // surface points; others are projected onto the faces.
                let size = mesh_a.bounds.diagonal().max(mesh_b.bounds.diagonal());
                let exact = la.eps[mesh_a.local[c.tri_here]] <= 1e-12 * size
                    && lb.eps[mesh_b.local[c.tri_there]] <= 1e-12 * size;
                // A closest mesh point on a mesh edge or vertex belongs to
                // every face meeting there: refinement tries each pair and
                // keeps the nearest feasible one (coaxial pin and bore tie
                // along the bore's rim, where the plate's top face cannot
                // carry the clearance and the bore wall can; PK-651).
                let refined = if exact {
                    Some((c.here, c.there, fa, fb))
                } else {
                    let faces_a = faces_at(mesh_a, la, &c.here, fa, size);
                    let faces_b = faces_at(mesh_b, lb, &c.there, fb, size);
                    faces_a
                        .iter()
                        .flat_map(|&x| faces_b.iter().map(move |&y| (x, y)))
                        .filter_map(|(x, y)| {
                            refine(
                                brep_a, x, mesh_a, la, brep_b, y, mesh_b, lb, c.here, c.there,
                            )
                            .map(|(p, q)| (p, q, x, y))
                        })
                        .min_by(|l, r| (l.0 - l.1).length().total_cmp(&(r.0 - r.1).length()))
                };
                let (distance, here, there, is_refined, fa, fb) = match refined {
                    Some((p, q, x, y)) => ((p - q).length(), p, q, true, x, y),
                    None => (c.distance, c.here, c.there, false, fa, fb),
                };
                if best.as_ref().is_none_or(|b| distance < b.0) {
                    let names = [
                        self.face_path(&a.path, a.definition, brep_a, fa),
                        self.face_path(&b.path, b.definition, brep_b, fb),
                    ];
                    best = Some((distance, is_refined, here, there, names));
                }
            }
        }
        let (distance, refined, here, there, names) = best?;
        let lower_bound = lower.min(distance);
        if lower_bound > clearance {
            return None;
        }
        let status = if lower_bound <= tolerance {
            PairStatus::Touching { distance }
        } else {
            PairStatus::Clear {
                distance,
                lower_bound,
                refined,
            }
        };
        Some(PairReport {
            a: a.path.clone(),
            b: b.path.clone(),
            status,
            closest: Some([[here.x, here.y, here.z], [there.x, there.y, there.z]]),
            faces: Some(names),
        })
    }

    /// `<occurrence>:<face's stable name>`, or `<occurrence>:face` when
    /// the face carries no name.
    fn face_path(
        &self,
        path: &str,
        definition: DefinitionId,
        brep: &BRep,
        face: FaceHandle,
    ) -> String {
        let spelled = self
            .definition(definition)
            .ok()
            .and_then(|d| d.model())
            .and_then(|model| {
                let name = brep.faces.get(face)?.name?;
                let entity = brep.name_pool.get(name)?;
                model.spell_entity(entity).ok()
            })
            .unwrap_or_else(|| "face".to_string());
        format!("{path}:{spelled}")
    }
}

/// Alternating projection of the closest mesh points onto the true
/// faces; `Some` only when both end points lie on their trimmed faces.
#[allow(clippy::too_many_arguments)]
fn refine(
    brep_a: &BRep,
    fa: FaceHandle,
    mesh_a: &WorldMesh,
    la: &LocalMesh,
    brep_b: &BRep,
    fb: FaceHandle,
    mesh_b: &WorldMesh,
    lb: &LocalMesh,
    here: Point3,
    there: Point3,
) -> Option<(Point3, Point3)> {
    let sa = &brep_a.faces.get(fa)?.surface.as_ref()?.geometry;
    let sb = &brep_b.faces.get(fb)?.surface.as_ref()?.geometry;
    let size = mesh_a
        .bounds
        .diagonal()
        .max(mesh_b.bounds.diagonal())
        .max(1e-12);
    let project = |s: &popkern_core::geometry::surface::SurfaceGeometry, p: &Point3| {
        let (u, v) = project_point_to_uv(s, p)?;
        Some(s.point_at(u, v))
    };
    let mut p = project(sa, &here)?;
    let mut q = project(sb, &there)?;
    for _ in 0..64 {
        let q_next = project(sb, &p)?;
        let p_next = project(sa, &q_next)?;
        let moved = (p_next - p).length() + (q_next - q).length();
        p = p_next;
        q = q_next;
        if moved <= 1e-14 * size {
            break;
        }
    }
    let on_face = |mesh: &WorldMesh, local: &LocalMesh, face: FaceHandle, x: &Point3| {
        let mut nearest = f64::INFINITY;
        let mut eps: f64 = 0.0;
        for (k, tri) in mesh.tris.iter().enumerate() {
            let i = mesh.local[k];
            if local.faces[i] != face {
                continue;
            }
            eps = eps.max(local.eps[i]);
            nearest = nearest.min((crate::tri_geometry::closest_on_triangle(x, tri) - *x).length());
        }
        nearest <= eps + 1e-9 * size
    };
    (on_face(mesh_a, la, fa, &p) && on_face(mesh_b, lb, fb, &q)).then_some((p, q))
}

/// The faces whose triangles hold mesh point `x` (within `1e-9` of
/// `size`), `own` first: more than one when `x` lies on a mesh edge or
/// vertex between faces.
fn faces_at(
    mesh: &WorldMesh,
    local: &LocalMesh,
    x: &Point3,
    own: FaceHandle,
    size: f64,
) -> Vec<FaceHandle> {
    let mut faces = vec![own];
    for (k, tri) in mesh.tris.iter().enumerate() {
        let face = local.faces[mesh.local[k]];
        if !faces.contains(&face)
            && (crate::tri_geometry::closest_on_triangle(x, tri) - *x).length() <= 1e-9 * size
        {
            faces.push(face);
        }
    }
    faces
}

fn bounds_diagonal(brep: &BRep) -> f64 {
    let mut lo = [f64::INFINITY; 3];
    let mut hi = [f64::NEG_INFINITY; 3];
    for v in brep.vertices.values() {
        for (k, x) in [v.position.x, v.position.y, v.position.z]
            .into_iter()
            .enumerate()
        {
            lo[k] = lo[k].min(x);
            hi[k] = hi[k].max(x);
        }
    }
    (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt()
}

fn rank(status: &PairStatus) -> f64 {
    match status {
        PairStatus::Interfering { .. } => -2.0,
        PairStatus::Touching { .. } => -1.0,
        PairStatus::Clear { distance, .. } => *distance,
    }
}
