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
use std::sync::Arc;

use popkern_core::geometry::Point3;
use popkern_core::tolerance::Tolerance;

use crate::assembly::{Assembly, DefinitionId};
use crate::error::AssemblyError;
use crate::interference_refine::{bounds_diagonal, faces_at, refine};
use crate::mesh_bvh::{LocalMesh, WorldMesh};

#[path = "interference_report.rs"]
mod report;
pub use report::{InterferenceOptions, InterferenceReport, PairReport, PairStatus};

/// One pair's check: its report (`None` beyond the clearance) and
/// whether the mesh resolution settled it (see `interference_certify`).
pub(crate) struct PairCheck {
    pub(crate) report: Option<PairReport>,
    pub(crate) settled: bool,
}

pub(crate) struct Occurrence {
    pub(crate) path: String,
    pub(crate) definition: DefinitionId,
    pub(crate) world: crate::placement::Placement,
    pub(crate) bodies: Vec<(
        Arc<crate::interference_occurrences::WorldBrep>,
        Arc<WorldMesh>,
        usize,
    )>,
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
        self.interference_where(options, &|_, _| true)
    }

    /// [`Assembly::interference`] over the occurrence pairs (paths)
    /// `keep` accepts: a swing checks only the pairs its motion changes.
    pub(crate) fn interference_where(
        &mut self,
        options: &InterferenceOptions,
        keep: &dyn Fn(&str, &str) -> bool,
    ) -> Result<InterferenceReport, AssemblyError> {
        let tolerance = options.tolerance.unwrap_or_else(|| {
            self.definitions
                .iter()
                .filter_map(|d| d.model())
                .map(|m| m.tolerance().linear())
                .fold(0.0, f64::max)
        });
        let selected = |path: &str| {
            options.among.is_empty()
                || options.among.iter().any(|p| {
                    path == p
                        || path
                            .strip_prefix(p.as_str())
                            .is_some_and(|rest| rest.starts_with('/'))
                })
        };
        let wanted: std::collections::HashSet<DefinitionId> = self
            .leaves()
            .iter()
            .filter(|r| r.is_part && selected(&r.path.to_string()))
            .map(|r| r.definition)
            .collect();
        let parts = self.cached_parts(&wanted)?;
        // One tessellation per definition body, kept across checks by
        // content (a sweep moves placements, not geometry).
        let mut local: HashMap<(DefinitionId, usize), Arc<LocalMesh>> = HashMap::new();
        for (id, cached) in &parts {
            for (i, part) in cached.parts.iter().enumerate() {
                let brep = part.brep();
                let size = bounds_diagonal(brep).max(1e-12);
                let chord = popkern_mesh::TessellationOptions {
                    chord_height: options.chord_fraction * size,
                    ..popkern_mesh::TessellationOptions::default()
                };
                let mesh = self
                    .meshes
                    .get_or_build(cached.contents[i], brep, &chord)
                    .map_err(|e| AssemblyError::Geometry {
                        instance: self
                            .definition(*id)
                            .map(|d| d.name.clone())
                            .unwrap_or_default(),
                        detail: format!("tessellation for the interference check: {e}"),
                    })?;
                local.insert((*id, i), mesh);
            }
        }
        // Broad phase on world bounds; only the occurrences of candidate
        // pairs are placed (`interference_occurrences`).
        let candidates_all = self.interference_candidates(&selected, &local);
        let mut candidate_pairs: Vec<(usize, usize)> = Vec::new();
        for i in 0..candidates_all.len() {
            for j in i + 1..candidates_all.len() {
                let (a, b) = (&candidates_all[i], &candidates_all[j]);
                if !keep(&a.path, &b.path) {
                    continue;
                }
                let near = a.bodies.iter().any(|(_, ba, ea)| {
                    b.bodies.iter().any(|(_, bb, eb)| {
                        ba.distance(bb) - ea - eb <= options.clearance + tolerance
                    })
                });
                if near {
                    candidate_pairs.push((i, j));
                }
            }
        }
        let mut placed: HashMap<usize, usize> = HashMap::new();
        let mut occurrences: Vec<Occurrence> = Vec::new();
        for &(i, j) in &candidate_pairs {
            for k in [i, j] {
                if let std::collections::hash_map::Entry::Vacant(e) = placed.entry(k) {
                    e.insert(occurrences.len());
                    occurrences.push(self.place_candidate(&candidates_all[k], &parts, &local));
                }
            }
        }
        let candidate_pairs: Vec<(usize, usize)> = candidate_pairs
            .iter()
            .map(|(i, j)| (placed[i], placed[j]))
            .collect();
        let mut pairs = Vec::new();
        // Finer world meshes per occurrence body and level, shared by
        // every pair that occurrence certifies against.
        let worlds: crate::interference_certify::WorldMeshes = Default::default();
        let candidates = candidate_pairs.len();
        let this: &Self = self;
        let check = |&(i, j): &(usize, usize)| -> Option<PairReport> {
            let (a, b) = (&occurrences[i], &occurrences[j]);
            let first = this.check_pair(a, b, &local, options, tolerance, false);
            if first.settled {
                return first.report;
            }
            this.certify_pair(a, b, &parts, &worlds, options, tolerance)
                .unwrap_or_else(|| {
                    this.check_pair(a, b, &local, options, tolerance, true)
                        .report
                })
        };
        pairs.extend(
            crate::parallel_map::ordered(&candidate_pairs, check)
                .into_iter()
                .flatten(),
        );
        pairs.sort_by(|x, y| {
            crate::interference_refine::rank(&x.status)
                .total_cmp(&crate::interference_refine::rank(&y.status))
        });
        Ok(InterferenceReport {
            pairs,
            occurrences: candidates_all.len(),
            candidate_pairs: candidates,
            tolerance,
            clearance: options.clearance,
        })
    }

    pub(crate) fn check_pair(
        &self,
        a: &Occurrence,
        b: &Occurrence,
        local: &HashMap<(DefinitionId, usize), Arc<LocalMesh>>,
        options: &InterferenceOptions,
        tolerance: f64,
        run_boolean: bool,
    ) -> PairCheck {
        let (clearance, volume) = (options.clearance, options.overlap_volume);
        // Overlap first: crossing meshes or containment, then the
        // boolean decides interfering versus touching.
        for (brep_a, mesh_a, _) in &a.bodies {
            for (brep_b, mesh_b, _) in &b.bodies {
                let crossing = mesh_a.crosses(mesh_b)
                    || brep_a.first_vertex().is_some_and(|p| mesh_b.contains(&p))
                    || brep_b.first_vertex().is_some_and(|p| mesh_a.contains(&p));
                if !crossing {
                    continue;
                }
                // Deep enough to prove from the meshes: no finer mesh,
                // no boolean, when the volume is not wanted.
                if !volume && crate::mesh_bvh::certain_overlap(mesh_a, mesh_b, tolerance) {
                    return PairCheck {
                        report: Some(PairReport {
                            a: a.path.clone(),
                            b: b.path.clone(),
                            status: PairStatus::Interfering {
                                volume: None,
                                decline: None,
                            },
                            closest: None,
                            faces: None,
                        }),
                        settled: true,
                    };
                }
                // Meshes crossing at a chord coarser than the gap prove
                // nothing: leave the pair unsettled for finer meshes and
                // build the overlap only where they still cross.
                if !run_boolean {
                    return PairCheck {
                        report: None,
                        settled: false,
                    };
                }
                let (Some(world_a), Some(world_b)) = (brep_a.get(), brep_b.get()) else {
                    continue;
                };
                let status = match popkern_bool::intersection_ref(
                    world_a,
                    world_b,
                    Tolerance::default().with_linear(tolerance),
                    None,
                ) {
                    // No overlap volume: the bodies are apart or in
                    // contact — the distance branch below measures which
                    // (a 0.02 mm fit reported "touching at 0" while its
                    // coarse meshes crossed, PK-747).
                    Ok(overlap) if overlap.faces.is_empty() => continue,
                    Ok(overlap) => {
                        match popkern_measure::body_volume(&overlap, &Tolerance::default()) {
                            Ok(volume) if volume.abs() > tolerance.powi(3) => {
                                PairStatus::Interfering {
                                    volume: Some(volume.abs()),
                                    decline: None,
                                }
                            }
                            // A zero-volume overlap (a degenerate sheet
                            // the boolean returns for near-contact
                            // bodies) is no overlap either: distance
                            // decides contact.
                            Ok(_) => continue,
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
                // A positive overlap volume settles the pair; a decline
                // may be a near-contact artefact of the coarse meshes.
                let settled = matches!(
                    status,
                    PairStatus::Interfering {
                        volume: Some(_),
                        ..
                    }
                );
                return PairCheck {
                    report: Some(PairReport {
                        a: a.path.clone(),
                        b: b.path.clone(),
                        status,
                        closest: None,
                        faces: None,
                    }),
                    settled,
                };
            }
        }
        // Apart: the closest body pair.
        let mut best: Option<(f64, bool, Point3, Point3, [String; 2])> = None;
        let mut lower = f64::INFINITY;
        for (brep_a, mesh_a, body_a) in &a.bodies {
            for (brep_b, mesh_b, body_b) in &b.bodies {
                let Some(c) = mesh_a.closest_within(mesh_b, clearance + tolerance) else {
                    continue;
                };
                lower = lower.min(c.lower_bound);
                // Certainly beyond the clearance: nothing to refine (the
                // pair is not reported).
                if c.lower_bound > clearance {
                    continue;
                }
                let (Some(brep_a), Some(brep_b)) = (brep_a.get(), brep_b.get()) else {
                    continue;
                };
                let (Some(la), Some(lb)) = (
                    local.get(&(a.definition, *body_a)),
                    local.get(&(b.definition, *body_b)),
                ) else {
                    continue;
                };
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
        let Some((distance, refined, here, there, names)) = best else {
            return PairCheck {
                report: None,
                settled: true,
            };
        };
        let lower_bound = lower.min(distance);
        if lower_bound > clearance {
            return PairCheck {
                report: None,
                settled: true,
            };
        }
        // Contact is certain only on refined surface points within the
        // tolerance; a mesh distance (crossing coarse meshes measure 0)
        // or a refined gap above it with a loose lower bound is not.
        let settled = lower_bound > tolerance || (refined && distance <= tolerance);
        let status = if lower_bound <= tolerance {
            PairStatus::Touching { distance }
        } else {
            PairStatus::Clear {
                distance,
                lower_bound,
                refined,
            }
        };
        PairCheck {
            report: Some(PairReport {
                a: a.path.clone(),
                b: b.path.clone(),
                status,
                closest: Some([[here.x, here.y, here.z], [there.x, there.y, there.z]]),
                faces: Some(names),
            }),
            settled,
        }
    }
}
