//! The mate solver (design §8.2): the placements of the root instances
//! the mates reference, solved together.
//!
//! - Unknowns: a 6-DoF local increment `(ω, t)` per free root instance,
//!   applied about the instance's current origin (`T(t) ∘ R(o, ω) ∘ P`).
//!   The ground (the named one, else the first root instance) never
//!   moves.
//! - Residuals: per mate, rows that vanish exactly when it holds
//!   (`mate_rows`); lengths are divided by the model size so every row
//!   is dimensionless.
//! - Steps: Gauss–Newton with the truncated-SVD pseudo-inverse (the
//!   minimum-norm step — redundant mates are harmless and unconstrained
//!   motions stay where the user placed them) and a backtracking line
//!   search. A stall short of the tolerance (a stationary point that is
//!   not a solution) retries from the start with Levenberg–Marquardt;
//!   a conflict therefore costs both before it is reported. The Jacobian is central-differenced (h = 1e-7): the
//!   acceptance test is on the residual itself, so Jacobian error only
//!   costs iterations, never accuracy.
//! - Diagnostics at the solution: free DoF per instance (the Jacobian's
//!   null space restricted to that instance's columns), redundant mates
//!   (rows adding no rank given the mates before them), and, when the
//!   residual will not vanish, the conflicting mates. A conflict leaves
//!   every placement as it was.

use popkern_math::PortableMath;
use std::collections::BTreeMap;

use popkern_core::geometry::{Transform3, Vec3};

use crate::mate_steps::max_abs;
use serde::Serialize;

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::history::AssemblyStepKind;
use crate::mates::{Mate, MateKind, Primitive, Resolved};
use crate::placement::Placement;

/// Rows within this (dimensionless) are satisfied.
pub const MATE_TOLERANCE: f64 = 1e-10;
pub(crate) const RANK_REL: f64 = 1e-9;

/// What a solve found.
#[derive(Debug, Clone, Serialize)]
pub struct MateReport {
    /// Every mate holds to [`MATE_TOLERANCE`].
    pub converged: bool,
    /// Gauss–Newton iterations taken.
    pub iterations: usize,
    /// The largest residual row at the end.
    pub residual: f64,
    /// The instance that stayed put.
    pub ground: Option<String>,
    /// Remaining free motions per instance the mates move (0: fully
    /// placed).
    pub free_dof: BTreeMap<String, usize>,
    /// Mates (indices) that add no constraint beyond the earlier ones.
    pub redundant: Vec<usize>,
    /// Mates (indices) that do not hold when the solve gives up.
    pub conflicting: Vec<usize>,
    /// The mates in words, by index.
    pub mates: Vec<String>,
    /// The clusters solved in turn, from the ground out: a root whose
    /// mates to already-placed roots fix it alone is its own cluster;
    /// coupled roots (a loop) solve together.
    pub clusters: Vec<ClusterReport>,
}

/// One cluster of the decomposed solve.
#[derive(Debug, Clone, Serialize)]
pub struct ClusterReport {
    /// Its root instances.
    pub roots: Vec<String>,
    /// Gauss–Newton iterations taken.
    pub iterations: usize,
    /// Its largest residual row at the end.
    pub residual: f64,
    /// It started at its predicted placement: carried with the roots
    /// fixing it (which all moved alike) or turned to its joint's angle.
    /// An unpredicted cluster started where it was.
    pub predicted: bool,
}

impl Assembly {
    /// The mates, in the order they were added.
    pub fn mates(&self) -> &[Mate] {
        &self.mates
    }

    /// The root instance the solve holds fixed, when one was set.
    pub fn ground(&self) -> Option<&str> {
        self.ground.as_deref()
    }

    /// Fix `name` (a root instance) as the solve's ground.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::UnknownInstance`].
    pub fn set_ground(&mut self, name: &str) -> Result<(), AssemblyError> {
        if !self.root.iter().any(|i| i.name == name) {
            return Err(AssemblyError::UnknownInstance {
                path: name.to_string(),
            });
        }
        self.ground = Some(name.to_string());
        Ok(())
    }

    /// Add `mate` and solve. A mate that conflicts with the others is
    /// removed again and nothing moves (transactional).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Mate`] for an entity that cannot be mated;
    /// [`AssemblyError::MateConflict`] (with the report) when the mates
    /// cannot all hold.
    pub fn add_mate(&mut self, mate: Mate) -> Result<MateReport, AssemblyError> {
        let before = self.root_placements();
        self.mates.push(mate);
        let solved = self.solve_in_place(true);
        match solved {
            Ok(_) => self.record_solve(AssemblyStepKind::Mate(self.mates.len() - 1), &before),
            Err(_) => {
                self.mates.pop();
            }
        }
        solved
    }

    /// Solve every mate, moving the free root instances.
    ///
    /// # Errors
    ///
    /// As [`Assembly::add_mate`].
    pub fn solve_mates(&mut self) -> Result<MateReport, AssemblyError> {
        let before = self.root_placements();
        let solved = self.solve_in_place(true)?;
        self.record_solve(AssemblyStepKind::Solve, &before);
        Ok(solved)
    }

    fn root_placements(&self) -> Vec<Placement> {
        self.root.iter().map(|i| i.placement).collect()
    }

    /// Solve every mate in place. `diagnose` adds the per-cluster rank
    /// diagnostics (free motions, redundant mates) — authoring solves
    /// want them; a slider re-solve changes no structure and skips
    /// them (they were 70% of a 32-root chain's re-solve). Conflicts
    /// are found either way, from the final residual of every mate.
    pub(crate) fn solve_in_place(&mut self, diagnose: bool) -> Result<MateReport, AssemblyError> {
        let ground = self
            .ground
            .clone()
            .or_else(|| self.root.first().map(|i| i.name.clone()));
        let mut resolved = Vec::new();
        let rows = self.tree();
        for mate in &self.mates.clone() {
            let about = match &mate.kind {
                MateKind::Angle { about: Some(r), .. } => Some(self.resolve_mate_ref_in(r, &rows)?),
                _ => None,
            };
            resolved.push((
                mate.clone(),
                self.resolve_mate_ref_in(&mate.a, &rows)?,
                self.resolve_mate_ref_in(&mate.b, &rows)?,
                about,
            ));
        }
        // The root instances each mate touches; free roots in
        // first-reference order.
        // A released joint's angle mate binds nothing (PK-808).
        let touches: Vec<Vec<usize>> = (resolved.iter().enumerate())
            .map(|(i, (_, a, b, about))| {
                if self.released.contains(&i) {
                    return Vec::new();
                }
                let mut t = vec![a.root, b.root];
                t.extend(about.as_ref().map(|r| r.root));
                t.sort_unstable();
                t.dedup();
                t
            })
            .collect();
        let mut free: Vec<usize> = Vec::new();
        for (_, a, b, about) in &resolved {
            for r in [Some(a), Some(b), about.as_ref()].into_iter().flatten() {
                let name = &self.root[r.root].name;
                if Some(name) != ground.as_ref() && !free.contains(&r.root) {
                    free.push(r.root);
                }
            }
        }
        let mut current: Vec<Placement> = self.root.iter().map(|i| i.placement).collect();
        let size = 1.0
            + resolved
                .iter()
                .flat_map(|(_, a, b, _)| [a, b])
                .map(|r| crate::mate_steps::coords_max(&r.local.placed(&current[r.root])))
                .fold(0.0, f64::max);
        let clusters = crate::mate_clusters::plan(
            &resolved,
            &touches,
            &free,
            &current,
            size,
            self.root.len(),
        )?;
        // A released joint's angle is the solve's answer, not a target.
        let joints: Vec<_> = (self.joints.iter())
            .filter(|j| !self.released.contains(&j.mates[2]))
            .cloned()
            .collect();
        let predicted = crate::mate_predict::predict(
            &joints,
            &resolved,
            &touches,
            &clusters,
            &mut current,
            size,
        );
        // The released joints' limits, for the solve to keep where a
        // target jump can choose its branch (PK-902).
        let bounds: Vec<JointBound> = (self.joints.iter())
            .filter(|j| self.released.contains(&j.mates[2]))
            .filter_map(|j| {
                let (_, za, zb, Some(axis)) = resolved.get(j.mates[2])? else {
                    return None;
                };
                let limits = self.params.iter().find(|p| p.name == j.name)?.bounds?;
                Some(JointBound {
                    za: za.clone(),
                    zb: zb.clone(),
                    axis: axis.clone(),
                    limits,
                })
            })
            .collect();
        let mut known: Vec<bool> = vec![true; self.root.len()];
        for &f in &free {
            known[f] = false;
        }
        let mut iterations = 0;
        let mut free_dof = BTreeMap::new();
        let mut redundant = Vec::new();
        let mut conflicting = Vec::new();
        let mut cluster_reports = Vec::new();
        for (cluster, &predicted) in clusters.iter().zip(&predicted) {
            let mate_ids: Vec<usize> = (0..resolved.len())
                .filter(|&i| {
                    !self.released.contains(&i)
                        && touches[i].iter().any(|r| cluster.contains(r))
                        && touches[i].iter().all(|r| known[*r] || cluster.contains(r))
                })
                .collect();
            let sub: Vec<_> = mate_ids.iter().map(|&i| resolved[i].clone()).collect();
            let problem = Problem {
                resolved: &sub,
                free: cluster,
                base: &current,
                size,
                limits: &bounds,
            };
            let (x, its) = crate::mate_steps::gauss_newton(&problem)?;
            iterations += its;
            let cluster_residual = if diagnose {
                let d = crate::mate_steps::diagnose(&problem, &x)?;
                for (k, &root) in cluster.iter().enumerate() {
                    free_dof.insert(self.root[root].name.clone(), d.free_dof[k]);
                }
                redundant.extend(d.redundant.iter().map(|&i| mate_ids[i]));
                conflicting.extend(d.conflicting.iter().map(|&i| mate_ids[i]));
                d.residual
            } else {
                max_abs(&problem.residual(&x)?)
            };
            cluster_reports.push(ClusterReport {
                roots: cluster.iter().map(|&r| self.root[r].name.clone()).collect(),
                iterations: its,
                residual: cluster_residual,
                predicted,
            });
            let moved: Vec<Placement> = cluster
                .iter()
                .enumerate()
                .map(|(k, &root)| problem.placement(&x, k, root))
                .collect();
            for (&root, p) in cluster.iter().zip(moved) {
                current[root] = p;
                known[root] = true;
            }
        }
        // Every mate at the final placements (mates among held
        // occurrences only are checked here, never solved).
        let all = Problem {
            resolved: &resolved,
            free: &[],
            base: &current,
            size,
            limits: &[],
        };
        let rows = all.rows(&[])?;
        let mut r = all.residual(&[])?;
        for (index, range) in rows.iter().enumerate() {
            if self.released.contains(&index) {
                r[range.clone()].iter_mut().for_each(|v| *v = 0.0);
            } else if r[range.clone()].iter().any(|v| v.abs() > MATE_TOLERANCE) {
                conflicting.push(index);
            }
        }
        redundant.sort_unstable();
        redundant.dedup();
        conflicting.sort_unstable();
        conflicting.dedup();
        let residual = max_abs(&r);
        let converged = residual <= MATE_TOLERANCE;
        let report = MateReport {
            converged,
            iterations,
            residual,
            ground,
            free_dof,
            redundant,
            conflicting,
            mates: self.mates.iter().map(Mate::describe).collect(),
            clusters: cluster_reports,
        };
        if !converged {
            return Err(AssemblyError::MateConflict {
                report: Box::new(report),
            });
        }
        let before = self.root_placements();
        for &root in &free {
            self.root[root].placement =
                Placement::from_matrix(*current[root].transform(), &self.root[root].name)?;
        }
        let (mates, params) = (self.mates.clone(), self.params.clone());
        if self.settle_released(&resolved, &before)? {
            // A released joint's followers turn with its solved angle:
            // solve again with their mates at it (the joint's own
            // cluster is converged, so this costs the followers' solve).
            return self.solve_in_place(diagnose).inspect_err(|_| {
                for (instance, placement) in self.root.iter_mut().zip(&before) {
                    instance.placement = *placement;
                }
                (self.mates, self.params) = (mates, params);
            });
        }
        Ok(report)
    }
}

#[path = "mate_problem.rs"]
mod problem;
pub(crate) use problem::{JointBound, Problem, limit_excess, mate_rows};
