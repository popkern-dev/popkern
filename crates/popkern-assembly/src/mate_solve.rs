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
//!   search. The Jacobian is central-differenced (h = 1e-7): the
//!   acceptance test is on the residual itself, so Jacobian error only
//!   costs iterations, never accuracy.
//! - Diagnostics at the solution: free DoF per instance (the Jacobian's
//!   null space restricted to that instance's columns), redundant mates
//!   (rows adding no rank given the mates before them), and, when the
//!   residual will not vanish, the conflicting mates. A conflict leaves
//!   every placement as it was.

use popkern_math::PortableMath;
use std::collections::BTreeMap;

use popkern_core::geometry::{Point3, Transform3, Vec3};
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
        let touches: Vec<Vec<usize>> = resolved
            .iter()
            .map(|(_, a, b, about)| {
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
                .map(|r| point_of(&r.local.placed(&current[r.root])).coords_max())
                .fold(0.0, f64::max);
        let clusters = crate::mate_clusters::plan(
            &resolved,
            &touches,
            &free,
            &current,
            size,
            self.root.len(),
        )?;
        let mut known: Vec<bool> = vec![true; self.root.len()];
        for &f in &free {
            known[f] = false;
        }
        let mut iterations = 0;
        let mut free_dof = BTreeMap::new();
        let mut redundant = Vec::new();
        let mut conflicting = Vec::new();
        let mut cluster_reports = Vec::new();
        for cluster in &clusters {
            let mate_ids: Vec<usize> = (0..resolved.len())
                .filter(|&i| {
                    touches[i].iter().any(|r| cluster.contains(r))
                        && touches[i].iter().all(|r| known[*r] || cluster.contains(r))
                })
                .collect();
            let sub: Vec<_> = mate_ids.iter().map(|&i| resolved[i].clone()).collect();
            let problem = Problem {
                resolved: &sub,
                free: cluster,
                base: &current,
                size,
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
        };
        let rows = all.rows(&[])?;
        let r = all.residual(&[])?;
        for (index, range) in rows.iter().enumerate() {
            if r[range.clone()].iter().any(|v| v.abs() > MATE_TOLERANCE) {
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
        for &root in &free {
            self.root[root].placement =
                Placement::from_matrix(*current[root].transform(), &self.root[root].name)?;
        }
        Ok(report)
    }
}

pub(crate) struct Problem<'a> {
    pub(crate) resolved: &'a [(Mate, Resolved, Resolved, Option<Resolved>)],
    pub(crate) free: &'a [usize],
    pub(crate) base: &'a [Placement],
    pub(crate) size: f64,
}

impl Problem<'_> {
    fn placement(&self, x: &[f64], k: usize, root: usize) -> Placement {
        let base = self.base[root];
        let w = Vec3::new(x[6 * k], x[6 * k + 1], x[6 * k + 2]);
        let t = Vec3::new(x[6 * k + 3], x[6 * k + 4], x[6 * k + 5]) * self.size;
        let angle = w.length();
        let mut out = *base.transform();
        if angle > 0.0
            && let Ok(turn) =
                Transform3::rotation_around_axis(base.origin(), w * (1.0 / angle), angle)
        {
            out = turn.compose(&out);
        }
        Placement::from_raw(Transform3::translation(t.x, t.y, t.z).compose(&out))
    }

    fn world(&self, x: &[f64], r: &Resolved) -> Primitive {
        match self.free.iter().position(|&f| f == r.root) {
            Some(k) => r.local.placed(&self.placement(x, k, r.root)),
            None => r.local.placed(&self.base[r.root]),
        }
    }

    /// Row ranges per mate (same layout as [`Problem::residual`]).
    pub(crate) fn rows(&self, x: &[f64]) -> Result<Vec<std::ops::Range<usize>>, AssemblyError> {
        let mut out = Vec::new();
        let mut at = 0;
        for (mate, a, b, about) in self.resolved {
            let n = mate_rows(
                mate,
                &self.world(x, a),
                &self.world(x, b),
                about.as_ref().map(|r| self.world(x, r)),
                self.size,
            )?
            .len();
            out.push(at..at + n);
            at += n;
        }
        Ok(out)
    }

    pub(crate) fn residual(&self, x: &[f64]) -> Result<Vec<f64>, AssemblyError> {
        let mut out = Vec::new();
        for (mate, a, b, about) in self.resolved {
            out.extend(mate_rows(
                mate,
                &self.world(x, a),
                &self.world(x, b),
                about.as_ref().map(|r| self.world(x, r)),
                self.size,
            )?);
        }
        Ok(out)
    }

    pub(crate) fn jacobian(&self, x: &[f64]) -> Result<Vec<Vec<f64>>, AssemblyError> {
        let h = 1e-7;
        let mut columns = Vec::with_capacity(x.len());
        for j in 0..x.len() {
            let mut plus = x.to_vec();
            let mut minus = x.to_vec();
            plus[j] += h;
            minus[j] -= h;
            let (rp, rm) = (self.residual(&plus)?, self.residual(&minus)?);
            columns.push(
                rp.iter()
                    .zip(&rm)
                    .map(|(p, m)| (p - m) / (2.0 * h))
                    .collect::<Vec<f64>>(),
            );
        }
        let m = columns.first().map_or(0, Vec::len);
        Ok((0..m)
            .map(|i| columns.iter().map(|c| c[i]).collect())
            .collect())
    }
}

/// The residual rows of one mate between world primitives.
fn mate_rows(
    mate: &Mate,
    a: &Primitive,
    b: &Primitive,
    about: Option<Primitive>,
    size: f64,
) -> Result<Vec<f64>, AssemblyError> {
    let refuse = |detail: &str| AssemblyError::Mate {
        mate: mate.describe(),
        detail: detail.to_string(),
    };
    let v3 = |v: Vec3| [v.x, v.y, v.z];
    let lin = |v: Vec3| [v.x / size, v.y / size, v.z / size];
    Ok(match (&mate.kind, a, b) {
        (MateKind::Coincident { aligned }, Primitive::Plane(pa, na), Primitive::Plane(pb, nb))
        | (
            MateKind::Distance { aligned, .. },
            Primitive::Plane(pa, na),
            Primitive::Plane(pb, nb),
        ) => {
            let s = if *aligned { -1.0 } else { 1.0 };
            let offset = match mate.kind {
                MateKind::Distance { distance, .. } => distance,
                _ => 0.0,
            };
            let mut rows = v3(*na + *nb * s).to_vec();
            rows.push(((*pb - *pa).dot(na) - offset) / size);
            rows
        }
        (
            MateKind::Coincident { .. } | MateKind::Concentric,
            Primitive::Axis(pa, da),
            Primitive::Axis(pb, db),
        ) => {
            let mut rows = v3(da.cross(db)).to_vec();
            rows.extend(lin((*pb - *pa).cross(da)));
            rows
        }
        (MateKind::AxialOffset { distance }, Primitive::Axis(pa, da), Primitive::Axis(pb, _)) => {
            vec![((*pb - *pa).dot(da) - distance) / size]
        }
        (MateKind::AxialOffset { .. }, _, _) => {
            return Err(refuse("an axial offset needs two axes"));
        }
        (MateKind::Coincident { .. }, Primitive::Point(pa), Primitive::Point(pb)) => {
            lin(*pb - *pa).to_vec()
        }
        (MateKind::Coincident { .. }, Primitive::Point(p), Primitive::Plane(o, n))
        | (MateKind::Coincident { .. }, Primitive::Plane(o, n), Primitive::Point(p)) => {
            vec![(*p - *o).dot(n) / size]
        }
        (MateKind::Distance { distance, .. }, Primitive::Point(p), Primitive::Plane(o, n))
        | (MateKind::Distance { distance, .. }, Primitive::Plane(o, n), Primitive::Point(p)) => {
            vec![((*p - *o).dot(n) - distance) / size]
        }
        (
            MateKind::Coincident { .. } | MateKind::Concentric,
            Primitive::Point(p),
            Primitive::Axis(o, d),
        )
        | (
            MateKind::Coincident { .. } | MateKind::Concentric,
            Primitive::Axis(o, d),
            Primitive::Point(p),
        ) => lin((*p - *o).cross(d)).to_vec(),
        (MateKind::Coincident { .. }, Primitive::Plane(o, n), Primitive::Axis(p, d))
        | (MateKind::Coincident { .. }, Primitive::Axis(p, d), Primitive::Plane(o, n)) => {
            vec![d.dot(n), (*p - *o).dot(n) / size]
        }
        (MateKind::Angle { degrees, .. }, _, _) => {
            let (Some(ua), Some(ub)) = (a.direction(), b.direction()) else {
                return Err(refuse(
                    "an angle needs two directions (plane normals or axes)",
                ));
            };
            let theta = degrees.to_radians();
            match about.and_then(|k| k.direction()) {
                Some(k) => {
                    let signed = ua.cross(&ub).dot(&k).patan2(ua.dot(&ub));
                    let mut d = signed - theta;
                    d = (d + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                        - std::f64::consts::PI;
                    vec![d]
                }
                None => vec![ua.dot(&ub) - theta.pcos()],
            }
        }
        (MateKind::Concentric, _, _) => {
            return Err(refuse(
                "concentric needs two axes (cylinder, cone or torus faces, or datum axes) or a point and an axis",
            ));
        }
        _ => return Err(refuse("these two entities cannot take this mate")),
    })
}

pub(crate) fn max_abs(r: &[f64]) -> f64 {
    r.iter().map(|v| v.abs()).fold(0.0, f64::max)
}

fn point_of(p: &Primitive) -> Point3 {
    match *p {
        Primitive::Plane(o, _) | Primitive::Axis(o, _) | Primitive::Point(o) => o,
    }
}

trait CoordsMax {
    fn coords_max(&self) -> f64;
}

impl CoordsMax for Point3 {
    fn coords_max(&self) -> f64 {
        self.x.abs().max(self.y.abs()).max(self.z.abs())
    }
}
