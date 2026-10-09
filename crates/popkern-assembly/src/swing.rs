//! Swinging one assembly parameter: where does the motion first
//! collide (PK-701; arm scale, six-axis-arm §5.4).
//!
//! A swing samples the parameter from `from` to `to`, checks
//! interference at each value, and bisects the first contact. It
//! reports what the MOTION does, so it checks only the pairs the motion
//! changes:
//!
//! - **Relative-pose restriction.** Before sampling, the swing probes
//!   two other values and groups the root instances that keep their
//!   relative placement (a rigid cluster: the arm's links upstream of
//!   the swung joint, or everything downstream of it). A pair inside
//!   one cluster has the same status at every value, so it is not
//!   checked: a fastener in its tap-drill bore, a pin in its seat.
//!   Run `interference` for those.
//! - **Follower exclusion.** The joints that follow the swung parameter
//!   (a reducer's rotor, discs, roller cage) turn inside their housing
//!   at their own ratios; their fits are a drive train's business,
//!   checked by its own sweep at the resolution the fits need. Their
//!   roots are left out and named in the report.
//! - **Contact.** An overlap, or a touch between occurrences that were
//!   apart at `from` (a constant contact is a fit, not a collision).
//!
//! `clearance` bounds the pairs measured for the report's closest
//! approach; unbounded (the default) measures every moving pair.
//! `exclude` names pairs the caller vouches for (a joint's own races).

use serde::Serialize;

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::interference::{InterferenceOptions, PairStatus};
use crate::joints::JointDrive;
use crate::params::param_error;
use crate::placement::Placement;

/// A swing's range and resolution.
#[derive(Debug, Clone)]
pub struct SwingOptions {
    /// The first value.
    pub from: f64,
    /// The last value.
    pub to: f64,
    /// The sampling step (positive).
    pub step: f64,
    /// The first contact is bisected to this (positive).
    pub resolution: f64,
    /// Pairs farther apart than this are not measured; `None`
    /// measures every moving pair.
    pub clearance: Option<f64>,
    /// Occurrence pairs left out, as path prefixes `[a, b]` (either
    /// order): fits the swing should not re-measure at every step, such
    /// as a joint's bearing races and its flange in its own housing,
    /// which turn about their common axis.
    pub exclude: Vec<[String; 2]>,
}

/// One value of a swing.
#[derive(Debug, Clone, Serialize)]
pub struct SwingSample {
    /// The parameter value.
    pub value: f64,
    /// The smallest distance between moving occurrences apart there.
    pub min_distance: f64,
    /// Whether occurrences collide there.
    pub contact: bool,
}

/// What [`Assembly::swing`] found.
#[derive(Debug, Clone, Serialize)]
pub struct SwingReport {
    /// The parameter swept.
    pub param: String,
    /// The first value where occurrences collide (refined), if any.
    pub first_contact: Option<f64>,
    /// The occurrences (paths) that collide there.
    pub first_contact_pair: Option<(String, String)>,
    /// The smallest clearance seen among collision-free samples.
    pub min_clearance: Option<f64>,
    /// Where that clearance was seen.
    pub min_clearance_at: Option<f64>,
    /// Every sample taken on the sweep (bisection probes excluded).
    pub samples: Vec<SwingSample>,
    /// Root instances grouped by the rigid clusters the motion keeps.
    pub clusters: Vec<Vec<String>>,
    /// Follower roots left out of the check (their drive train's).
    pub followers: Vec<String>,
}

/// Which pairs a swing checks: roots in different rigid clusters,
/// neither a follower of the swung parameter.
struct Scope {
    cluster: std::collections::HashMap<String, usize>,
    followers: std::collections::HashSet<String>,
    exclude: Vec<[String; 2]>,
}

impl Scope {
    fn checks(&self, a: &str, b: &str) -> bool {
        if self.exclude.iter().any(|[x, y]| {
            (a.starts_with(x.as_str()) && b.starts_with(y.as_str()))
                || (a.starts_with(y.as_str()) && b.starts_with(x.as_str()))
        }) {
            return false;
        }
        let root = |p: &str| p.split('/').next().unwrap_or(p).to_string();
        let (ra, rb) = (root(a), root(b));
        if self.followers.contains(&ra) || self.followers.contains(&rb) {
            return false;
        }
        match (self.cluster.get(&ra), self.cluster.get(&rb)) {
            (Some(x), Some(y)) => x != y,
            _ => true,
        }
    }
}

struct Sample {
    min_distance: f64,
    contact: Option<(String, String)>,
    touching: Vec<(String, String)>,
}

/// Probe values (fractions of the range) that decide the clusters: two
/// generic values, so a full turn (whose ends coincide) and a single
/// coincidental return are both caught.
const PROBES: [f64; 2] = [0.37, 0.81];

/// Placements equal to this, translation relative to the assembly's
/// scale, are one rigid relation.
const RIGID: f64 = 1e-9;

/// The swing's base mesh chord, as a fraction of each part's size: ten
/// times the interference default. A swing measures the same pairs at
/// every step; near pairs are certified on finer meshes anyway
/// (`interference_certify`), so the base mesh only has to sort the far
/// pairs from the near — and the default chord on a threaded screw
/// (0.017 mm on a 17 mm part) made its distance queries the swing's
/// whole cost.
const SWING_CHORD_FRACTION: f64 = 1e-2;

impl Assembly {
    /// Sweep assembly parameter `param` from `from` to `to` in steps of
    /// `step`, refining the first contact by bisection to `resolution`
    /// (see the module docs for which pairs it checks). The parameter
    /// is restored afterwards.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Param`] for an unknown parameter or a bad range;
    /// the interference check's and the solve's errors (the parameter is
    /// restored first).
    pub fn swing(
        &mut self,
        param: &str,
        from: f64,
        to: f64,
        step: f64,
        resolution: f64,
    ) -> Result<SwingReport, AssemblyError> {
        self.swing_with(
            param,
            &SwingOptions {
                from,
                to,
                step,
                resolution,
                clearance: None,
                exclude: Vec::new(),
            },
        )
    }

    /// [`Assembly::swing`] with every option.
    ///
    /// # Errors
    ///
    /// As [`Assembly::swing`]; a negative clearance is a bad range.
    pub fn swing_with(
        &mut self,
        param: &str,
        options: &SwingOptions,
    ) -> Result<SwingReport, AssemblyError> {
        let Some(original) = self
            .params
            .iter()
            .find(|p| p.name == param)
            .map(|p| p.display_value())
        else {
            return Err(param_error(param, "no assembly parameter with this name"));
        };
        let o = options;
        if !(o.from.is_finite()
            && o.to.is_finite()
            && o.step > 0.0
            && o.resolution > 0.0
            && o.clearance.is_none_or(|c| c >= 0.0))
        {
            return Err(param_error(
                param,
                "the swing needs finite ends, positive steps and a non-negative clearance",
            ));
        }
        let result = self.swing_inner(param, o);
        let restore = self.set_assembly_param(param, original);
        let report = result?;
        restore?;
        Ok(report)
    }

    fn swing_inner(&mut self, param: &str, o: &SwingOptions) -> Result<SwingReport, AssemblyError> {
        let scope = self.swing_scope(param, o)?;
        let first = self.swing_sample(param, o.from, None, &scope, o.clearance)?;
        let baseline: Vec<(String, String)> = first.touching.clone();
        let mut samples = vec![SwingSample {
            value: o.from,
            min_distance: first.min_distance,
            contact: first.contact.is_some(),
        }];
        let mut first_contact = first.contact.map(|pair| (o.from, pair));
        let dir = if o.to >= o.from { 1.0 } else { -1.0 };
        let count = ((o.to - o.from).abs() / o.step).ceil() as usize;
        let mut last_clear = o.from;
        for k in 1..=count {
            if first_contact.is_some() {
                break;
            }
            let v = if k == count {
                o.to
            } else {
                o.from + dir * o.step * k as f64
            };
            let s = self.swing_sample(param, v, Some(&baseline), &scope, o.clearance)?;
            samples.push(SwingSample {
                value: v,
                min_distance: s.min_distance,
                contact: s.contact.is_some(),
            });
            match s.contact {
                Some(pair) => {
                    // Bisect between the last clear value and this one.
                    let (mut clear, mut hit, mut at) = (last_clear, v, pair);
                    while (hit - clear).abs() > o.resolution {
                        let mid = 0.5 * (clear + hit);
                        let probe =
                            self.swing_sample(param, mid, Some(&baseline), &scope, o.clearance)?;
                        match probe.contact {
                            Some(pair) => {
                                hit = mid;
                                at = pair;
                            }
                            None => clear = mid,
                        }
                    }
                    first_contact = Some((hit, at));
                }
                None => last_clear = v,
            }
        }
        let closest = samples
            .iter()
            .filter(|s| !s.contact && s.min_distance.is_finite())
            .min_by(|a, b| a.min_distance.total_cmp(&b.min_distance));
        let mut clusters: Vec<Vec<String>> = Vec::new();
        let mut members: Vec<(&String, &usize)> = scope.cluster.iter().collect();
        members.sort();
        for (root, &c) in members {
            if clusters.len() <= c {
                clusters.resize(c + 1, Vec::new());
            }
            clusters[c].push(root.clone());
        }
        let mut followers: Vec<String> = scope.followers.into_iter().collect();
        followers.sort();
        Ok(SwingReport {
            param: param.to_string(),
            first_contact: first_contact.as_ref().map(|c| c.0),
            first_contact_pair: first_contact.map(|c| c.1),
            min_clearance: closest.map(|s| s.min_distance),
            min_clearance_at: closest.map(|s| s.value),
            samples,
            clusters,
            followers,
        })
    }

    /// The pairs a swing of `param` checks: root clusters from the
    /// placements at `from` and two probe values, and the followers.
    fn swing_scope(&mut self, param: &str, o: &SwingOptions) -> Result<Scope, AssemblyError> {
        let followers: std::collections::HashSet<String> = self
            .joints
            .iter()
            .filter(|j| matches!(&j.drive, JointDrive::Follow { leader, .. } if leader == param))
            .filter_map(|j| self.mates.get(j.mates[0]))
            .map(|m| m.b.instance.split('/').next().unwrap_or("").to_string())
            .collect();
        let roots = |a: &Assembly| -> Vec<(String, Placement)> {
            a.root
                .iter()
                .map(|i| (i.name.clone(), i.placement))
                .collect()
        };
        self.set_assembly_param(param, o.from)?;
        let at_from = roots(self);
        let scale = at_from
            .iter()
            .map(|(_, p)| {
                let o = p.origin();
                (o.x * o.x + o.y * o.y + o.z * o.z).sqrt()
            })
            .fold(1.0, f64::max);
        let mut deltas: Vec<Vec<Placement>> = vec![Vec::new(); at_from.len()];
        for f in PROBES {
            self.set_assembly_param(param, o.from + f * (o.to - o.from))?;
            for (k, (_, p)) in roots(self).iter().enumerate() {
                if let Some((_, p0)) = at_from.get(k) {
                    deltas[k].push(p0.inverse().then(p));
                }
            }
        }
        let same = |a: &[Placement], b: &[Placement]| {
            a.iter().zip(b).all(|(x, y)| {
                let (mx, my) = (&x.transform().matrix, &y.transform().matrix);
                (0..3).all(|i| {
                    (0..4).all(|j| {
                        let s = if j == 3 { scale } else { 1.0 };
                        (mx[i][j] - my[i][j]).abs() <= RIGID * s
                    })
                })
            })
        };
        let mut cluster = std::collections::HashMap::new();
        let mut representatives: Vec<usize> = Vec::new();
        for k in 0..at_from.len() {
            let found = representatives
                .iter()
                .position(|&r| same(&deltas[r], &deltas[k]));
            let c = found.unwrap_or_else(|| {
                representatives.push(k);
                representatives.len() - 1
            });
            cluster.insert(at_from[k].0.clone(), c);
        }
        Ok(Scope {
            cluster,
            followers,
            exclude: o.exclude.clone(),
        })
    }

    fn swing_sample(
        &mut self,
        param: &str,
        value: f64,
        baseline: Option<&[(String, String)]>,
        scope: &Scope,
        clearance: Option<f64>,
    ) -> Result<Sample, AssemblyError> {
        self.set_assembly_param(param, value)?;
        // Unbounded: every moving pair the bounds bring within reach
        // (the clearance the size of the assembly; the BVH prunes).
        let report = self.interference_where(
            &InterferenceOptions {
                clearance: clearance.unwrap_or(f64::MAX / 4.0),
                chord_fraction: SWING_CHORD_FRACTION,
                overlap_volume: false,
                ..InterferenceOptions::default()
            },
            &|a, b| scope.checks(a, b),
        )?;
        let mut min_distance = f64::INFINITY;
        let mut contact = None;
        let mut touching = Vec::new();
        for pair in &report.pairs {
            let key = (pair.a.clone(), pair.b.clone());
            match &pair.status {
                PairStatus::Clear { distance, .. } => min_distance = min_distance.min(*distance),
                PairStatus::Touching { .. } => {
                    if contact.is_none() && baseline.is_some_and(|b| !b.contains(&key)) {
                        contact = Some(key.clone());
                    }
                    touching.push(key);
                }
                PairStatus::Interfering { .. } => {
                    contact.get_or_insert(key);
                }
            }
        }
        Ok(Sample {
            min_distance,
            contact,
            touching,
        })
    }
}
