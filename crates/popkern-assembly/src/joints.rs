//! Revolute joints (design `docs/architecture/six-axis-arm.md` §4.3):
//! one named object for "this occurrence turns about that axis", made
//! of three mates the solver already understands —
//!
//! - **concentric** between the two axes;
//! - **axial offset** between the axes' origins (the axial seat);
//! - **angle** between two zero directions about the first axis,
//!   driven either by a new assembly parameter (the joint's slider) or,
//!   for a *follower*, by another parameter times a ratio (a reducer's
//!   rotor turning `−(N−1)×` its joint, a bearing cage at `½×`).
//!
//! The joint is recorded beside its mates so show-as-code emits one
//! verb and reports can name it. A follower creates no parameter, so a
//! follower can never be followed: ratio chains go through parameters
//! and cannot cycle by construction.

use popkern_core::geometry::Vec3;
use popkern_history::unit::Unit;
use popkern_math::PortableMath;
use serde::{Deserialize, Serialize};

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::mate_solve::MateReport;
use crate::mates::{Mate, MateEntity, MateKind, MateRef, Primitive};

/// A zero direction must be perpendicular to its joint axis to this
/// cosine: the angle mate measures full directions, so a tilted zero
/// would make the slider angle a non-linear function of the turn.
const PERPENDICULAR_COSINE: f64 = 1e-6;

/// What drives a joint's angle.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum JointDrive {
    /// The joint's own assembly parameter (degrees), named after it.
    Param,
    /// Another assembly parameter: `angle = phase + ratio × leader`.
    Follow {
        /// The leading parameter.
        leader: String,
        /// Turns of this joint per turn of the leader.
        ratio: f64,
        /// The angle at leader zero, degrees.
        phase: f64,
    },
}

/// A recorded joint.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Joint {
    /// Its name (and, when it drives itself, its parameter's name).
    pub name: String,
    /// Indices of its mates: concentric, axial offset, angle.
    pub mates: [usize; 3],
    /// What drives the angle.
    pub drive: JointDrive,
}

/// How [`Assembly::add_revolute`] sets up a joint.
#[derive(Debug, Clone, Default)]
pub struct RevoluteOptions {
    /// Directions whose signed angle about the first axis is the joint
    /// angle (plane normals or axes, perpendicular to the joint axis).
    /// `None`: a direction perpendicular to each axis in its own
    /// occurrence frame.
    pub zero: Option<(MateRef, MateRef)>,
    /// The starting angle (degrees); `None` keeps the current one.
    pub angle: Option<f64>,
    /// Slider limits (degrees).
    pub limits: Option<[f64; 2]>,
    /// The second axis origin's offset along the first axis.
    pub offset: f64,
    /// Follow another parameter at a ratio instead of creating one.
    pub follow: Option<(String, f64)>,
    /// A follower's angle at leader zero; `None` keeps the current pose.
    pub phase: Option<f64>,
}

/// What [`Assembly::add_revolute`] made.
#[derive(Debug, Clone, Serialize)]
pub struct RevoluteReport {
    /// The joint's name.
    pub joint: String,
    /// Its angle after the solve (degrees).
    pub angle: f64,
    /// The parameter it created, when it drives itself.
    pub param: Option<String>,
    /// The solve after the last mate.
    pub solve: MateReport,
}

fn joint_entity(joint: &str, role: &'static str, found: &str) -> AssemblyError {
    AssemblyError::JointEntity {
        joint: joint.to_string(),
        role,
        found: found.to_string(),
    }
}

fn kind_name(p: &Primitive) -> &'static str {
    match p {
        Primitive::Plane(..) => "a plane",
        Primitive::Axis(..) => "an axis",
        Primitive::Point(..) => "a point",
    }
}

/// A unit direction perpendicular to `d`, deterministic.
fn perpendicular(d: Vec3) -> Vec3 {
    let x = Vec3::new(1.0, 0.0, 0.0);
    let y = Vec3::new(0.0, 1.0, 0.0);
    let c = if d.cross(&x).length() > 0.5 {
        d.cross(&x)
    } else {
        d.cross(&y)
    };
    c.normalized().unwrap_or(y)
}

/// The signed angle (degrees) from `ua` to `ub` about `k`.
fn signed_degrees(ua: Vec3, ub: Vec3, k: Vec3) -> f64 {
    ua.cross(&ub).dot(&k).patan2(ua.dot(&ub)).to_degrees()
}

impl Assembly {
    /// The joints, in the order they were added.
    pub fn joints(&self) -> &[Joint] {
        &self.joints
    }

    /// The joint whose mates include mate `index`.
    pub fn joint_of_mate(&self, index: usize) -> Option<&Joint> {
        self.joints.iter().find(|j| j.mates.contains(&index))
    }

    /// Add a revolute joint `name`: `axis_b`'s occurrence turns about
    /// `axis_a` (see the module docs). Without `follow`, a new assembly
    /// parameter `name` (degrees) drives the angle; it starts at
    /// `angle`, or at the current angle so nothing moves.
    ///
    /// Transactional: on any error every mate it added is removed and
    /// every placement restored.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::JointEntity`] for an axis that is not an axis or
    /// a zero entity with no direction perpendicular to its axis;
    /// [`AssemblyError::Param`] for a name clash, an unknown leader, a
    /// zero or non-finite ratio, or limits that exclude the angle; the
    /// solve's errors.
    pub fn add_revolute(
        &mut self,
        name: &str,
        axis_a: MateRef,
        axis_b: MateRef,
        options: &RevoluteOptions,
    ) -> Result<RevoluteReport, AssemblyError> {
        let param_error = crate::params::param_error;
        if self.joints.iter().any(|j| j.name == name) {
            return Err(param_error(name, "a joint with this name already exists"));
        }
        if options.follow.is_none() && self.params.iter().any(|p| p.name == name) {
            return Err(param_error(
                name,
                "an assembly parameter with this name already exists (a joint names its slider)",
            ));
        }
        let leader = match &options.follow {
            Some((leader, ratio)) => {
                if !ratio.is_finite() || *ratio == 0.0 {
                    return Err(param_error(
                        name,
                        &format!(
                            "a follower's ratio must be finite and non-zero, got {ratio} \
                             (a joint that never turns is a fixed mate)"
                        ),
                    ));
                }
                let Some(p) = self.params.iter().find(|p| &p.name == leader) else {
                    return Err(param_error(
                        name,
                        &format!("it follows '{leader}', which is no assembly parameter"),
                    ));
                };
                Some((leader.clone(), *ratio, p.display_value()))
            }
            None => None,
        };
        let ra = self.resolve_mate_ref(&axis_a)?;
        let rb = self.resolve_mate_ref(&axis_b)?;
        let Primitive::Axis(_, da_local) = ra.local else {
            return Err(joint_entity(name, "axis", kind_name(&ra.local)));
        };
        let Primitive::Axis(_, db_local) = rb.local else {
            return Err(joint_entity(name, "axis", kind_name(&rb.local)));
        };
        // Zero directions: given, or perpendicular to each axis in its
        // root instance's frame (a literal there resolves unchanged).
        let (zero_a, zero_b) = match &options.zero {
            Some((za, zb)) => (za.clone(), zb.clone()),
            None => {
                let literal = |root: usize, d: Vec3| -> MateRef {
                    let p = perpendicular(d);
                    MateRef {
                        instance: self.root[root].name.clone(),
                        entity: MateEntity::Axis {
                            origin: [0.0, 0.0, 0.0],
                            direction: [p.x, p.y, p.z],
                        },
                    }
                };
                (literal(ra.root, da_local), literal(rb.root, db_local))
            }
        };
        let za = self.resolve_mate_ref(&zero_a)?;
        let zb = self.resolve_mate_ref(&zero_b)?;
        for (z, axis) in [(&za, &ra), (&zb, &rb)] {
            let Some(zd) = z.local.direction() else {
                return Err(joint_entity(name, "zero", kind_name(&z.local)));
            };
            // Perpendicularity is a property of one rigid body: check it
            // where the zero rides the same root as its axis.
            if let (true, Some(ad)) = (z.root == axis.root, axis.local.direction())
                && zd.dot(&ad).abs() > PERPENDICULAR_COSINE
            {
                return Err(joint_entity(
                    name,
                    "zero",
                    "a direction not perpendicular to the joint axis",
                ));
            }
        }
        // The current angle, world frame.
        let world = |r: &crate::mates::Resolved| r.local.placed(&self.root[r.root].placement);
        let (Some(ka), Some(ua), Some(ub)) = (
            world(&ra).direction(),
            world(&za).direction(),
            world(&zb).direction(),
        ) else {
            return Err(joint_entity(name, "zero", "an entity with no direction"));
        };
        let current = signed_degrees(ua, ub, ka);
        let (value, drive, created) = match leader {
            Some((leader, ratio, lead_value)) => {
                let phase = options.phase.unwrap_or(current - ratio * lead_value);
                (
                    phase + ratio * lead_value,
                    JointDrive::Follow {
                        leader,
                        ratio,
                        phase,
                    },
                    None,
                )
            }
            None => (
                options.angle.unwrap_or(current),
                JointDrive::Param,
                Some(name),
            ),
        };
        let mates_before = self.mates.len();
        let params_before = self.params.clone();
        let history_before = self.history.clone();
        let placements: Vec<_> = self.root.iter().map(|i| i.placement).collect();
        let restore = |a: &mut Assembly| {
            a.mates.truncate(mates_before);
            a.params = params_before.clone();
            a.history = history_before.clone();
            for (instance, p) in a.root.iter_mut().zip(&placements) {
                instance.placement = *p;
            }
        };
        let result = (|| -> Result<MateReport, AssemblyError> {
            self.add_mate(Mate {
                kind: MateKind::Concentric,
                a: axis_a.clone(),
                b: axis_b.clone(),
                drive: None,
            })?;
            self.add_mate(Mate {
                kind: MateKind::AxialOffset {
                    distance: options.offset,
                },
                a: axis_a.clone(),
                b: axis_b.clone(),
                drive: None,
            })?;
            if let Some(param) = created {
                self.add_param(param, value, Unit::Degrees)?;
                if let Some([lo, hi]) = options.limits {
                    self.set_param_bounds(param, lo, hi)?;
                }
            }
            let report = self.add_mate(Mate {
                kind: MateKind::Angle {
                    degrees: value,
                    about: Some(axis_a.clone()),
                },
                a: zero_a.clone(),
                b: zero_b.clone(),
                drive: None,
            })?;
            let index = self.mates.len() - 1;
            let (param, scale, offset) = match &drive {
                JointDrive::Param => (name.to_string(), 1.0, 0.0),
                JointDrive::Follow {
                    leader,
                    ratio,
                    phase,
                } => (leader.clone(), *ratio, *phase),
            };
            self.drive_mate(&param, index, scale, offset)?;
            Ok(report)
        })();
        match result {
            Ok(solve) => {
                self.joints.push(Joint {
                    name: name.to_string(),
                    mates: [mates_before, mates_before + 1, mates_before + 2],
                    drive,
                });
                Ok(RevoluteReport {
                    joint: name.to_string(),
                    angle: value,
                    param: created.map(str::to_string),
                    solve,
                })
            }
            Err(error) => {
                restore(self);
                Err(error)
            }
        }
    }
}
