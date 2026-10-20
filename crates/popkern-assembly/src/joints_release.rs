//! Released joints and their limits (PK-808, six-axis-arm §11 E4).
//!
//! A released joint's angle mate leaves the solve: the joint turns
//! freely, its angle decided by the other mates (a flange mated to a
//! target), and its parameter reads back the solved angle. Driven again,
//! the angle mate rejoins at the angle the joint has, so nothing jumps.
//! After every solve each released joint's angle — measured through
//! its zero directions about its axis, as the angle mate measures it —
//! must lie within its parameter's limits; one outside refuses the
//! solve transactionally ([`AssemblyError::JointLimit`], nothing moves).

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::joints::JointDrive;
use crate::mates::{MateKind, Resolved};
use crate::params::param_error;
use crate::placement::Placement;

impl Assembly {
    /// Release joint `name` (it must drive itself): its angle is left to
    /// the other mates and read back into its parameter.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Param`] for an unknown joint or one following
    /// another parameter.
    pub fn release_joint(&mut self, name: &str) -> Result<(), AssemblyError> {
        let joint = self.own_joint(name)?;
        self.released.insert(joint.mates[2]);
        Ok(())
    }

    /// Drive joint `name` by its parameter again, from the angle it has.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Param`] for an unknown joint or one following
    /// another parameter; the solve's errors (the joint stays released).
    pub fn drive_joint(&mut self, name: &str) -> Result<(), AssemblyError> {
        let joint = self.own_joint(name)?;
        let index = joint.mates[2];
        if !self.released.contains(&index) {
            return Ok(());
        }
        let angle = self
            .params
            .iter()
            .find(|p| p.name == name)
            .map(|p| p.display_value())
            .ok_or_else(|| param_error(name, "the joint's parameter is gone"))?;
        let kind = self.mates[index].kind.clone();
        if let MateKind::Angle { degrees, .. } = &mut self.mates[index].kind {
            *degrees = angle;
        }
        self.released.remove(&index);
        if let Err(error) = self.apply_mate_drives(name) {
            self.mates[index].kind = kind;
            self.released.insert(index);
            return Err(error);
        }
        Ok(())
    }

    /// The released joints, in joint order.
    pub fn released_joints(&self) -> Vec<String> {
        self.joints
            .iter()
            .filter(|j| self.released.contains(&j.mates[2]))
            .map(|j| j.name.clone())
            .collect()
    }

    /// Whether `param` is a released joint's (it reads back the solve
    /// and cannot be set).
    pub(crate) fn is_released_param(&self, param: &str) -> bool {
        self.joints
            .iter()
            .any(|j| j.name == param && self.released.contains(&j.mates[2]))
    }

    fn own_joint(&self, name: &str) -> Result<crate::joints::Joint, AssemblyError> {
        let joint = self
            .joints
            .iter()
            .find(|j| j.name == name)
            .ok_or_else(|| param_error(name, "no joint with this name"))?;
        if joint.drive != JointDrive::Param {
            return Err(param_error(
                name,
                "the joint follows another parameter; release its leader",
            ));
        }
        Ok(joint.clone())
    }

    /// After a solve: every released joint's angle read back into its
    /// parameter, or, when one lies outside its limits, the placements
    /// restored to `before` and the solve refused. Returns whether a
    /// mate the read-back angles drive (a follower: a reducer's rotor,
    /// discs and cage) changed and must be solved again.
    pub(crate) fn settle_released(
        &mut self,
        resolved: &[(crate::mates::Mate, Resolved, Resolved, Option<Resolved>)],
        before: &[Placement],
    ) -> Result<bool, AssemblyError> {
        if self.released.is_empty() {
            return Ok(false);
        }
        let mut angles = Vec::new();
        for joint in &self.joints {
            let index = joint.mates[2];
            if !self.released.contains(&index) {
                continue;
            }
            let Some((_, za, zb, Some(axis))) = resolved.get(index) else {
                continue;
            };
            let world = |r: &Resolved| r.local.placed(&self.root[r.root].placement).direction();
            let (Some(k), Some(ua), Some(ub)) = (world(axis), world(za), world(zb)) else {
                continue;
            };
            let angle = crate::joints::signed_degrees(ua, ub, k);
            let limits = self
                .params
                .iter()
                .find(|p| p.name == joint.name)
                .and_then(|p| p.bounds);
            if let Some([lo, hi]) = limits
                && !(lo..=hi).contains(&angle)
            {
                for (instance, placement) in self.root.iter_mut().zip(before) {
                    instance.placement = *placement;
                }
                return Err(AssemblyError::JointLimit {
                    joint: joint.name.clone(),
                    angle,
                    limits: [lo, hi],
                });
            }
            angles.push((joint.name.clone(), index, angle));
        }
        let mut followers_moved = false;
        for (name, index, angle) in angles {
            if let Some(param) = self.params.iter_mut().find(|p| p.name == name) {
                param.value = param.unit.to_si(angle);
            }
            if let MateKind::Angle { degrees, .. } = &mut self.mates[index].kind {
                *degrees = angle;
            }
            for (i, mate) in self.mates.iter_mut().enumerate() {
                let Some(drive) = mate.drive.as_ref().filter(|d| d.param == name) else {
                    continue;
                };
                let target = drive.value(angle);
                let value = match &mut mate.kind {
                    MateKind::Angle { degrees, .. } if i != index => degrees,
                    MateKind::Distance { distance, .. } => distance,
                    _ => continue,
                };
                // Degrees and model units: a micro-unit change is no move.
                if (*value - target).abs() > 1e-9 {
                    *value = target;
                    followers_moved = true;
                }
            }
        }
        Ok(followers_moved)
    }
}
