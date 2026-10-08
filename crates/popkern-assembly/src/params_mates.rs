//! Assembly parameters that drive mates (PK-701): a hinge's open angle
//! as a slider. A mate bound to a parameter takes
//! `offset + scale × value` as its angle (degrees) or distance; setting
//! the parameter re-solves the mates and moves occurrences — no part
//! replays, no history step (a slider is not authoring).

use crate::assembly::Assembly;
use crate::error::AssemblyError;
use crate::mates::{MateDrive, MateKind};
use crate::params::param_error;

/// The widest step a driven angle takes in one solve (degrees); larger
/// moves are continued through intermediate solutions (§5.3).
pub(crate) const MAX_ANGLE_STEP: f64 = 30.0;

impl Assembly {
    /// Bind mate `index` (an angle or distance mate) to assembly
    /// parameter `param` (`value = offset + scale × param`) and solve.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Param`] for an unknown parameter or mate, a mate
    /// with no value to drive, or a mate already driven; the solve's
    /// errors (the binding is then undone).
    pub fn drive_mate(
        &mut self,
        param: &str,
        index: usize,
        scale: f64,
        offset: f64,
    ) -> Result<(), AssemblyError> {
        if !self.params.iter().any(|p| p.name == param) {
            return Err(param_error(param, "no assembly parameter with this name"));
        }
        let Some(mate) = self.mates.get(index) else {
            return Err(param_error(param, &format!("there is no mate {index}")));
        };
        if !matches!(
            mate.kind,
            MateKind::Angle { .. } | MateKind::Distance { .. }
        ) {
            return Err(param_error(
                param,
                "only angle and distance mates have a value a parameter can drive",
            ));
        }
        if mate.drive.is_some() {
            return Err(param_error(
                param,
                &format!("mate {index} is already driven"),
            ));
        }
        if !(scale.is_finite() && offset.is_finite()) {
            return Err(param_error(param, "scale and offset must be finite"));
        }
        let kind = mate.kind.clone();
        self.mates[index].drive = Some(MateDrive {
            param: param.to_string(),
            scale,
            offset,
        });
        if let Err(error) = self.apply_mate_drives(param) {
            self.mates[index].drive = None;
            self.mates[index].kind = kind;
            return Err(error);
        }
        Ok(())
    }

    /// Limit assembly parameter `name` to `[min, max]` (display unit).
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Param`] for an unknown parameter, an empty range,
    /// or a current value outside it.
    pub fn set_param_bounds(
        &mut self,
        name: &str,
        min: f64,
        max: f64,
    ) -> Result<(), AssemblyError> {
        let Some(param) = self.params.iter_mut().find(|p| p.name == name) else {
            return Err(param_error(name, "no assembly parameter with this name"));
        };
        if !(min.is_finite() && max.is_finite() && min <= max) {
            return Err(param_error(name, &format!("[{min}, {max}] is not a range")));
        }
        let current = param.display_value();
        if !(min..=max).contains(&current) {
            return Err(param_error(
                name,
                &format!("its value {current} is outside [{min}, {max}]"),
            ));
        }
        param.bounds = Some([min, max]);
        Ok(())
    }

    /// Whether any mate is driven by `param`.
    pub fn drives_mates(&self, param: &str) -> bool {
        self.mates
            .iter()
            .any(|m| m.drive.as_ref().is_some_and(|d| d.param == param))
    }

    /// Set every mate `param` drives from its current value and re-solve
    /// (transactional: a failed solve restores the mate values and the
    /// placements).
    pub(crate) fn apply_mate_drives(&mut self, param: &str) -> Result<(), AssemblyError> {
        if !self.drives_mates(param) {
            return Ok(());
        }
        let Some(value) = self
            .params
            .iter()
            .find(|p| p.name == param)
            .map(|p| p.display_value())
        else {
            return Err(param_error(param, "no assembly parameter with this name"));
        };
        let mates = self.mates.clone();
        let placements: Vec<_> = self.root.iter().map(|i| i.placement).collect();
        // Continuation: an angle residual is an atan2 and wraps, so a
        // large jump solved in one go can land on another branch (or
        // not converge). Step every driven angle by at most
        // `MAX_ANGLE_STEP` from the solution before.
        let targets: Vec<Option<(f64, f64)>> = self
            .mates
            .iter()
            .map(|mate| {
                let drive = mate.drive.as_ref().filter(|d| d.param == param)?;
                let from = match mate.kind {
                    MateKind::Angle { degrees, .. } => degrees,
                    MateKind::Distance { distance, .. } => distance,
                    _ => return None,
                };
                Some((from, drive.value(value)))
            })
            .collect();
        let widest = self
            .mates
            .iter()
            .zip(&targets)
            .filter_map(|(mate, t)| match (mate.kind.clone(), t) {
                (MateKind::Angle { .. }, Some((from, to))) => Some((to - from).abs()),
                _ => None,
            })
            .fold(0.0_f64, f64::max);
        let steps = ((widest / MAX_ANGLE_STEP).ceil() as usize).max(1);
        for k in 1..=steps {
            let f = k as f64 / steps as f64;
            for (mate, t) in self.mates.iter_mut().zip(&targets) {
                let Some((from, to)) = *t else { continue };
                let v = from + (to - from) * f;
                match &mut mate.kind {
                    MateKind::Angle { degrees, .. } => *degrees = v,
                    MateKind::Distance { distance, .. } => *distance = v,
                    _ => {}
                }
            }
            if let Err(error) = self.solve_in_place(false) {
                self.mates = mates;
                for (instance, placement) in self.root.iter_mut().zip(placements) {
                    instance.placement = placement;
                }
                return Err(error);
            }
        }
        Ok(())
    }
}
