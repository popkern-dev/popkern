//! A hinge as an assembly (PK-701, design `docs/architecture/hinges.md`
//! §4.4): its leaves (and pin) become part definitions placed in the
//! hinge frame and mated about the pin axis, with an assembly parameter
//! `<name>_angle` driving the open angle — a slider that moves leaf B
//! without rebuilding any part — and a swing sweep that finds where the
//! motion first collides.

use popkern_history::feature_hinge::{HingeComponentDesc, HingeKindDesc};
use popkern_history::unit::Unit;
use popkern_script::Model;
use popkern_script::body_hinge::HingeOptions;
use popkern_script::scalar::Scalar;
use serde::Serialize;

use crate::assembly::{Assembly, DefinitionId, InstanceRef, Parent};
use crate::error::AssemblyError;
use crate::interference::{InterferenceOptions, PairStatus};
use crate::mates::{Datum, Mate, MateEntity, MateKind, MateRef};
use crate::params::param_error;
use crate::placement::Placement;

/// What [`Assembly::add_hinge`] made.
#[derive(Debug, Clone)]
pub struct HingeParts {
    /// Each component's part definition, occurrence name and occurrence.
    pub parts: Vec<(HingeComponentDesc, DefinitionId, String, InstanceRef)>,
    /// The open-angle assembly parameter.
    pub angle: String,
}

fn datum(instance: &str, datum: Datum) -> MateRef {
    MateRef {
        instance: instance.to_string(),
        entity: MateEntity::Datum { datum },
    }
}

impl Assembly {
    /// Add a hinge called `name`: a part per component
    /// (`<name>_leaf_a`, `<name>_leaf_b`, `<name>_pin`) placed in the
    /// hinge frame; leaf B and the pin mated to leaf A (concentric pin
    /// axes, coincident bases), leaf B turned by an angle mate driven by
    /// the new parameter `<name>_angle` (degrees, limited to [0, 360];
    /// 180 flat). `options.angle` (a literal) sets its starting value;
    /// the leaves themselves are built flat.
    ///
    /// # Errors
    ///
    /// The hinge's typed declines (as [`AssemblyError::Script`]), name
    /// clashes, and the solve's errors.
    pub fn add_hinge(
        &mut self,
        name: &str,
        kind: HingeKindDesc,
        length: Scalar,
        options: &HingeOptions,
    ) -> Result<HingeParts, AssemblyError> {
        let start = match &options.angle {
            None => 180.0,
            Some(Scalar::Literal(v)) => *v,
            Some(_) => {
                return Err(param_error(
                    name,
                    "an assembly hinge's starting angle is a number (the slider drives it)",
                ));
            }
        };
        let mut components = vec![HingeComponentDesc::LeafA, HingeComponentDesc::LeafB];
        if kind == HingeKindDesc::Butt {
            components.push(HingeComponentDesc::Pin);
        }
        let mut parts = Vec::new();
        for component in components {
            let part_name = format!("{name}_{}", component.spelling());
            let mut model = Model::new();
            let flat = HingeOptions {
                angle: Some(Scalar::Literal(180.0)),
                only: Some(component),
                ..options.clone()
            };
            model
                .hinge(kind, length.clone(), &flat, &part_name)
                .map_err(|source| AssemblyError::Script {
                    definition: part_name.clone(),
                    source,
                })?;
            let id = self.add_part(&part_name, model)?;
            let at = self.place(Parent::Root, id, Some(&part_name), Placement::identity())?;
            parts.push((component, id, part_name, at));
        }
        let leaf_a = format!("{name}_leaf_a");
        for (component, _, part, _) in &parts {
            if *component == HingeComponentDesc::LeafA {
                continue;
            }
            self.add_mate(Mate {
                kind: MateKind::Concentric,
                a: datum(part, Datum::Z),
                b: datum(&leaf_a, Datum::Z),
                drive: None,
            })?;
            self.add_mate(Mate {
                kind: MateKind::Coincident { aligned: true },
                a: datum(part, Datum::Xy),
                b: datum(&leaf_a, Datum::Xy),
                drive: None,
            })?;
        }
        let angle = format!("{name}_angle");
        self.add_param(&angle, start, Unit::Degrees)?;
        self.set_param_bounds(&angle, 0.0, 360.0)?;
        // Leaf B built flat sits at 180°; turning it by (180° − θ) about
        // the pin axis opens the hinge to θ.
        self.add_mate(Mate {
            kind: MateKind::Angle {
                degrees: 180.0 - start,
                about: Some(datum(&leaf_a, Datum::Z)),
            },
            a: datum(&leaf_a, Datum::X),
            b: datum(&format!("{name}_leaf_b"), Datum::X),
            drive: None,
        })?;
        let index = self.mates.len() - 1;
        self.drive_mate(&angle, index, -1.0, 180.0)?;
        Ok(HingeParts { parts, angle })
    }

    /// Sweep assembly parameter `param` from `from` to `to` in steps of
    /// `step`, checking interference at each value, and refine the first
    /// contact by bisection to `resolution`. Contact is an overlap, or a
    /// touch between occurrences that were apart at `from`. The
    /// parameter is restored afterwards.
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
        let Some(original) = self
            .params
            .iter()
            .find(|p| p.name == param)
            .map(|p| p.display_value())
        else {
            return Err(param_error(param, "no assembly parameter with this name"));
        };
        if !(from.is_finite() && to.is_finite() && step > 0.0 && resolution > 0.0) {
            return Err(param_error(
                param,
                "the swing needs finite ends and positive steps",
            ));
        }
        let result = self.swing_inner(param, from, to, step, resolution);
        let restore = self.set_assembly_param(param, original);
        let report = result?;
        restore?;
        Ok(report)
    }

    fn swing_inner(
        &mut self,
        param: &str,
        from: f64,
        to: f64,
        step: f64,
        resolution: f64,
    ) -> Result<SwingReport, AssemblyError> {
        let first = self.swing_sample(param, from, None)?;
        let baseline: Vec<(String, String)> = first.touching.clone();
        let mut samples = vec![SwingSample {
            value: from,
            min_distance: first.min_distance,
            contact: first.contact,
        }];
        let mut first_contact = first.contact.then_some(from);
        let dir = if to >= from { 1.0 } else { -1.0 };
        let count = ((to - from).abs() / step).ceil() as usize;
        let mut last_clear = from;
        for k in 1..=count {
            if first_contact.is_some() {
                break;
            }
            let v = if k == count {
                to
            } else {
                from + dir * step * k as f64
            };
            let s = self.swing_sample(param, v, Some(&baseline))?;
            samples.push(SwingSample {
                value: v,
                min_distance: s.min_distance,
                contact: s.contact,
            });
            if s.contact {
                // Bisect between the last clear value and this one.
                let (mut clear, mut hit) = (last_clear, v);
                while (hit - clear).abs() > resolution {
                    let mid = 0.5 * (clear + hit);
                    if self.swing_sample(param, mid, Some(&baseline))?.contact {
                        hit = mid;
                    } else {
                        clear = mid;
                    }
                }
                first_contact = Some(hit);
            } else {
                last_clear = v;
            }
        }
        let closest = samples
            .iter()
            .filter(|s| !s.contact)
            .min_by(|a, b| a.min_distance.total_cmp(&b.min_distance));
        Ok(SwingReport {
            param: param.to_string(),
            first_contact,
            min_clearance: closest.map(|s| s.min_distance),
            min_clearance_at: closest.map(|s| s.value),
            samples,
        })
    }

    fn swing_sample(
        &mut self,
        param: &str,
        value: f64,
        baseline: Option<&[(String, String)]>,
    ) -> Result<Sample, AssemblyError> {
        self.set_assembly_param(param, value)?;
        // Report every pair the bounds bring within reach: clearance the
        // size of the assembly's bounds is generous; the BVH prunes.
        let report = self.interference(&InterferenceOptions {
            clearance: f64::MAX / 4.0,
            ..InterferenceOptions::default()
        })?;
        let mut min_distance = f64::INFINITY;
        let mut contact = false;
        let mut touching = Vec::new();
        for pair in &report.pairs {
            match &pair.status {
                PairStatus::Clear { distance, .. } => min_distance = min_distance.min(*distance),
                PairStatus::Touching { .. } => {
                    let key = (pair.a.clone(), pair.b.clone());
                    if baseline.is_some_and(|b| !b.contains(&key)) {
                        contact = true;
                    }
                    touching.push(key);
                }
                PairStatus::Interfering { .. } => contact = true,
            }
        }
        Ok(Sample {
            min_distance,
            contact,
            touching,
        })
    }
}

struct Sample {
    min_distance: f64,
    contact: bool,
    touching: Vec<(String, String)>,
}

/// One value of a swing.
#[derive(Debug, Clone, Serialize)]
pub struct SwingSample {
    /// The parameter value.
    pub value: f64,
    /// The smallest distance between occurrences apart there.
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
    /// The smallest clearance seen among collision-free samples.
    pub min_clearance: Option<f64>,
    /// Where that clearance was seen.
    pub min_clearance_at: Option<f64>,
    /// Every sample taken on the sweep (bisection probes excluded).
    pub samples: Vec<SwingSample>,
}
