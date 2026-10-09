//! A hinge as an assembly (PK-701, design `docs/architecture/hinges.md`
//! §4.4): its leaves (and pin) become part definitions placed in the
//! hinge frame and mated about the pin axis, with an assembly parameter
//! `<name>_angle` driving the open angle — a slider that moves leaf B
//! without rebuilding any part. Its swing sweep (where the motion first
//! collides) is `swing.rs`.

use popkern_history::feature_hinge::{HingeComponentDesc, HingeKindDesc};
use popkern_history::unit::Unit;
use popkern_script::Model;
use popkern_script::body_hinge::HingeOptions;
use popkern_script::scalar::Scalar;

use crate::assembly::{Assembly, DefinitionId, InstanceRef, Parent};
use crate::error::AssemblyError;
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
}
