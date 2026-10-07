//! A seam hinge between two parts (PK-712): a box and its lid placed in
//! an assembly, hinged along the seam where they meet, so the lid opens
//! by an assembly parameter — a slider — and swings like any mated part.
//!
//! Each part's document records its own knuckles (`Feature::SeamHinge`,
//! owning that part's copy of the hinge parameters), with the seam
//! given in the part's frame. The lid is mated to the box about the pin
//! axis (a literal axis entity: no face or origin datum carries it),
//! held along it by a plane across it, and turned by an angle mate the
//! open-angle parameter drives.

use popkern_core::geometry::{Point3, Vec3};
use popkern_history::feature_hinge::{HingeComponentDesc, HingeKindDesc, SeamDesc};
use popkern_history::unit::Unit;
use popkern_script::Model;
use popkern_script::body_seam_hinge::SeamHingeOptions;
use popkern_script::scalar::Scalar;

use crate::assembly::{Assembly, InstanceRef, Parent};
use crate::error::AssemblyError;
use crate::mates::{Mate, MateEntity, MateKind, MateRef};
use crate::params::param_error;
use crate::placement::Placement;

/// What [`Assembly::add_seam_hinge`] made.
#[derive(Debug, Clone)]
pub struct SeamHingeParts {
    /// The open-angle parameter (degrees, 0 = closed).
    pub angle: String,
    /// A butt hinge's pin part occurrence.
    pub pin: Option<InstanceRef>,
}

fn arr(p: Point3) -> [f64; 3] {
    [p.x, p.y, p.z]
}

fn arrv(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}

impl Assembly {
    /// The single body of the part at `at` (a root instance), its
    /// placement, and the world centre of its bounds.
    fn hinge_part(
        &mut self,
        at: InstanceRef,
    ) -> Result<(String, crate::DefinitionId, Placement, Point3), AssemblyError> {
        let instance = self.instance(at)?.clone();
        if at.parent != Parent::Root {
            return Err(param_error(
                &instance.name,
                "a seam hinge joins top-level occurrences",
            ));
        }
        let model = self.part_model_mut(instance.definition)?;
        let bodies = model.alive_bodies();
        let [body] = bodies[..] else {
            return Err(param_error(
                &instance.name,
                &format!(
                    "a seam hinge's part has one body (this one has {})",
                    bodies.len()
                ),
            ));
        };
        let brep = model
            .brep_of(body)
            .map_err(|source| AssemblyError::Script {
                definition: instance.name.clone(),
                source,
            })?;
        let mut lo = [f64::MAX; 3];
        let mut hi = [f64::MIN; 3];
        for v in brep.vertices.values() {
            let p = [v.position.x, v.position.y, v.position.z];
            for k in 0..3 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        let centre = Point3::new(
            0.5 * (lo[0] + hi[0]),
            0.5 * (lo[1] + hi[1]),
            0.5 * (lo[2] + hi[2]),
        );
        Ok((
            instance.name.clone(),
            instance.definition,
            instance.placement,
            instance.placement.apply_point(&centre),
        ))
    }

    /// Hinge the part at `a` (the box) and the part at `b` (the lid)
    /// along the seam line through `origin` in `direction` (world),
    /// `length` long from `origin`. Records each part's knuckles in its
    /// own document, adds `<name>_angle` (degrees, 0 closed, limits
    /// [0, 300]) and mates the lid to the box about the pin axis; a butt
    /// hinge's pin becomes part `<name>_pin`.
    ///
    /// # Errors
    ///
    /// [`AssemblyError::Param`] for parts with other than one body, a
    /// zero direction or bodies along the seam; the hinge's typed
    /// declines (as [`AssemblyError::Script`]); the solve's errors.
    #[allow(clippy::too_many_arguments)]
    pub fn add_seam_hinge(
        &mut self,
        name: &str,
        kind: HingeKindDesc,
        a: InstanceRef,
        b: InstanceRef,
        origin: Point3,
        direction: Vec3,
        length: Scalar,
        options: &SeamHingeOptions,
        start_angle: f64,
    ) -> Result<SeamHingeParts, AssemblyError> {
        let (name_a, def_a, place_a, centre_a) = self.hinge_part(a)?;
        let (name_b, def_b, place_b, centre_b) = self.hinge_part(b)?;
        let axis = direction
            .normalized()
            .ok_or_else(|| param_error(name, "the seam direction has zero length"))?;
        let across = centre_b - centre_a;
        let up = (across - axis * across.dot(&axis))
            .normalized()
            .ok_or_else(|| param_error(name, "the two parts' centres lie along the seam"))?;
        let mut outward = axis.cross(&up);
        let mid = Point3::new(
            0.5 * (centre_a.x + centre_b.x),
            0.5 * (centre_a.y + centre_b.y),
            0.5 * (centre_a.z + centre_b.z),
        );
        if (origin - mid).dot(&outward) < 0.0 {
            outward = -outward;
        }
        // The seam in a part's own frame.
        let local = |p: &Placement| {
            let inv = p.inverse();
            let t = inv.transform();
            SeamDesc {
                origin: arr(t.apply_point(&origin)),
                axis: arrv(t.apply_direction(&axis)),
                outward: arrv(t.apply_direction(&outward)),
                up: arrv(t.apply_direction(&up)),
            }
        };
        let record = |asm: &mut Assembly,
                      def,
                      part: &str,
                      component,
                      seam: SeamDesc|
         -> Result<(), AssemblyError> {
            let model = asm.part_model_mut(def)?;
            let body = model.alive_bodies().first().copied();
            model
                .seam_hinge_owned(kind, component, body, seam, length.clone(), options, name)
                .map_err(|source| AssemblyError::Script {
                    definition: part.to_string(),
                    source,
                })?;
            Ok(())
        };
        record(
            self,
            def_a,
            &name_a,
            HingeComponentDesc::LeafA,
            local(&place_a),
        )?;
        record(
            self,
            def_b,
            &name_b,
            HingeComponentDesc::LeafB,
            local(&place_b),
        )?;
        // The pin axis sits outward of the seam by the knuckle radius
        // plus the clearance (the values just recorded in the box).
        let offset = {
            let model = self.part_model_mut(def_a)?;
            let value = |role: &str| {
                model
                    .params()
                    .into_iter()
                    .find(|p| p.name == format!("{name}.{role}"))
                    .map(|p| p.value)
            };
            let missing = || param_error(name, "the recorded hinge has no parameters");
            let (l, c) = (
                value("length").ok_or_else(missing)?,
                value("clearance").ok_or_else(missing)?,
            );
            let r = match value("knuckle_diameter") {
                Some(d) => d / 2.0,
                None => {
                    value("pin_diameter").ok_or_else(missing)? / 2.0
                        + c
                        + value("wall").ok_or_else(missing)?
                }
            };
            r + c.max(l * 1e-4)
        };
        let pin_point = origin + outward * offset;
        // Literal entities in each part's frame.
        let axis_of = |p: &Placement, part: &str| {
            let t = p.inverse();
            let t = t.transform();
            MateRef {
                instance: part.to_string(),
                entity: MateEntity::Axis {
                    origin: arr(t.apply_point(&pin_point)),
                    direction: arrv(t.apply_direction(&axis)),
                },
            }
        };
        let plane_of = |p: &Placement, part: &str, normal: Vec3| {
            let t = p.inverse();
            let t = t.transform();
            MateRef {
                instance: part.to_string(),
                entity: MateEntity::Plane {
                    origin: arr(t.apply_point(&origin)),
                    normal: arrv(t.apply_direction(&normal)),
                },
            }
        };
        let mut pin = None;
        if kind == HingeKindDesc::Butt {
            let pin_name = format!("{name}_pin");
            let mut model = Model::new();
            let seam = local(&Placement::identity());
            model
                .seam_hinge_owned(
                    kind,
                    HingeComponentDesc::Pin,
                    None,
                    seam,
                    length.clone(),
                    options,
                    name,
                )
                .map_err(|source| AssemblyError::Script {
                    definition: pin_name.clone(),
                    source,
                })?;
            let id = self.add_part(&pin_name, model)?;
            let at = self.place(Parent::Root, id, Some(&pin_name), Placement::identity())?;
            let ident = Placement::identity();
            self.add_mate(Mate {
                kind: MateKind::Concentric,
                a: axis_of(&ident, &pin_name),
                b: axis_of(&place_a, &name_a),
                drive: None,
            })?;
            self.add_mate(Mate {
                kind: MateKind::Coincident { aligned: true },
                a: plane_of(&ident, &pin_name, axis),
                b: plane_of(&place_a, &name_a, axis),
                drive: None,
            })?;
            pin = Some(at);
        }
        self.add_mate(Mate {
            kind: MateKind::Concentric,
            a: axis_of(&place_b, &name_b),
            b: axis_of(&place_a, &name_a),
            drive: None,
        })?;
        self.add_mate(Mate {
            kind: MateKind::Coincident { aligned: true },
            a: plane_of(&place_b, &name_b, axis),
            b: plane_of(&place_a, &name_a, axis),
            drive: None,
        })?;
        let angle = format!("{name}_angle");
        self.add_param(&angle, start_angle, Unit::Degrees)?;
        self.set_param_bounds(&angle, 0.0, 300.0)?;
        // Opening turns the lid so the side away from the hinge rises:
        // the rotation taking −outward toward up, about the axis.
        let sense = if (-outward).cross(&up).dot(&axis) > 0.0 {
            1.0
        } else {
            -1.0
        };
        self.add_mate(Mate {
            kind: MateKind::Angle {
                degrees: sense * start_angle,
                about: Some(axis_of(&place_a, &name_a)),
            },
            a: plane_of(&place_a, &name_a, up),
            b: plane_of(&place_b, &name_b, up),
            drive: None,
        })?;
        let index = self.mates.len() - 1;
        self.drive_mate(&angle, index, sense, 0.0)?;
        Ok(SeamHingeParts { angle, pin })
    }
}
