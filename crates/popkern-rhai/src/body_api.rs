//! The `Body` value: every feature op of spec §5 as a method that
//! records into the shared model and returns the served body, so
//! calls chain in timeline order.

use popkern_history::feature::{BoolOpKind, OperationMode};
use popkern_script::body_sweep::HelixParams;
use popkern_script::{Body, EdgeSelector, End, FaceSelector, RevolveEnd, Scalar};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::error::{IntoScriptResult, argument_error};
use crate::feature_api::ScriptFeature;
use crate::model_api::ScriptModel;
use crate::profile_api::{ScriptAxis, ScriptEnd};
use crate::recipe_api::{chamfer_options, fillet_options, shell_options};

/// A body of the document, bound to its model.
#[derive(Debug, Clone)]
pub struct ScriptBody {
    pub(crate) model: ScriptModel,
    pub(crate) body: Body,
}

const EXTRUDE_OPTIONS: &[&str] = &["label", "mode"];
const LABEL_ONLY: &[&str] = &["label"];
const SWEEP_OPTIONS: &[&str] = &["label", "mode", "twist"];
const HELICAL_OPTIONS: &[&str] = &["label", "mode", "taper_to", "start_angle", "left_hand"];

impl ScriptBody {
    pub(crate) fn new(model: ScriptModel, body: Body) -> ScriptBody {
        ScriptBody { model, body }
    }

    fn served(&self, body: Body) -> ScriptBody {
        ScriptBody::new(self.model.clone(), body)
    }

    fn extrude(&self, end: ScriptEnd, options: &Dynamic) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("extrude", options, EXTRUDE_OPTIONS)?;
        let mode = options
            .with("mode", args::mode)?
            .unwrap_or(OperationMode::NewBody);
        let end = match end {
            ScriptEnd::Extrude(end) => end,
            ScriptEnd::FullCircle => {
                return Err(argument_error("end", "full_circle() is a revolve end"));
            }
        };
        let body = self
            .body
            .extrude_with(&mut self.model.borrow_mut(), end, mode, &options.label()?)
            .script()?;
        Ok(self.served(body))
    }

    fn revolve(
        &self,
        axis: ScriptAxis,
        end: ScriptEnd,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("revolve", options, LABEL_ONLY)?;
        let end = match end {
            ScriptEnd::FullCircle => RevolveEnd::FullCircle,
            ScriptEnd::Extrude(End::Blind(angle)) => RevolveEnd::Blind(angle),
            ScriptEnd::Extrude(End::MidPlane(half)) => RevolveEnd::MidPlane(half),
            ScriptEnd::Extrude(_) => {
                return Err(argument_error(
                    "end",
                    "revolve ends are blind(angle), mid_plane(half_angle) or full_circle()",
                ));
            }
        };
        let body = self
            .body
            .revolve_with(
                &mut self.model.borrow_mut(),
                axis.spec(),
                end,
                &options.label()?,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn fillet(
        &self,
        edges: EdgeSelector,
        radius: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let radius = args::scalar("radius", &radius)?;
        let (mut recipe, law, label) = fillet_options(options)?;
        recipe.end_setbacks = self
            .body
            .setback_ends(
                &mut self.model.borrow_mut(),
                crate::recipe_api::setback_end_selectors(options)?,
            )
            .script()?;
        let body = self
            .body
            .fillet_law(
                &mut self.model.borrow_mut(),
                edges,
                radius,
                law,
                recipe,
                &label,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn chamfer(
        &self,
        edges: EdgeSelector,
        distance: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let distance = args::scalar("distance", &distance)?;
        let (kind, label) = chamfer_options(options)?;
        let body = self
            .body
            .chamfer_with(&mut self.model.borrow_mut(), edges, distance, kind, &label)
            .script()?;
        Ok(self.served(body))
    }

    fn shell(
        &self,
        faces: FaceSelector,
        thickness: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let thickness = args::scalar("thickness", &thickness)?;
        let (shell, label) = shell_options(options)?;
        let body = self
            .body
            .shell_with(
                &mut self.model.borrow_mut(),
                faces,
                thickness,
                shell,
                &label,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn boolean(
        &self,
        tool: ScriptBody,
        op: BoolOpKind,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("boolean", options, LABEL_ONLY)?;
        let body = self
            .body
            .boolean_with(
                &mut self.model.borrow_mut(),
                op,
                tool.body,
                &options.label()?,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn translate(
        &self,
        dx: Dynamic,
        dy: Dynamic,
        dz: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("translate", options, LABEL_ONLY)?;
        let body = self
            .body
            .translate_with(
                &mut self.model.borrow_mut(),
                args::scalar("dx", &dx)?,
                args::scalar("dy", &dy)?,
                args::scalar("dz", &dz)?,
                &options.label()?,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn rotate(
        &self,
        axis: ScriptAxis,
        angle: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("rotate", options, LABEL_ONLY)?;
        let angle = args::scalar("angle", &angle)?;
        let body = self
            .body
            .rotate_with(
                &mut self.model.borrow_mut(),
                axis.spec(),
                angle,
                &options.label()?,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn sweep(
        &self,
        path_body: ScriptBody,
        path: EdgeSelector,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("sweep", options, SWEEP_OPTIONS)?;
        let mode = options
            .with("mode", args::mode)?
            .unwrap_or(OperationMode::NewBody);
        let twist = options
            .with("twist", args::scalar)?
            .unwrap_or_else(|| Scalar::from(0.0));
        let body = self
            .body
            .sweep_with(
                &mut self.model.borrow_mut(),
                path_body.body,
                path,
                twist,
                mode,
                &options.label()?,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn pipe(
        &self,
        path: EdgeSelector,
        radius: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let radius = args::scalar("radius", &radius)?;
        let options = Options::parse("pipe", options, EXTRUDE_OPTIONS)?;
        let mode = options
            .with("mode", args::mode)?
            .unwrap_or(OperationMode::NewBody);
        let body = self
            .body
            .pipe_with(
                &mut self.model.borrow_mut(),
                path,
                radius,
                mode,
                &options.label()?,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn helical_sweep(
        &self,
        axis: ScriptAxis,
        radius: Dynamic,
        pitch: Dynamic,
        turns: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let options = Options::parse("helical_sweep", options, HELICAL_OPTIONS)?;
        let mode = options
            .with("mode", args::mode)?
            .unwrap_or(OperationMode::NewBody);
        let helix = HelixParams {
            radius: args::scalar("radius", &radius)?,
            taper_to: options.with("taper_to", args::scalar)?,
            pitch: args::scalar("pitch", &pitch)?,
            turns: args::scalar("turns", &turns)?,
            start_angle: options
                .with("start_angle", args::scalar)?
                .unwrap_or_else(|| Scalar::from(0.0)),
            left_handed: options.with("left_hand", args::boolean)?.unwrap_or(false),
        };
        let body = self
            .body
            .helical_sweep_with(
                &mut self.model.borrow_mut(),
                axis.spec(),
                helix,
                mode,
                &options.label()?,
            )
            .script()?;
        Ok(self.served(body))
    }

    fn last(&self) -> Result<ScriptFeature, Box<EvalAltResult>> {
        let handle = self
            .body
            .last(&self.model.borrow_mut())
            .ok_or_else(|| argument_error("body", "body has no features yet"))?;
        Ok(ScriptFeature::new(self.model.clone(), handle))
    }
}

/// Register the body value.
pub fn register(engine: &mut Engine) {
    engine
        .register_type_with_name::<ScriptBody>("Body")
        .register_get("id", |b: &mut ScriptBody| b.body.id.0 as i64)
        .register_get("label", |b: &mut ScriptBody| {
            b.model
                .borrow_mut()
                .body_label(b.body)
                .map(Dynamic::from)
                .unwrap_or(Dynamic::UNIT)
        })
        .register_fn("is_alive", |b: &mut ScriptBody| {
            b.body.is_alive(&b.model.borrow_mut())
        })
        .register_fn("extrude", |b: &mut ScriptBody, end: ScriptEnd| {
            b.extrude(end, &Dynamic::UNIT)
        })
        .register_fn(
            "extrude",
            |b: &mut ScriptBody, end: ScriptEnd, options: Dynamic| b.extrude(end, &options),
        )
        .register_fn(
            "revolve",
            |b: &mut ScriptBody, axis: ScriptAxis, end: ScriptEnd| {
                b.revolve(axis, end, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "revolve",
            |b: &mut ScriptBody, axis: ScriptAxis, end: ScriptEnd, options: Dynamic| {
                b.revolve(axis, end, &options)
            },
        )
        .register_fn(
            "fillet",
            |b: &mut ScriptBody, edges: EdgeSelector, radius: Dynamic| {
                b.fillet(edges, radius, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "fillet",
            |b: &mut ScriptBody, edges: EdgeSelector, radius: Dynamic, options: Dynamic| {
                b.fillet(edges, radius, &options)
            },
        )
        .register_fn("fillet_each", |b: &mut ScriptBody, groups: rhai::Array| {
            b.fillet_each(groups, &Dynamic::UNIT)
        })
        .register_fn(
            "fillet_each",
            |b: &mut ScriptBody, groups: rhai::Array, options: Dynamic| {
                b.fillet_each(groups, &options)
            },
        )
        .register_fn(
            "chamfer",
            |b: &mut ScriptBody, edges: EdgeSelector, d: Dynamic| {
                b.chamfer(edges, d, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "chamfer",
            |b: &mut ScriptBody, edges: EdgeSelector, d: Dynamic, options: Dynamic| {
                b.chamfer(edges, d, &options)
            },
        )
        .register_fn(
            "shell",
            |b: &mut ScriptBody, faces: FaceSelector, t: Dynamic| b.shell(faces, t, &Dynamic::UNIT),
        )
        .register_fn(
            "shell",
            |b: &mut ScriptBody, faces: FaceSelector, t: Dynamic, options: Dynamic| {
                b.shell(faces, t, &options)
            },
        )
        .register_fn("union", |b: &mut ScriptBody, tool: ScriptBody| {
            b.boolean(tool, BoolOpKind::Union, &Dynamic::UNIT)
        })
        .register_fn(
            "union",
            |b: &mut ScriptBody, tool: ScriptBody, o: Dynamic| {
                b.boolean(tool, BoolOpKind::Union, &o)
            },
        )
        .register_fn("subtract", |b: &mut ScriptBody, tool: ScriptBody| {
            b.boolean(tool, BoolOpKind::Subtraction, &Dynamic::UNIT)
        })
        .register_fn(
            "subtract",
            |b: &mut ScriptBody, tool: ScriptBody, o: Dynamic| {
                b.boolean(tool, BoolOpKind::Subtraction, &o)
            },
        )
        .register_fn("intersect", |b: &mut ScriptBody, tool: ScriptBody| {
            b.boolean(tool, BoolOpKind::Intersection, &Dynamic::UNIT)
        })
        .register_fn(
            "intersect",
            |b: &mut ScriptBody, tool: ScriptBody, o: Dynamic| {
                b.boolean(tool, BoolOpKind::Intersection, &o)
            },
        )
        .register_fn(
            "translate",
            |b: &mut ScriptBody, dx: Dynamic, dy: Dynamic, dz: Dynamic| {
                b.translate(dx, dy, dz, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "translate",
            |b: &mut ScriptBody, dx: Dynamic, dy: Dynamic, dz: Dynamic, o: Dynamic| {
                b.translate(dx, dy, dz, &o)
            },
        )
        .register_fn(
            "rotate",
            |b: &mut ScriptBody, axis: ScriptAxis, angle: Dynamic| {
                b.rotate(axis, angle, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "rotate",
            |b: &mut ScriptBody, axis: ScriptAxis, angle: Dynamic, o: Dynamic| {
                b.rotate(axis, angle, &o)
            },
        )
        .register_fn(
            "sweep",
            |b: &mut ScriptBody, path_body: ScriptBody, path: EdgeSelector| {
                b.sweep(path_body, path, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "sweep",
            |b: &mut ScriptBody, path_body: ScriptBody, path: EdgeSelector, o: Dynamic| {
                b.sweep(path_body, path, &o)
            },
        )
        .register_fn(
            "pipe",
            |b: &mut ScriptBody, path: EdgeSelector, radius: Dynamic| {
                b.pipe(path, radius, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "pipe",
            |b: &mut ScriptBody, path: EdgeSelector, radius: Dynamic, o: Dynamic| {
                b.pipe(path, radius, &o)
            },
        )
        .register_fn(
            "helical_sweep",
            |b: &mut ScriptBody,
             axis: ScriptAxis,
             radius: Dynamic,
             pitch: Dynamic,
             turns: Dynamic| {
                b.helical_sweep(axis, radius, pitch, turns, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "helical_sweep",
            |b: &mut ScriptBody,
             axis: ScriptAxis,
             radius: Dynamic,
             pitch: Dynamic,
             turns: Dynamic,
             o: Dynamic| { b.helical_sweep(axis, radius, pitch, turns, &o) },
        )
        .register_fn("last", |b: &mut ScriptBody| b.last())
        .register_fn("to_string", |b: &mut ScriptBody| {
            format!("Body({})", b.body.id.0)
        });
}

#[path = "body_api_groups.rs"]
mod body_api_groups;
