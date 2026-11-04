//! `thicken`, `offset_solid` and `draft` on a body (spec §5 preview):
//!
//! ```rhai
//! box = box.offset_solid(1.0, #{ join: "intersection" });
//! sheet = sheet.thicken(0.5);
//! box = box.draft(faces().facing("x"), 5.0, xy(), #{ label: "taper" });
//! ```
//!
//! `draft`'s neutral plane is a face query (`faces()…` or
//! `on_face(…)`: its outward normal pulls) or a plane (`xy()`, `xz()`,
//! `yz()`, `xy_at(z)` and the like with a literal offset, `plane_at(o,
//! n)`: its normal pulls).

use popkern_history::feature::BasePlane;
use popkern_script::{FaceSelector, NeutralSpec, OffsetOptions, Plane, Scalar};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::{IntoScriptResult, argument_error};

const OFFSET_OPTIONS: &[&str] = &["label", "join", "trim"];
const DRAFT_OPTIONS: &[&str] = &["label"];

fn offset_options(
    op: &'static str,
    value: &Dynamic,
) -> Result<(OffsetOptions, String), Box<EvalAltResult>> {
    let options = Options::parse(op, value, OFFSET_OPTIONS)?;
    let mut offset = OffsetOptions::default();
    if let Some(j) = options.with("join", args::join_style)? {
        offset.join = j;
    }
    if let Some(t) = options.with("trim", args::trim_policy)? {
        offset.trim = t;
    }
    Ok((offset, options.label()?))
}

/// A draft's neutral plane from a face query or a plane value.
fn neutral(value: &Dynamic) -> Result<NeutralSpec, Box<EvalAltResult>> {
    if let Some(faces) = value.clone().try_cast::<FaceSelector>() {
        return Ok(NeutralSpec::Face(faces));
    }
    let Some(plane) = value.clone().try_cast::<Plane>() else {
        return Err(argument_error(
            "neutral",
            "expects a face query (faces()…, on_face(…)) or a plane (xy(), plane_at(o, n), …)",
        ));
    };
    let base = |base: BasePlane| match base {
        BasePlane::XY => [0.0, 0.0, 1.0],
        BasePlane::XZ => [0.0, 1.0, 0.0],
        BasePlane::YZ => [1.0, 0.0, 0.0],
    };
    let at = |normal: [f64; 3], d: f64| NeutralSpec::Plane {
        origin: normal.map(|c| c * d),
        normal,
    };
    match plane {
        Plane::XY => Ok(at(base(BasePlane::XY), 0.0)),
        Plane::XZ => Ok(at(base(BasePlane::XZ), 0.0)),
        Plane::YZ => Ok(at(base(BasePlane::YZ), 0.0)),
        Plane::Offset(b, Scalar::Literal(d)) => Ok(at(base(b), d)),
        Plane::OnFace(faces) => Ok(NeutralSpec::Face(faces)),
        Plane::Custom { origin, normal } => Ok(NeutralSpec::Plane { origin, normal }),
        Plane::Offset(..) | Plane::OnFaceOffset(..) => Err(argument_error(
            "neutral",
            "a draft's neutral plane is a face, a base plane at a literal offset, or plane_at(o, n)",
        )),
    }
}

impl ScriptBody {
    fn offset_solid(
        &self,
        distance: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let distance = args::scalar("distance", &distance)?;
        let (offset, label) = offset_options("offset_solid", options)?;
        let body = self
            .body
            .offset_solid_with(&mut self.model.borrow_mut(), distance, offset, &label)
            .script()?;
        Ok(self.served(body))
    }

    fn thicken(
        &self,
        thickness: Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let thickness = args::scalar("thickness", &thickness)?;
        let (offset, label) = offset_options("thicken", options)?;
        let body = self
            .body
            .thicken_with(&mut self.model.borrow_mut(), thickness, offset, &label)
            .script()?;
        Ok(self.served(body))
    }

    fn draft(
        &self,
        faces: FaceSelector,
        angle: Dynamic,
        neutral_plane: &Dynamic,
        options: &Dynamic,
    ) -> Result<ScriptBody, Box<EvalAltResult>> {
        let angle = args::scalar("angle", &angle)?;
        let neutral_plane = neutral(neutral_plane)?;
        let options = Options::parse("draft", options, DRAFT_OPTIONS)?;
        let body = self
            .body
            .draft_with(
                &mut self.model.borrow_mut(),
                faces,
                angle,
                neutral_plane,
                &options.label()?,
            )
            .script()?;
        Ok(self.served(body))
    }
}

/// Register `offset_solid`, `thicken` and `draft` on `Body`.
pub fn register(engine: &mut Engine) {
    engine
        .register_fn("offset_solid", |b: &mut ScriptBody, d: Dynamic| {
            b.offset_solid(d, &Dynamic::UNIT)
        })
        .register_fn(
            "offset_solid",
            |b: &mut ScriptBody, d: Dynamic, o: Dynamic| b.offset_solid(d, &o),
        )
        .register_fn("thicken", |b: &mut ScriptBody, t: Dynamic| {
            b.thicken(t, &Dynamic::UNIT)
        })
        .register_fn("thicken", |b: &mut ScriptBody, t: Dynamic, o: Dynamic| {
            b.thicken(t, &o)
        })
        .register_fn(
            "draft",
            |b: &mut ScriptBody, faces: FaceSelector, angle: Dynamic, plane: Dynamic| {
                b.draft(faces, angle, &plane, &Dynamic::UNIT)
            },
        )
        .register_fn(
            "draft",
            |b: &mut ScriptBody,
             faces: FaceSelector,
             angle: Dynamic,
             plane: Dynamic,
             o: Dynamic| { b.draft(faces, angle, &plane, &o) },
        );
}
