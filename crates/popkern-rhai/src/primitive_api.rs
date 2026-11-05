//! Primitive verbs (PK-1083): `m.box(x, y, z)`, `m.cylinder(r, h)`,
//! `m.cone(r1, r2, h)`, `m.sphere(r)`, `m.prism(sides, r, h)`,
//! `m.pyramid(sides, r1, r2, h)`, `m.wedge(x, y, z, top_x)` and
//! `m.torus(major, minor)`, each with an options map: `label`, the
//! placement `at`, `axis`, `ref`, and per shape `centred` (box), `lat`
//! (sphere: `[lo, hi]` degrees), `angle` (a sector of a cylinder,
//! cone, sphere or torus, degrees), `flats` (prism and pyramid: the
//! radius runs to the flats).

use popkern_script::{Placement, RegularPolygon, Scalar};
use rhai::{Dynamic, Engine, EvalAltResult};

use crate::args::{self, Options};
use crate::body_api::ScriptBody;
use crate::error::IntoScriptResult;
use crate::model_api::ScriptModel;

const BOX_OPTIONS: &[&str] = &["label", "at", "axis", "ref", "centred"];
const SECTOR_OPTIONS: &[&str] = &["label", "at", "axis", "ref", "angle"];
const SPHERE_OPTIONS: &[&str] = &["label", "at", "axis", "ref", "lat", "angle"];
const POLYGON_OPTIONS: &[&str] = &["label", "at", "axis", "ref", "flats"];
const PLACED: &[&str] = &["label", "at", "axis", "ref"];

fn regular_polygon(
    sides: &Dynamic,
    options: &Options,
) -> Result<RegularPolygon, Box<EvalAltResult>> {
    let sides = sides
        .as_int()
        .ok()
        .and_then(|n| u32::try_from(n).ok())
        .filter(|n| *n >= 3)
        .ok_or_else(|| crate::error::argument_error("sides", "an integer of at least 3"))?;
    let flats = options.with("flats", args::boolean)?.unwrap_or(false);
    Ok(RegularPolygon { sides, flats })
}

fn prism(
    m: &mut ScriptModel,
    sides: &Dynamic,
    r: &Dynamic,
    h: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("prism", o, POLYGON_OPTIONS)?;
    let polygon = regular_polygon(sides, &options)?;
    let [r, h] = scalars(["radius", "height"], [r, h])?;
    let body = m
        .borrow_mut()
        .prism(polygon, r, h, placement(&options)?, &options.label()?);
    served(m, body)
}

fn pyramid(
    m: &mut ScriptModel,
    sides: &Dynamic,
    radii: [&Dynamic; 2],
    h: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("pyramid", o, POLYGON_OPTIONS)?;
    let polygon = regular_polygon(sides, &options)?;
    let radii = scalars(["r1", "r2"], radii)?;
    let [h] = scalars(["height"], [h])?;
    let body = m
        .borrow_mut()
        .pyramid(polygon, radii, h, placement(&options)?, &options.label()?);
    served(m, body)
}

fn wedge(
    m: &mut ScriptModel,
    size: [&Dynamic; 3],
    top_x: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("wedge", o, PLACED)?;
    let size = scalars(["x", "y", "z"], size)?;
    let [t] = scalars(["top_x"], [top_x])?;
    let body = m
        .borrow_mut()
        .wedge(size, t, placement(&options)?, &options.label()?);
    served(m, body)
}

fn placement(options: &Options) -> Result<Placement, Box<EvalAltResult>> {
    let mut p = Placement::default();
    if let Some(at) = options.with("at", args::vec3)? {
        p.at = at;
    }
    if let Some(axis) = options.with("axis", args::vec3)? {
        p.axis = axis;
    }
    if let Some(reference) = options.with("ref", args::vec3)? {
        p.reference = reference;
    }
    Ok(p)
}

fn scalars<const N: usize>(
    names: [&str; N],
    values: [&Dynamic; N],
) -> Result<[Scalar; N], Box<EvalAltResult>> {
    let mut out = Vec::with_capacity(N);
    for (n, v) in names.iter().zip(values) {
        out.push(args::scalar(n, v)?);
    }
    out.try_into()
        .map_err(|_| crate::error::argument_error("dimensions", "wrong count"))
}

fn served(
    m: &ScriptModel,
    body: Result<popkern_script::Body, popkern_script::ScriptError>,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    Ok(ScriptBody::new(m.clone(), body.script()?))
}

fn cuboid(
    m: &mut ScriptModel,
    x: &Dynamic,
    y: &Dynamic,
    z: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("box", o, BOX_OPTIONS)?;
    let size = scalars(["x", "y", "z"], [x, y, z])?;
    let centred = options.with("centred", args::boolean)?.unwrap_or(false);
    let body = m
        .borrow_mut()
        .cuboid(size, centred, placement(&options)?, &options.label()?);
    served(m, body)
}

fn cylinder(
    m: &mut ScriptModel,
    r: &Dynamic,
    h: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("cylinder", o, SECTOR_OPTIONS)?;
    let [r, h] = scalars(["radius", "height"], [r, h])?;
    let angle = options.with("angle", args::scalar)?;
    let body = m
        .borrow_mut()
        .cylinder(r, h, angle, placement(&options)?, &options.label()?);
    served(m, body)
}

fn cone(
    m: &mut ScriptModel,
    r1: &Dynamic,
    r2: &Dynamic,
    h: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("cone", o, SECTOR_OPTIONS)?;
    let [r1, r2, h] = scalars(["r1", "r2", "height"], [r1, r2, h])?;
    let angle = options.with("angle", args::scalar)?;
    let body = m
        .borrow_mut()
        .cone(r1, r2, h, angle, placement(&options)?, &options.label()?);
    served(m, body)
}

fn sphere(m: &mut ScriptModel, r: &Dynamic, o: &Dynamic) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("sphere", o, SPHERE_OPTIONS)?;
    let [r] = scalars(["radius"], [r])?;
    let lat = options.with("lat", |a, v| {
        let pair = args::point(a, v, 2)?;
        Ok([Scalar::Literal(pair[0]), Scalar::Literal(pair[1])])
    })?;
    let angle = options.with("angle", args::scalar)?;
    let body = m
        .borrow_mut()
        .sphere(r, lat, angle, placement(&options)?, &options.label()?);
    served(m, body)
}

fn torus(
    m: &mut ScriptModel,
    major: &Dynamic,
    minor: &Dynamic,
    o: &Dynamic,
) -> Result<ScriptBody, Box<EvalAltResult>> {
    let options = Options::parse("torus", o, SECTOR_OPTIONS)?;
    let [major, minor] = scalars(["major", "minor"], [major, minor])?;
    let angle = options.with("angle", args::scalar)?;
    let body = m
        .borrow_mut()
        .torus(major, minor, angle, placement(&options)?, &options.label()?);
    served(m, body)
}

/// Register the primitive verbs.
pub fn register(engine: &mut Engine) {
    let none = Dynamic::UNIT;
    let n = none.clone();
    engine.register_fn(
        "box",
        move |m: &mut ScriptModel, x: Dynamic, y: Dynamic, z: Dynamic| cuboid(m, &x, &y, &z, &n),
    );
    engine.register_fn(
        "box",
        |m: &mut ScriptModel, x: Dynamic, y: Dynamic, z: Dynamic, o: Dynamic| {
            cuboid(m, &x, &y, &z, &o)
        },
    );
    let n = none.clone();
    engine.register_fn(
        "cylinder",
        move |m: &mut ScriptModel, r: Dynamic, h: Dynamic| cylinder(m, &r, &h, &n),
    );
    engine.register_fn(
        "cylinder",
        |m: &mut ScriptModel, r: Dynamic, h: Dynamic, o: Dynamic| cylinder(m, &r, &h, &o),
    );
    let n = none.clone();
    engine.register_fn(
        "cone",
        move |m: &mut ScriptModel, r1: Dynamic, r2: Dynamic, h: Dynamic| cone(m, &r1, &r2, &h, &n),
    );
    engine.register_fn(
        "cone",
        |m: &mut ScriptModel, r1: Dynamic, r2: Dynamic, h: Dynamic, o: Dynamic| {
            cone(m, &r1, &r2, &h, &o)
        },
    );
    let n = none.clone();
    engine.register_fn("sphere", move |m: &mut ScriptModel, r: Dynamic| {
        sphere(m, &r, &n)
    });
    engine.register_fn("sphere", |m: &mut ScriptModel, r: Dynamic, o: Dynamic| {
        sphere(m, &r, &o)
    });
    let n = none.clone();
    engine.register_fn(
        "prism",
        move |m: &mut ScriptModel, s: Dynamic, r: Dynamic, h: Dynamic| prism(m, &s, &r, &h, &n),
    );
    engine.register_fn(
        "prism",
        |m: &mut ScriptModel, s: Dynamic, r: Dynamic, h: Dynamic, o: Dynamic| {
            prism(m, &s, &r, &h, &o)
        },
    );
    let n = none.clone();
    engine.register_fn(
        "pyramid",
        move |m: &mut ScriptModel, s: Dynamic, a: Dynamic, b: Dynamic, h: Dynamic| {
            pyramid(m, &s, [&a, &b], &h, &n)
        },
    );
    engine.register_fn(
        "pyramid",
        |m: &mut ScriptModel, s: Dynamic, a: Dynamic, b: Dynamic, h: Dynamic, o: Dynamic| {
            pyramid(m, &s, [&a, &b], &h, &o)
        },
    );
    let n = none.clone();
    engine.register_fn(
        "wedge",
        move |m: &mut ScriptModel, x: Dynamic, y: Dynamic, z: Dynamic, t: Dynamic| {
            wedge(m, [&x, &y, &z], &t, &n)
        },
    );
    engine.register_fn(
        "wedge",
        |m: &mut ScriptModel, x: Dynamic, y: Dynamic, z: Dynamic, t: Dynamic, o: Dynamic| {
            wedge(m, [&x, &y, &z], &t, &o)
        },
    );
    let n = none;
    engine.register_fn(
        "torus",
        move |m: &mut ScriptModel, a: Dynamic, b: Dynamic| torus(m, &a, &b, &n),
    );
    engine.register_fn(
        "torus",
        |m: &mut ScriptModel, a: Dynamic, b: Dynamic, o: Dynamic| torus(m, &a, &b, &o),
    );
}
