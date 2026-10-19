//! The assembly values (design `docs/architecture/assemblies.md` §6.1,
//! `preview`): `assembly(name?)`, part definitions recorded by a
//! closure on a fresh model, sub-assemblies, and named instances placed
//! by chained rigid motions.
//!
//! ```rhai
//! let a = assembly("bearing");
//! let ball = a.part("ball", |m| { m.profile(xz()) … .revolve(…); });
//! for i in 0..8 { a.place(ball, "ball_" + i).translate(0.7 * cos(t), 0.7 * sin(t), 0.0); }
//! let stage = a.subassembly("stage", |s| { s.place(ball); });
//! a
//! ```

use std::cell::{RefCell, RefMut};
use std::rc::Rc;

use popkern_assembly::{Assembly, AssemblyError, DefinitionId, InstanceRef, Parent, Placement};
use popkern_core::geometry::{Point3, Vec3};
use popkern_history::deferred_certification::ContinuityCertification;
use popkern_script::Model;
use rhai::{Dynamic, Engine, EvalAltResult, FnPtr, Map, NativeCallContext};

use crate::args;
use crate::error::{argument_error, raise};
use crate::model_api::ScriptModel;
use crate::profile_api::ScriptAxis;

/// A parametric part family declared with `a.part(name, defaults, |m, p| …)`.
pub(crate) struct Family {
    pub(crate) name: String,
    pub(crate) defaults: std::collections::BTreeMap<String, popkern_assembly::ParamValue>,
    pub(crate) builder: FnPtr,
}

pub(crate) struct State {
    pub(crate) assembly: Assembly,
    pub(crate) certification: ContinuityCertification,
    pub(crate) families: Vec<Family>,
}

/// A shared handle to an assembly document, placing into `parent`
/// (the root, or the sub-assembly a `subassembly` closure fills).
#[derive(Clone)]
pub struct ScriptAssembly {
    pub(crate) state: Rc<RefCell<State>>,
    pub(crate) parent: Parent,
}

impl std::fmt::Debug for ScriptAssembly {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("Assembly")
    }
}

impl ScriptAssembly {
    /// A new, empty assembly whose part models certify under
    /// `certification`.
    pub fn new(name: &str, certification: ContinuityCertification) -> ScriptAssembly {
        ScriptAssembly {
            state: Rc::new(RefCell::new(State {
                assembly: Assembly::new(name),
                certification,
                families: Vec::new(),
            })),
            parent: Parent::Root,
        }
    }

    /// Wrap a loaded assembly.
    pub fn from_assembly(
        assembly: Assembly,
        certification: ContinuityCertification,
    ) -> ScriptAssembly {
        ScriptAssembly {
            state: Rc::new(RefCell::new(State {
                assembly,
                certification,
                families: Vec::new(),
            })),
            parent: Parent::Root,
        }
    }

    pub(crate) fn borrow_mut(&self) -> RefMut<'_, State> {
        self.state.borrow_mut()
    }

    /// Take the document out (the handle keeps an empty one).
    pub fn into_inner(self) -> Assembly {
        let mut state = self.state.borrow_mut();
        std::mem::replace(&mut state.assembly, Assembly::new(""))
    }

    /// Whether this handle is the document's root (not a sub-assembly
    /// handle a closure was given).
    pub fn is_root(&self) -> bool {
        self.parent == Parent::Root
    }
}

/// A part or sub-assembly definition handle.
#[derive(Clone)]
pub struct ScriptDefinition {
    pub(crate) owner: ScriptAssembly,
    pub(crate) id: DefinitionId,
    /// The family this definition is a variant of.
    pub(crate) family: Option<usize>,
}

impl std::fmt::Debug for ScriptDefinition {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Definition#{}", self.id.0)
    }
}

/// A placed instance handle.
#[derive(Clone)]
pub struct ScriptInstance {
    pub(crate) owner: ScriptAssembly,
    pub(crate) at: InstanceRef,
}

impl std::fmt::Debug for ScriptInstance {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Instance[{}]", self.at.index)
    }
}

/// An assembly error as a script runtime error (the shared payload).
pub(crate) fn assembly_error(error: AssemblyError) -> Box<EvalAltResult> {
    raise(error.payload())
}

pub(crate) trait IntoAssemblyResult<T> {
    fn assembly(self) -> Result<T, Box<EvalAltResult>>;
}

impl<T> IntoAssemblyResult<T> for Result<T, AssemblyError> {
    fn assembly(self) -> Result<T, Box<EvalAltResult>> {
        self.map_err(assembly_error)
    }
}

/// An error raised inside a part's closure, unwrapped from the
/// function-call frames to the raising call (so its line is the
/// closure's) and tagged with the part.
pub(crate) fn in_part(part: &str, error: Box<EvalAltResult>) -> Box<EvalAltResult> {
    match *crate::error::innermost(error) {
        EvalAltResult::ErrorRuntime(value, position) => {
            let value = match value.clone().try_cast::<Map>() {
                Some(mut map) => {
                    map.insert("part".into(), part.into());
                    if let Some(message) = map
                        .get("message")
                        .and_then(|m| m.clone().into_string().ok())
                    {
                        map.insert(
                            "message".into(),
                            format!("in part '{part}': {message}").into(),
                        );
                    }
                    Dynamic::from_map(map)
                }
                None => value,
            };
            Box::new(EvalAltResult::ErrorRuntime(value, position))
        }
        other => Box::new(other),
    }
}

/// Run `build` on a fresh model and return the model it recorded.
pub(crate) fn record_part(
    ctx: &NativeCallContext,
    owner: &ScriptAssembly,
    name: &str,
    build: &FnPtr,
    params: Option<Map>,
) -> Result<Model, Box<EvalAltResult>> {
    let certification = owner.borrow_mut().certification;
    let mut model = Model::new();
    model.set_continuity_certification(certification);
    let handle = ScriptModel::new(model);
    let result: Result<Dynamic, _> = match params {
        None => build.call_within_context(ctx, (handle.clone(),)),
        Some(p) => build.call_within_context(ctx, (handle.clone(), Dynamic::from_map(p))),
    };
    let _: Dynamic = result.map_err(|e| in_part(name, e))?;
    Ok(handle.into_inner())
}

fn part(
    ctx: NativeCallContext,
    a: &mut ScriptAssembly,
    name: &str,
    build: FnPtr,
) -> Result<ScriptDefinition, Box<EvalAltResult>> {
    let model = record_part(&ctx, a, name, &build, None)?;
    let id = a.borrow_mut().assembly.add_part(name, model).assembly()?;
    Ok(ScriptDefinition {
        owner: a.clone(),
        id,
        family: None,
    })
}

fn subassembly(
    ctx: NativeCallContext,
    a: &mut ScriptAssembly,
    name: &str,
    fill: FnPtr,
) -> Result<ScriptDefinition, Box<EvalAltResult>> {
    let id = a.borrow_mut().assembly.add_subassembly(name).assembly()?;
    let inner = ScriptAssembly {
        state: Rc::clone(&a.state),
        parent: Parent::Definition(id),
    };
    let _: Dynamic = fill
        .call_within_context(&ctx, (inner,))
        .map_err(|e| in_part(name, e))?;
    Ok(ScriptDefinition {
        owner: a.clone(),
        id,
        family: None,
    })
}

fn place(
    a: &mut ScriptAssembly,
    definition: ScriptDefinition,
    name: Option<&str>,
) -> Result<ScriptInstance, Box<EvalAltResult>> {
    if !Rc::ptr_eq(&a.state, &definition.owner.state) {
        return Err(argument_error(
            "definition",
            "this definition belongs to another assembly",
        ));
    }
    let at = a
        .borrow_mut()
        .assembly
        .place(a.parent, definition.id, name, Placement::identity())
        .assembly()?;
    Ok(ScriptInstance {
        owner: a.clone(),
        at,
    })
}

impl ScriptInstance {
    fn moved(&self, delta: Placement) -> Result<ScriptInstance, Box<EvalAltResult>> {
        let mut state = self.owner.borrow_mut();
        let current = state.assembly.instance(self.at).assembly()?.placement;
        state
            .assembly
            .set_placement(self.at, current.then(&delta))
            .assembly()?;
        Ok(self.clone())
    }

    fn name(&self) -> Result<String, Box<EvalAltResult>> {
        let state = self.owner.borrow_mut();
        Ok(state.assembly.instance(self.at).assembly()?.name.clone())
    }
}

fn translate(
    i: &mut ScriptInstance,
    x: Dynamic,
    y: Dynamic,
    z: Dynamic,
) -> Result<ScriptInstance, Box<EvalAltResult>> {
    i.moved(Placement::translation(
        args::float("dx", &x)?,
        args::float("dy", &y)?,
        args::float("dz", &z)?,
    ))
}

fn rotate(
    i: &mut ScriptInstance,
    axis: ScriptAxis,
    degrees: Dynamic,
) -> Result<ScriptInstance, Box<EvalAltResult>> {
    let ScriptAxis::Literal { origin, direction } = axis else {
        return Err(argument_error(
            "axis",
            "an instance turns about a literal axis_at(point, direction): an edge selection names an edge of a body, not of the assembly",
        ));
    };
    let turn = Placement::rotation(
        Point3::new(origin[0], origin[1], origin[2]),
        Vec3::new(direction[0], direction[1], direction[2]),
        args::float("angle", &degrees)?,
    )
    .map_err(|_| {
        argument_error(
            "axis",
            "the axis direction must be non-zero and the angle finite",
        )
    })?;
    i.moved(turn)
}

fn frame(
    i: &mut ScriptInstance,
    origin: Dynamic,
    z: Dynamic,
    x: Dynamic,
) -> Result<ScriptInstance, Box<EvalAltResult>> {
    let (o, z, x) = (
        args::vec3("origin", &origin)?,
        args::vec3("z", &z)?,
        args::vec3("x", &x)?,
    );
    let frame = Placement::frame(
        Point3::new(o[0], o[1], o[2]),
        Vec3::new(z[0], z[1], z[2]),
        Vec3::new(x[0], x[1], x[2]),
    )
    .map_err(|_| argument_error("frame", "z and x must be non-zero and not parallel"))?;
    i.moved(frame)
}

/// `instance.seat(axis_at(origin, direction))`: place a fastener on a
/// hole cut along that axis, in the hole's insertion frame (z down the
/// axis from its origin, x the fastener family's seed), so its thread
/// phase matches the hole's (PK-897).
fn seat(
    i: &mut ScriptInstance,
    axis: crate::profile_api::ScriptAxis,
) -> Result<ScriptInstance, Box<EvalAltResult>> {
    let crate::profile_api::ScriptAxis::Literal { origin, direction } = axis else {
        return Err(argument_error(
            "seat",
            "the hole axis must be axis_at(origin, direction)",
        ));
    };
    let z = Vec3::new(direction[0], direction[1], direction[2]);
    let frame = Placement::frame(
        Point3::new(origin[0], origin[1], origin[2]),
        z,
        popkern_fastener::holes::insertion_x_seed(z),
    )
    .map_err(|_| argument_error("seat", "the axis direction must be non-zero"))?;
    i.moved(frame)
}

fn assert_watertight(a: &mut ScriptAssembly) -> Result<(), Box<EvalAltResult>> {
    a.borrow_mut()
        .assembly
        .assert_watertight()
        .map_err(assembly_error)
}

/// `a.interference(#{ clearance, tolerance, among })`: every pair of
/// part occurrences within `clearance` (default 0: contact and overlap
/// only), as the report map (design §8.1); `among` (instance paths)
/// limits the check to those occurrences and the ones below them.
fn interference(a: &mut ScriptAssembly, options: Dynamic) -> Result<Dynamic, Box<EvalAltResult>> {
    let options = crate::args::Options::parse(
        "interference",
        &options,
        &["clearance", "tolerance", "among"],
    )?;
    let mut check = popkern_assembly::InterferenceOptions::default();
    if let Some(c) = options.with("clearance", args::float)? {
        check.clearance = c;
    }
    check.tolerance = options.with("tolerance", args::float)?;
    if let Some(paths) = options.with("among", |_, v: &Dynamic| {
        v.clone()
            .into_typed_array::<rhai::ImmutableString>()
            .map(|a| a.into_iter().map(|s| s.to_string()).collect::<Vec<_>>())
            .map_err(|_| argument_error("among", "expected an array of instance paths"))
    })? {
        check.among = paths;
    }
    let report = a.borrow_mut().assembly.interference(&check).assembly()?;
    rhai::serde::to_dynamic(&report).map_err(|e| argument_error("interference", e.to_string()))
}

/// Register `assembly()` / `assembly(name)` under `certification`.
pub fn register_constructors(engine: &mut Engine, certification: ContinuityCertification) {
    engine
        .register_fn("assembly", move || {
            ScriptAssembly::new("assembly", certification)
        })
        .register_fn("assembly", move |name: &str| {
            ScriptAssembly::new(name, certification)
        });
}

/// Register the assembly, definition and instance values.
pub fn register(engine: &mut Engine) {
    register_constructors(engine, ContinuityCertification::Inline);
    engine
        .register_type_with_name::<ScriptAssembly>("Assembly")
        .register_type_with_name::<ScriptDefinition>("Definition")
        .register_type_with_name::<ScriptInstance>("Instance")
        .register_fn("part", part)
        .register_fn("part", crate::assembly_family_api::family)
        .register_fn("variant", crate::assembly_family_api::variant)
        .register_fn("subassembly", subassembly)
        .register_fn("place", |a: &mut ScriptAssembly, d: ScriptDefinition| {
            place(a, d, None)
        })
        .register_fn(
            "place",
            |a: &mut ScriptAssembly, d: ScriptDefinition, name: &str| place(a, d, Some(name)),
        )
        .register_fn("translate", translate)
        .register_fn("rotate", rotate)
        .register_fn("frame", frame)
        .register_fn("seat", seat)
        .register_fn("hide", |i: &mut ScriptInstance| {
            let mut state = i.owner.borrow_mut();
            state.assembly.set_visible(i.at, false).assembly()
        })
        .register_get("name", |i: &mut ScriptInstance| i.name())
        .register_get("name", |d: &mut ScriptDefinition| {
            let state = d.owner.borrow_mut();
            Ok::<_, Box<EvalAltResult>>(state.assembly.definition(d.id).assembly()?.name.clone())
        })
        .register_fn("assert_watertight", assert_watertight)
        .register_fn("interference", |a: &mut ScriptAssembly| {
            interference(a, Dynamic::UNIT)
        })
        .register_fn("interference", interference)
        .register_fn("save", |a: &mut ScriptAssembly, path: &str| {
            a.borrow_mut().assembly.save(path).assembly()
        })
        .register_fn("to_step", |a: &mut ScriptAssembly, path: &str| {
            a.borrow_mut().assembly.to_step(path).assembly()
        });
}
