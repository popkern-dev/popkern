//! "Show as code" for an assembly (design §6.1): the recorded documents
//! are the truth, so every part definition is emitted as a
//! `a.part(name, |m| { … })` holding its own canonical Rhai (a family
//! variant becomes its own part — the closure that made it is not
//! recorded), sub-assemblies as `a.subassembly(name, |s| { … })` after
//! the definitions they place, and every root instance as a `place`
//! whose placement is spelled as nothing (identity), `.translate(…)`
//! (no rotation) or `.frame(origin, z, x)`. Running the result rebuilds
//! the same tree with every placement to 1e-12.

use popkern_assembly::{Assembly, DefinitionContent, DefinitionId, Document, Instance, Placement};

use crate::emit::{EmitError, to_rhai};
use crate::emit_spelling::{float, quoted};

/// The canonical Rhai of a document of either kind.
///
/// # Errors
///
/// [`EmitError`] when a part has no spelling (as [`to_rhai`]), or a
/// part's build tolerance differs from the default (a part closure
/// records into a default-tolerance model).
pub fn document_to_rhai(document: &mut Document) -> Result<String, EmitError> {
    match document {
        Document::Part(model) => to_rhai(model),
        Document::Assembly(assembly) => assembly_to_rhai(assembly),
    }
}

/// The canonical Rhai of an assembly.
///
/// # Errors
///
/// As [`document_to_rhai`].
pub fn assembly_to_rhai(assembly: &mut Assembly) -> Result<String, EmitError> {
    let mut out = format!("let a = assembly({});\n", quoted(assembly.name()));
    // A joint's own parameter comes back with its `revolute` verb.
    let joint_params: Vec<String> = assembly
        .joints()
        .iter()
        .filter(|j| j.drive == popkern_assembly::JointDrive::Param)
        .map(|j| j.name.clone())
        .collect();
    for param in assembly.params() {
        if joint_params.contains(&param.name) {
            continue;
        }
        out.push_str(&format!(
            "a.param({}, {}, {});\n",
            quoted(&param.name),
            float(param.display_value()),
            quoted(crate::args::unit_name(param.unit))
        ));
        if let Some([lo, hi]) = param.bounds {
            out.push_str(&format!(
                "a.limits({}, {}, {});\n",
                quoted(&param.name),
                float(lo),
                float(hi)
            ));
        }
    }
    let count = assembly.definitions().len();
    // Parts first (definition order), then sub-assemblies children
    // before parents, so every handle exists before it is placed.
    for i in 0..count {
        let id = DefinitionId(i as u32);
        let Ok(model) = assembly.part_model_mut(id) else {
            continue;
        };
        let default = popkern_core::tolerance::Tolerance::default().linear();
        if model.tolerance().linear() != default {
            return Err(EmitError::NotApplicable {
                feature: part_name(assembly, id),
                operation: "part".into(),
                detail: "a part with a non-default build tolerance has no assembly spelling: \
                         a part closure records into a default-tolerance model"
                    .into(),
            });
        }
        let code = to_rhai(model)?;
        let name = part_name(assembly, id);
        out.push_str(&format!("let d{i} = a.part({}, |m| {{\n", quoted(&name)));
        for line in body_lines(&code) {
            out.push_str(&format!("    {line}\n"));
        }
        out.push_str("});\n");
    }
    for drive in assembly.drives().to_vec() {
        if let Some(id) = assembly.definition_named(&drive.part) {
            out.push_str(&format!(
                "a.drive({}, d{}, {});\n",
                quoted(&drive.param),
                id.0,
                quoted(&drive.part_param)
            ));
        }
    }
    let mut emitted: Vec<DefinitionId> = (0..count)
        .map(|i| DefinitionId(i as u32))
        .filter(|id| assembly.part_model_mut(*id).is_ok())
        .collect();
    while emitted.len() < count {
        let before = emitted.len();
        for i in 0..count {
            let id = DefinitionId(i as u32);
            if emitted.contains(&id) {
                continue;
            }
            let Ok(definition) = assembly.definition(id) else {
                continue;
            };
            let DefinitionContent::Assembly(instances) = &definition.content else {
                continue;
            };
            if !instances
                .iter()
                .all(|inst| emitted.contains(&inst.definition))
            {
                continue;
            }
            out.push_str(&format!(
                "let d{i} = a.subassembly({}, |s| {{\n",
                quoted(&definition.name)
            ));
            for instance in instances {
                out.push_str(&format!("    {};\n", place("s", instance)));
            }
            out.push_str("});\n");
            emitted.push(id);
        }
        if emitted.len() == before {
            break;
        }
    }
    let root: Vec<Instance> = assembly
        .instances(popkern_assembly::Parent::Root)
        .map(<[Instance]>::to_vec)
        .unwrap_or_default();
    for instance in &root {
        out.push_str(&format!("{};\n", place("a", instance)));
    }
    // The placements above already satisfy the mates, so declaring them
    // moves nothing; they come back for the solver and the sliders
    // (PK-711).
    if let Some(ground) = assembly.ground() {
        out.push_str(&format!("a.ground(a.occurrence({}));\n", quoted(ground)));
    }
    for (index, mate) in assembly.mates().iter().enumerate() {
        match assembly.joint_of_mate(index) {
            Some(joint) if joint.mates[0] == index => {
                out.push_str(&joint_statement(assembly, joint));
            }
            Some(_) => {}
            None => out.push_str(&mate_statement(mate)),
        }
    }
    out.push_str("a\n");
    Ok(out)
}

fn part_name(assembly: &Assembly, id: DefinitionId) -> String {
    assembly
        .definition(id)
        .map(|d| d.name.clone())
        .unwrap_or_default()
}

/// A part's canonical script minus its `let m = model(…);` header and
/// trailing `m`: the body of the part closure.
fn body_lines(code: &str) -> Vec<&str> {
    let lines: Vec<&str> = code.lines().collect();
    let start = usize::from(
        lines
            .first()
            .is_some_and(|l| l.starts_with("let m = model(")),
    );
    let end = if lines.last().is_some_and(|l| l.trim() == "m") {
        lines.len() - 1
    } else {
        lines.len()
    };
    lines[start..end.max(start)].to_vec()
}

fn place(owner: &str, instance: &Instance) -> String {
    let mut text = format!(
        "{owner}.place(d{}, {}){}",
        instance.definition.0,
        quoted(&instance.name),
        placement(&instance.placement)
    );
    if !instance.visible {
        text.push_str(".hide()");
    }
    text
}

fn placement(p: &Placement) -> String {
    if p.is_identity() {
        return String::new();
    }
    let o = p.origin();
    if p.is_translation() {
        return format!(".translate({}, {}, {})", float(o.x), float(o.y), float(o.z));
    }
    let (z, x) = (p.z_axis(), p.x_axis());
    format!(
        ".frame([{}, {}, {}], [{}, {}, {}], [{}, {}, {}])",
        float(o.x),
        float(o.y),
        float(o.z),
        float(z.x),
        float(z.y),
        float(z.z),
        float(x.x),
        float(x.y),
        float(x.z)
    )
}

/// `a.occurrence("path").datum("z")` / `.face("name")`.
fn mate_ref(r: &popkern_assembly::MateRef) -> String {
    let entity = match &r.entity {
        popkern_assembly::MateEntity::Datum { datum } => {
            let d = match datum {
                popkern_assembly::Datum::Xy => "xy",
                popkern_assembly::Datum::Yz => "yz",
                popkern_assembly::Datum::Zx => "zx",
                popkern_assembly::Datum::X => "x",
                popkern_assembly::Datum::Y => "y",
                popkern_assembly::Datum::Z => "z",
            };
            format!("datum({})", quoted(d))
        }
        popkern_assembly::MateEntity::Face { name } => format!("face({})", quoted(name)),
        popkern_assembly::MateEntity::Axis { origin, direction } => {
            format!("axis({}, {})", triple(origin), triple(direction))
        }
        popkern_assembly::MateEntity::Plane { origin, normal } => {
            format!("plane({}, {})", triple(origin), triple(normal))
        }
    };
    format!("a.occurrence({}).{entity}", quoted(&r.instance))
}

/// A joint as its `revolute` verb, from its three mates: the axes of
/// the concentric, the axial offset, the angle mate's zero entities and
/// the driving parameter's value and limits (or the follow ratio and
/// phase).
fn joint_statement(assembly: &Assembly, joint: &popkern_assembly::Joint) -> String {
    use popkern_assembly::{JointDrive, MateKind};
    let mates = assembly.mates();
    let (Some(axes), Some(seat), Some(angle)) = (
        mates.get(joint.mates[0]),
        mates.get(joint.mates[1]),
        mates.get(joint.mates[2]),
    ) else {
        return String::new();
    };
    let mut opts = vec![format!(
        "zero: [{}, {}]",
        mate_ref(&angle.a),
        mate_ref(&angle.b)
    )];
    if let MateKind::AxialOffset { distance } = seat.kind
        && distance != 0.0
    {
        opts.push(format!("offset: {}", float(distance)));
    }
    match &joint.drive {
        JointDrive::Param => {
            if let Some(param) = assembly.params().iter().find(|p| p.name == joint.name) {
                opts.push(format!("angle: {}", float(param.display_value())));
                if let Some([lo, hi]) = param.bounds {
                    opts.push(format!("limits: [{}, {}]", float(lo), float(hi)));
                }
            }
        }
        JointDrive::Follow {
            leader,
            ratio,
            phase,
        } => {
            opts.push(format!("follow: {}", quoted(leader)));
            opts.push(format!("ratio: {}", float(*ratio)));
            opts.push(format!("phase: {}", float(*phase)));
        }
    }
    format!(
        "a.revolute({}, {}, {}, #{{ {} }});\n",
        quoted(&joint.name),
        mate_ref(&axes.a),
        mate_ref(&axes.b),
        opts.join(", ")
    )
}

/// One mate as its verb; a driven value as `p("param")` with its scale
/// and offset when they are not the identity.
fn mate_statement(mate: &popkern_assembly::Mate) -> String {
    use popkern_assembly::MateKind;
    let (a, b) = (mate_ref(&mate.a), mate_ref(&mate.b));
    let value = |v: f64| match &mate.drive {
        Some(d) => format!("p({})", quoted(&d.param)),
        None => float(v),
    };
    let mut pairs: Vec<String> = Vec::new();
    if let Some(d) = &mate.drive {
        if d.scale != 1.0 {
            pairs.push(format!("scale: {}", float(d.scale)));
        }
        if d.offset != 0.0 {
            pairs.push(format!("offset: {}", float(d.offset)));
        }
    }
    let opts = |mut extra: Vec<String>, pairs: &[String]| {
        extra.extend(pairs.iter().cloned());
        if extra.is_empty() {
            String::new()
        } else {
            format!(", #{{ {} }}", extra.join(", "))
        }
    };
    match &mate.kind {
        MateKind::Concentric => format!("a.concentric({a}, {b});\n"),
        MateKind::Coincident { aligned } => {
            let extra = if *aligned {
                vec!["aligned: true".to_string()]
            } else {
                vec![]
            };
            format!("a.coincident({a}, {b}{});\n", opts(extra, &[]))
        }
        MateKind::Distance { distance, aligned } => {
            let extra = if *aligned {
                vec!["aligned: true".to_string()]
            } else {
                vec![]
            };
            format!(
                "a.distance({a}, {b}, {}{});\n",
                value(*distance),
                opts(extra, &pairs)
            )
        }
        // Only joints make axial offsets, and a joint is emitted whole
        // (`joint_statement`).
        MateKind::AxialOffset { .. } => String::new(),
        MateKind::Angle { degrees, about } => {
            let extra = about
                .as_ref()
                .map(|r| vec![format!("about: {}", mate_ref(r))])
                .unwrap_or_default();
            format!(
                "a.angle({a}, {b}, {}{});\n",
                value(*degrees),
                opts(extra, &pairs)
            )
        }
    }
}

fn triple(v: &[f64; 3]) -> String {
    format!("[{}, {}, {}]", float(v[0]), float(v[1]), float(v[2]))
}
