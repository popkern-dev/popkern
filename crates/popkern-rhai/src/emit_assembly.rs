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
    for param in assembly.params() {
        out.push_str(&format!(
            "a.param({}, {}, {});\n",
            quoted(&param.name),
            float(param.display_value()),
            quoted(crate::args::unit_name(param.unit))
        ));
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
