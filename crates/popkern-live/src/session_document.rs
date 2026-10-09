//! Describing an authored document of either kind for the page: the
//! meshes (per alive body; for an assembly, per alive body of each part
//! definition, shared by its instances), the parameter rows (qualified
//! `part:name` in an assembly), and the assembly view the tree panel and
//! the instanced renderer read.

use popkern_assembly::{Assembly, DefinitionContent, Document};
use popkern_script::Model;

use crate::envelope::{AssemblyView, BodyMesh, DefinitionView, InstanceView, ParamRow};
use crate::mesh_buffers::{MeshPurpose, body_mesh};

/// Meshes, the assembly view (assemblies only) and the first meshing
/// error of `doc`.
pub(crate) fn mesh_document(
    doc: &mut Document,
    purpose: MeshPurpose,
    cache: Option<&mut DisplayCache>,
) -> (Vec<BodyMesh>, Option<AssemblyView>, Option<String>) {
    match doc {
        Document::Part(model) => {
            let (bodies, error) = mesh_bodies(model, purpose, None);
            (bodies, None, error)
        }
        Document::Assembly(assembly) => {
            let (bodies, view, error) = mesh_assembly(assembly, purpose, cache);
            (bodies, Some(view), error)
        }
    }
}

/// An assembly's part meshes from the last describe, by definition
/// name, with the signature they were built under (PK-725): a joint
/// slider moves placements only, and re-meshing every part made one
/// slider step of the six-axis arm take 30 s. The signature (feature
/// count and parameter values) changes when a driven part parameter
/// does; a mutable borrow of the document clears the cache.
pub(crate) type DisplayCache = std::collections::HashMap<String, (u64, Vec<BodyMesh>)>;

fn signature(model: &Model) -> u64 {
    use std::hash::{Hash, Hasher};
    let mut h = std::hash::DefaultHasher::new();
    model.features().len().hash(&mut h);
    format!("{:?}", model.params()).hash(&mut h);
    h.finish()
}

fn mesh_assembly(
    assembly: &mut Assembly,
    purpose: MeshPurpose,
    mut cache: Option<&mut DisplayCache>,
) -> (Vec<BodyMesh>, AssemblyView, Option<String>) {
    let mut bodies = Vec::new();
    let mut error = None;
    let mut definitions = Vec::new();
    let count = assembly.definitions().len();
    for i in 0..count {
        let id = popkern_assembly::DefinitionId(i as u32);
        let Ok(definition) = assembly.definition(id) else {
            continue;
        };
        let name = definition.name.clone();
        let is_part = matches!(definition.content, DefinitionContent::Part(_));
        let mut indices = Vec::new();
        let mut features = 0;
        if is_part && let Ok(model) = assembly.part_model_mut(id) {
            features = model.features().len();
            let key = signature(model);
            let cached = cache
                .as_deref()
                .and_then(|c| c.get(&name))
                .filter(|(k, _)| *k == key)
                .map(|(_, meshes)| meshes.clone());
            let meshes = match cached {
                Some(meshes) => meshes,
                None => {
                    let (meshes, e) = mesh_bodies(model, purpose, Some(&name));
                    if let Some(e) = e {
                        error.get_or_insert(format!("in part '{name}': {e}"));
                    } else if let Some(c) = cache.as_deref_mut() {
                        c.insert(name.clone(), (key, meshes.clone()));
                    }
                    meshes
                }
            };
            for mesh in meshes {
                indices.push(bodies.len());
                bodies.push(mesh);
            }
        }
        definitions.push(DefinitionView {
            name,
            kind: if is_part { "part" } else { "assembly" },
            bodies: indices,
            features,
        });
    }
    let instances = assembly.tree().iter().map(instance_view).collect();
    let view = AssemblyView {
        name: assembly.name().to_string(),
        definitions,
        instances,
    };
    (bodies, view, error)
}

/// One tree row as the page's tree and instanced renderer see it.
pub(crate) fn instance_view(row: &popkern_assembly::TreeRow) -> InstanceView {
    InstanceView {
        name: row.path.0.last().cloned().unwrap_or_default(),
        depth: row.path.depth(),
        path: row.path.to_string(),
        definition: row.definition.0,
        world: row.world.column_major_f32(),
        visible: row.visible,
        is_part: row.is_part,
    }
}

/// Mesh every live body of `model` for `purpose`: the buffers, and the
/// first meshing error (a body that fails to mesh is skipped). In an
/// assembly the body label is the part's name.
pub(crate) fn mesh_bodies(
    model: &mut Model,
    purpose: MeshPurpose,
    part: Option<&str>,
) -> (Vec<BodyMesh>, Option<String>) {
    let mut bodies = Vec::new();
    let mut error = None;
    for body in model.alive_bodies() {
        let brep = match model.brep_of(body) {
            Ok(brep) => brep,
            Err(_) => continue,
        };
        let label = match part {
            Some(part) => Some(part.to_string()),
            None => model.body_label(body),
        };
        let label_face = |handle| {
            let name = brep.faces.get(handle)?.name?;
            let entity = brep.name_pool.get(name)?;
            model.spell_entity(entity).ok()
        };
        match body_mesh(&brep, body.id.0 as u64, label, label_face, purpose) {
            Ok(mut mesh) => {
                mesh.appearance = model.appearance();
                bodies.push(mesh);
            }
            Err(e) => {
                error.get_or_insert_with(|| e.to_string());
            }
        }
    }
    (bodies, error)
}

/// The parameter rows of `doc`: a part document's own, an assembly's
/// qualified by part (`planet:width`).
///
/// An assembly's own parameters come first (unqualified, PK-638); a part
/// parameter one of them drives is listed read-only with the driving
/// parameter as its expression (the assembly parameter is its single
/// source).
pub(crate) fn document_params(doc: &mut Document) -> Vec<ParamRow> {
    let mut rows = Vec::new();
    let mut driven: Vec<(String, String, String)> = Vec::new();
    if let Document::Assembly(assembly) = doc {
        rows = assembly_param_rows(assembly);
        driven = assembly
            .drives()
            .iter()
            .map(|d| (d.part.clone(), d.part_param.clone(), d.param.clone()))
            .collect();
    }
    for (part, model) in doc.models_mut() {
        for p in model.params() {
            let by = part.as_ref().and_then(|part| {
                driven
                    .iter()
                    .find(|(dp, pp, _)| dp == part && *pp == p.name)
                    .map(|(_, _, a)| a.clone())
            });
            rows.push(ParamRow {
                name: match &part {
                    Some(part) => format!("{part}:{}", p.name),
                    None => p.name,
                },
                value: p.value,
                unit: serde_json::to_value(p.unit).unwrap_or(serde_json::Value::Null),
                expression: by.or(p.expression),
                bounds: p.bounds.map(|(lo, hi)| [lo, hi]),
                read_only: p.read_only,
                explicit: p.explicit,
            });
        }
    }
    rows
}

/// An assembly's own parameter rows (unqualified, in declaration order).
pub(crate) fn assembly_param_rows(assembly: &Assembly) -> Vec<ParamRow> {
    assembly
        .params()
        .iter()
        .map(|p| ParamRow {
            name: p.name.clone(),
            value: p.display_value(),
            unit: serde_json::to_value(p.unit).unwrap_or(serde_json::Value::Null),
            expression: None,
            bounds: p.bounds,
            read_only: false,
            explicit: true,
        })
        .collect()
}
