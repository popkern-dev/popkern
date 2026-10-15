//! What a film script may name in the built document
//! (docs/architecture/presentation-films.md §3 invariant 4): occurrence
//! paths or bodies, the placement-only joints, and the unit scale that
//! turns document lengths into the SI lengths the page uses.

use popkern_assembly::Document;
use serde::{Deserialize, Serialize};

/// A read-only snapshot of a built document for film scripts.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct FilmContext {
    /// `"part"` or `"assembly"`.
    pub kind: String,
    /// The document's name (the assembly's name, or `"part"`).
    pub name: String,
    /// Metres per document length unit (a mm part: 0.001; assemblies 1).
    pub unit_scale: f64,
    /// Occurrence paths of an assembly, depth-first (empty for a part).
    pub paths: Vec<String>,
    /// Alive bodies of a part document (0 for an assembly).
    pub bodies: usize,
    /// Assembly parameters that drive no part: they only move
    /// occurrences, so a film may sweep them.
    pub joints: Vec<String>,
    /// Every other parameter (rebuilding ones; named in errors).
    pub params: Vec<String>,
    /// Current values of the joints and parameters, in their display
    /// units (`doc.value("j1")`: a film starts a sweep from the pose the
    /// document is in).
    #[serde(default)]
    pub values: std::collections::BTreeMap<String, f64>,
}

impl FilmContext {
    /// The snapshot of `document`.
    pub fn from_document(document: &mut Document) -> FilmContext {
        match document {
            Document::Part(model) => {
                let unit = model.length_unit();
                let scale = unit.to_si(1.0);
                let infos = model.params();
                let params = infos.iter().map(|p| p.name.clone()).collect();
                let values = infos
                    .iter()
                    .map(|p| (p.name.clone(), p.unit.from_si(p.value)))
                    .collect();
                FilmContext {
                    kind: "part".into(),
                    name: "part".into(),
                    unit_scale: if scale.is_finite() && scale > 0.0 {
                        scale
                    } else {
                        1.0
                    },
                    paths: Vec::new(),
                    bodies: model.alive_bodies().len(),
                    joints: Vec::new(),
                    params,
                    values,
                }
            }
            Document::Assembly(assembly) => {
                let paths = assembly
                    .tree()
                    .iter()
                    .map(|row| row.path.to_string())
                    .collect();
                let driving: Vec<&str> =
                    assembly.drives().iter().map(|d| d.param.as_str()).collect();
                let (joints, params) = assembly
                    .params()
                    .iter()
                    .map(|p| p.name.clone())
                    .partition(|name| !driving.contains(&name.as_str()));
                let values = assembly
                    .params()
                    .iter()
                    .map(|p| (p.name.clone(), p.display_value()))
                    .collect();
                FilmContext {
                    kind: "assembly".into(),
                    name: assembly.name().to_string(),
                    unit_scale: 1.0,
                    paths,
                    bodies: 0,
                    joints,
                    params,
                    values,
                }
            }
        }
    }

    /// Whether `target` selects something: `"all"`, `"body N"` of a
    /// part, or an occurrence path (or a path prefix: everything under
    /// it).
    pub fn selects(&self, target: &str) -> bool {
        if target == "all" {
            return true;
        }
        if let Some(n) = target.strip_prefix("body ") {
            return self.kind == "part" && n.parse::<usize>().is_ok_and(|n| n < self.bodies);
        }
        self.paths.iter().any(|p| {
            p == target
                || p.strip_prefix(target)
                    .is_some_and(|rest| rest.starts_with('/'))
        })
    }

    /// What may be named, for an error message.
    pub fn known_targets(&self) -> String {
        let mut names: Vec<String> = vec!["all".into()];
        if self.kind == "part" {
            names.extend((0..self.bodies).map(|i| format!("body {i}")));
        } else {
            names.extend(self.paths.iter().take(24).cloned());
            if self.paths.len() > 24 {
                names.push(format!("… {} more", self.paths.len() - 24));
            }
        }
        names.join(", ")
    }
}
