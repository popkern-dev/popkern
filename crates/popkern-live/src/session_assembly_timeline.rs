//! The construction timeline of an assembly document
//! (live-projects-and-history.md §7, PK-666): one row per authored
//! place, move, mate or solve, and the occurrences after any step —
//! placements only, because the page already holds every definition's
//! meshes from the run.

use serde::Serialize;

use popkern_assembly::history::AssemblyStepKind;
use popkern_assembly::{Assembly, Document, InstanceRef, Parent, TreeRow};
use popkern_rhai::Outcome;
use popkern_script::source::SourceSpan;

use crate::envelope::InstanceView;
use crate::session::Session;

/// One construction step.
#[derive(Debug, Clone, Serialize)]
pub struct AssemblyTimelineStep {
    /// 1-based: the assembly after this many steps.
    pub index: usize,
    /// `place`, `move`, `mate` or `solve`.
    pub kind: &'static str,
    /// Words for the step (`place pin_0 (pin)`, the mate described).
    pub label: String,
    /// The occurrence paths the step placed or moved: an instance
    /// inside a sub-assembly definition is listed once per occurrence of
    /// that definition (`s1/a`, `s2/a`).
    pub instances: Vec<String>,
    /// The project migration that authored it, in a project.
    pub migration: Option<String>,
    /// Where in a script it was authored, when tracked.
    pub source: Option<SourceSpan>,
}

/// The construction timeline of an assembly.
#[derive(Debug, Clone, Serialize)]
pub struct AssemblyTimeline {
    /// Steps in authoring order.
    pub steps: Vec<AssemblyTimelineStep>,
    /// False for an assembly opened from a file: its steps are the
    /// stored instances at their stored placements, with no mate steps.
    pub recorded: bool,
}

/// Every occurrence path of instance `at` in the expanded `tree`.
fn occurrence_paths(assembly: &Assembly, tree: &[TreeRow], at: InstanceRef) -> Vec<String> {
    let Ok(instance) = assembly.instance(at) else {
        return Vec::new();
    };
    match at.parent {
        Parent::Root => vec![instance.name.clone()],
        Parent::Definition(container) => tree
            .iter()
            .filter(|row| row.definition == container)
            .map(|row| format!("{}/{}", row.path, instance.name))
            .collect(),
    }
}

impl Session {
    /// The assembly's construction timeline; `None` for a part
    /// document.
    pub fn assembly_timeline(&self) -> Result<Option<AssemblyTimeline>, Outcome> {
        let Document::Assembly(assembly) = self
            .document
            .as_ref()
            .ok_or_else(crate::session::no_model)?
        else {
            return Ok(None);
        };
        let history = assembly.history();
        let tree = assembly.tree();
        let steps = history
            .steps
            .iter()
            .enumerate()
            .map(|(i, step)| AssemblyTimelineStep {
                index: i + 1,
                kind: match step.kind {
                    AssemblyStepKind::Place => "place",
                    AssemblyStepKind::Move => "move",
                    AssemblyStepKind::Mate(_) => "mate",
                    AssemblyStepKind::Solve => "solve",
                },
                label: step.label.clone(),
                instances: step
                    .placements
                    .iter()
                    .flat_map(|(at, _)| occurrence_paths(assembly, &tree, *at))
                    .collect(),
                migration: self
                    .assembly_migrations
                    .iter()
                    .find(|(_, after)| i < *after)
                    .map(|(name, _)| name.clone()),
                source: step.source.clone(),
            })
            .collect();
        Ok(Some(AssemblyTimeline {
            steps,
            recorded: history.recorded,
        }))
    }

    /// The occurrences after the first `index` construction steps, in
    /// the run's tree order and definition indices (an occurrence not
    /// yet placed is absent).
    pub fn assembly_step(&self, index: usize) -> Result<Vec<InstanceView>, Outcome> {
        let Some(Document::Assembly(assembly)) = self.document.as_ref() else {
            return Err(Outcome::GateFailure {
                step: "assembly_step".into(),
                detail: "construction steps belong to an assembly; this document is one part"
                    .into(),
            });
        };
        Ok(assembly
            .tree_after(index)
            .iter()
            .map(crate::session_document::instance_view)
            .collect())
    }
}
