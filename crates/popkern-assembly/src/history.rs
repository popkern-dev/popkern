//! The assembly's construction history (live-projects-and-history.md
//! §7, PK-666): one step per authoring operation that changes what is
//! placed or where — placing an instance, moving one, adding a mate,
//! an explicit solve — each recording the placements it set.
//!
//! Mates are solved incrementally from the current placements, so a
//! later replay from scratch need not land where authoring did: the
//! step states are therefore *recorded*, not recomputed. The state
//! after step `k` is the steps' placement deltas applied in order, and
//! the last state equals the document's placements by construction
//! ([`Assembly::check_history`] asserts it).
//!
//! A placement set right after its own place or move (`place(…)
//! .translate(…).rotate(…)`) folds into that step, so a step is one
//! authored intent, not one call. An assembly opened from `.pkn` or
//! STEP has no authoring record: its log is the stored instances in
//! order at their stored (solved) placements, and
//! [`AssemblyHistory::recorded`] is false so viewers can say so.

use std::collections::HashMap;

use popkern_script::source::SourceSpan;

use crate::assembly::{Assembly, InstanceRef, Parent};
use crate::placement::Placement;

/// What one step did.
#[derive(Debug, Clone, PartialEq)]
pub enum AssemblyStepKind {
    /// An instance was placed.
    Place,
    /// An instance was moved.
    Move,
    /// A mate was added (the index into [`Assembly::mates`]) and solved.
    Mate(usize),
    /// The mates were re-solved.
    Solve,
}

/// One construction step.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyStep {
    /// What it did.
    pub kind: AssemblyStepKind,
    /// Words for a timeline (`place pin_0 (pin)`, the mate described).
    pub label: String,
    /// The placements the step set, in order: the placed or moved
    /// instance, or every root instance a solve moved.
    pub placements: Vec<(InstanceRef, Placement)>,
    /// Where in a script it was authored, when tracked.
    pub source: Option<SourceSpan>,
}

/// The construction steps of an assembly.
#[derive(Debug, Clone, PartialEq)]
pub struct AssemblyHistory {
    /// The steps in authoring order.
    pub steps: Vec<AssemblyStep>,
    /// Whether the steps were recorded while authoring (false for an
    /// assembly opened from a file: placements only, at their stored
    /// positions, and no mate steps).
    pub recorded: bool,
}

impl Default for AssemblyHistory {
    fn default() -> Self {
        AssemblyHistory {
            steps: Vec::new(),
            recorded: true,
        }
    }
}

/// Every instance present after a step, with its placement.
pub type StepPlacements = Vec<(InstanceRef, Placement)>;

impl Assembly {
    /// The construction history.
    pub fn history(&self) -> &AssemblyHistory {
        &self.history
    }

    /// Mark the history as not recorded (an assembly opened from a
    /// file).
    pub(crate) fn mark_history_unrecorded(&mut self) {
        self.history.recorded = false;
    }

    pub(crate) fn record_place(&mut self, at: InstanceRef, placement: Placement) {
        let label = match self.instance(at) {
            Ok(instance) => {
                let definition = self
                    .definition(instance.definition)
                    .map(|d| d.name.clone())
                    .unwrap_or_default();
                format!(
                    "place {} ({definition})",
                    self.step_path(at, &instance.name)
                )
            }
            Err(_) => "place".to_string(),
        };
        self.push_step(AssemblyStepKind::Place, label, vec![(at, placement)]);
    }

    pub(crate) fn record_move(&mut self, at: InstanceRef, placement: Placement) {
        if let Some(last) = self.history.steps.last_mut()
            && matches!(last.kind, AssemblyStepKind::Place | AssemblyStepKind::Move)
            && last.placements.len() == 1
            && last.placements[0].0 == at
        {
            last.placements[0].1 = placement;
            return;
        }
        let name = self
            .instance(at)
            .map(|i| i.name.clone())
            .unwrap_or_default();
        let label = format!("move {}", self.step_path(at, &name));
        self.push_step(AssemblyStepKind::Move, label, vec![(at, placement)]);
    }

    /// Record a solve that started from root placements `before`.
    pub(crate) fn record_solve(&mut self, kind: AssemblyStepKind, before: &[Placement]) {
        let label = match &kind {
            AssemblyStepKind::Mate(index) => self
                .mates
                .get(*index)
                .map(|m| format!("mate: {}", m.describe()))
                .unwrap_or_else(|| "mate".to_string()),
            _ => "solve mates".to_string(),
        };
        let moved = self
            .root
            .iter()
            .enumerate()
            .filter(|(i, instance)| before.get(*i) != Some(&instance.placement))
            .map(|(index, instance)| {
                (
                    InstanceRef {
                        parent: Parent::Root,
                        index,
                    },
                    instance.placement,
                )
            })
            .collect();
        self.push_step(kind, label, moved);
    }

    fn push_step(&mut self, kind: AssemblyStepKind, label: String, placements: StepPlacements) {
        // Consecutive re-solves (a parameter dragged, `solve()` called
        // twice) are one step: merge their moves.
        if kind == AssemblyStepKind::Solve {
            if placements.is_empty() {
                return;
            }
            if let Some(last) = self.history.steps.last_mut()
                && last.kind == AssemblyStepKind::Solve
            {
                for (at, placement) in placements {
                    match last.placements.iter_mut().find(|(r, _)| *r == at) {
                        Some(slot) => slot.1 = placement,
                        None => last.placements.push((at, placement)),
                    }
                }
                return;
            }
        }
        self.history.steps.push(AssemblyStep {
            kind,
            label,
            placements,
            source: popkern_script::source::current(),
        });
    }

    fn step_path(&self, at: InstanceRef, name: &str) -> String {
        match at.parent {
            Parent::Root => name.to_string(),
            Parent::Definition(_) => format!("{}/{name}", self.parent_name(at.parent)),
        }
    }

    /// The instances present after the first `count` steps, with their
    /// placements then (`count` past the end is the final state).
    pub fn placements_after(&self, count: usize) -> StepPlacements {
        let mut present: StepPlacements = Vec::new();
        let mut slot: HashMap<InstanceRef, usize> = HashMap::new();
        for step in self.history.steps.iter().take(count) {
            for (at, placement) in &step.placements {
                match slot.get(at) {
                    Some(&i) => present[i].1 = *placement,
                    None => {
                        slot.insert(*at, present.len());
                        present.push((*at, *placement));
                    }
                }
            }
        }
        present
    }

    /// Whether replaying the whole history reproduces every instance's
    /// placement (the invariant the step states rest on). Returns the
    /// first instance that disagrees, by parent and index.
    pub fn check_history(&self) -> Result<(), InstanceRef> {
        let last = self.placements_after(self.history.steps.len());
        let mut parents = vec![Parent::Root];
        parents.extend(
            (0..self.definitions.len()).map(|i| Parent::Definition(crate::DefinitionId(i as u32))),
        );
        for parent in parents {
            let Ok(instances) = self.instances(parent) else {
                continue;
            };
            for (index, instance) in instances.iter().enumerate() {
                let at = InstanceRef { parent, index };
                let replayed = last.iter().find(|(r, _)| *r == at).map(|(_, p)| *p);
                if replayed != Some(instance.placement) {
                    return Err(at);
                }
            }
        }
        Ok(())
    }
}
