//! The recorder session (spec §2): owns the feature tree, the
//! document tolerance, and the active body.

use popkern_core::naming::EntityName;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_history::FeatureTree;
use popkern_history::node::{BodyId, FeatureId};
use popkern_history::replay_report::{FeatureReport, ReplayReport};
use popkern_history::unit::Unit;
use popkern_part::Part;

use crate::body::Body;
use crate::error::ScriptError;
use crate::profile::ProfileBuilder;
use crate::selection::EntityRef;

/// A recorded feature: its id and the label it received.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FeatureHandle {
    /// Tree id (never changes).
    pub id: FeatureId,
    /// Document-unique label (auto or user-supplied).
    pub label: String,
    /// Operation identifier.
    pub operation: String,
}

/// A recording session over one document.
pub struct Model {
    pub(crate) tree: FeatureTree,
    pub(crate) tolerance: Tolerance,
    pub(crate) unit: Unit,
    pub(crate) active: Option<BodyId>,
    pub(crate) pending_counter: u64,
}

impl Model {
    /// Begin a new document.
    pub fn new() -> Model {
        Model::with_tolerance(Tolerance::default())
    }

    /// Begin a new document under an explicit build tolerance.
    pub fn with_tolerance(tolerance: Tolerance) -> Model {
        let mut tree = FeatureTree::new();
        tree.set_build_tolerance(Some(tolerance));
        Model {
            tree,
            tolerance,
            unit: Unit::Meters,
            active: None,
            pending_counter: 0,
        }
    }

    /// Set the length unit new literals are recorded in.
    pub fn unit(mut self, unit: Unit) -> Model {
        self.unit = unit;
        self
    }

    /// Open a `.pkn` document. No replay runs until requested.
    pub fn open(path: impl AsRef<std::path::Path>) -> Result<Model, ScriptError> {
        let path = path.as_ref();
        let data = std::fs::read(path).map_err(|e| ScriptError::Io {
            path: path.display().to_string(),
            detail: e.to_string(),
        })?;
        Model::from_pkn_bytes(&data)
    }

    /// Save the document at the current format version.
    pub fn save(&self, path: impl AsRef<std::path::Path>) -> Result<(), ScriptError> {
        let path = path.as_ref();
        let data = self.to_pkn_bytes()?;
        std::fs::write(path, data).map_err(|e| ScriptError::Io {
            path: path.display().to_string(),
            detail: e.to_string(),
        })
    }

    /// The document as `.pkn` bytes at the current format version —
    /// exactly what [`Model::save`] writes. Hosts without a filesystem
    /// (the browser playground) hand these to a download.
    pub fn to_pkn_bytes(&self) -> Result<Vec<u8>, ScriptError> {
        Ok(self.tree.save()?)
    }

    /// Open a document from `.pkn` bytes (the in-memory form of
    /// [`Model::open`]).
    pub fn from_pkn_bytes(data: &[u8]) -> Result<Model, ScriptError> {
        let tree = FeatureTree::load(data)?;
        let tolerance = tree.effective_build_tolerance();
        let active = tree
            .timeline()
            .last()
            .and_then(|&id| tree.feature(id))
            .map(|n| n.body);
        Ok(Model {
            tree,
            tolerance,
            unit: Unit::Meters,
            active,
            pending_counter: 0,
        })
    }

    /// The document's build tolerance.
    pub fn tolerance(&self) -> Tolerance {
        self.tolerance
    }

    /// The underlying feature tree (read access for migrations and tests).
    pub fn tree(&self) -> &FeatureTree {
        &self.tree
    }

    /// Mutable access to the tree for migration verbs the surface does
    /// not wrap yet. The caller owns replay afterwards.
    pub fn tree_mut(&mut self) -> &mut FeatureTree {
        &mut self.tree
    }

    /// Full replay under the document tolerance.
    pub fn replay(&mut self) -> Result<ReplayReport, ScriptError> {
        let result = self.tree.replay(self.tolerance)?;
        Ok(result.report)
    }

    /// The last replay report.
    pub fn report(&self) -> &ReplayReport {
        self.tree.last_report()
    }

    /// A feature's report row from the last replay.
    pub fn evidence(&self, feature: &FeatureHandle) -> Option<&FeatureReport> {
        self.tree.last_report().feature(feature.id)
    }

    /// Start an analytic profile on `plane`.
    pub fn profile(&mut self, plane: crate::plane::Plane) -> ProfileBuilder<'_> {
        ProfileBuilder::new(self, plane)
    }

    /// Record a tool model: everything the closure records starts a
    /// new body (sketch on a base plane), and the body is returned for
    /// a boolean.
    pub fn tool(
        &mut self,
        build: impl FnOnce(&mut Model) -> Result<Body, ScriptError>,
    ) -> Result<Body, ScriptError> {
        let previous = self.active;
        let body = build(self)?;
        self.active = previous.or(Some(body.id));
        Ok(body)
    }

    /// The active body (the one new features attach to).
    pub fn active(&self) -> Result<Body, ScriptError> {
        self.active.map(Body::new).ok_or(ScriptError::NoActiveBody)
    }

    /// Bodies alive at the end of the timeline.
    pub fn bodies(&self) -> Vec<Body> {
        self.tree.bodies().into_iter().map(Body::new).collect()
    }

    /// Look a feature up by label.
    pub fn feature(&self, label: &str) -> Result<FeatureHandle, ScriptError> {
        let id = self
            .tree
            .feature_by_label(label)
            .ok_or_else(|| popkern_history::HistoryError::LabelNotFound(label.to_string()))?;
        self.handle(id)
    }

    /// All features in timeline order.
    pub fn features(&self) -> Vec<FeatureHandle> {
        self.tree
            .timeline()
            .iter()
            .filter_map(|&id| self.handle(id).ok())
            .collect()
    }

    pub(crate) fn handle(&self, id: FeatureId) -> Result<FeatureHandle, ScriptError> {
        let node = self
            .tree
            .feature(id)
            .ok_or(popkern_history::HistoryError::FeatureNotFound)?;
        Ok(FeatureHandle {
            id,
            label: node.label.clone(),
            operation: node.feature.operation_label().to_string(),
        })
    }

    /// Set a parameter by qualified name (`extrude-1.distance`, `wall`).
    pub fn set_param(&mut self, name: &str, value: f64) -> Result<(), ScriptError> {
        let id =
            self.tree
                .parameters()
                .resolve_name(name)
                .ok_or_else(|| ScriptError::Argument {
                    argument: "parameter".into(),
                    reason: format!("unknown parameter '{name}'"),
                })?;
        let unit = self
            .tree
            .parameters()
            .get(id)
            .map(|p| p.unit)
            .unwrap_or(self.unit);
        self.tree.set_parameter_value_display(id, value, unit)?;
        Ok(())
    }

    /// Declare a global parameter (`value` in `unit`; stored in SI).
    pub fn param(&mut self, name: &str, value: f64, unit: Unit) -> Result<(), ScriptError> {
        self.tree.add_global_parameter_display(name, value, unit)?;
        Ok(())
    }

    /// Rename a feature (strict: `DuplicateLabel` / `EmptyLabel`).
    pub fn rename(&mut self, feature: &FeatureHandle, label: &str) -> Result<(), ScriptError> {
        Ok(self.tree.rename_feature(feature.id, label)?)
    }

    /// Bind an entity label (`DuplicateLabel` if taken).
    pub fn label(&mut self, name: &str, entity: EntityName) -> Result<(), ScriptError> {
        Ok(self.tree.label_entity(name, entity)?)
    }

    /// Bind an entity label to the single edge a selection resolves
    /// to on the active body (spec §6.3 `label(name, selection)`).
    pub fn label_edge(
        &mut self,
        name: &str,
        selection: &crate::selection::EdgeSelector,
    ) -> Result<EntityName, ScriptError> {
        let body = self.active()?;
        let brep = self.body_brep(body.id)?;
        let labels = |r: EntityRef| self.resolve_entity_ref(r);
        let snapshot = crate::selection::resolve_edges(&brep, selection, &labels)?;
        let entity = crate::selection::exactly_one(snapshot, crate::error::SelectionKind::Edges)?;
        self.tree.label_entity(name, entity.clone())?;
        Ok(entity)
    }

    /// Bind an entity label to the single face a selection resolves to.
    pub fn label_face(
        &mut self,
        name: &str,
        selection: &crate::selection::FaceSelector,
    ) -> Result<EntityName, ScriptError> {
        let body = self.active()?;
        let brep = self.body_brep(body.id)?;
        let labels = |r: EntityRef| self.resolve_entity_ref(r);
        let snapshot = crate::selection::resolve_faces(&brep, selection, &labels)?;
        let entity = crate::selection::exactly_one(snapshot, crate::error::SelectionKind::Faces)?;
        self.tree.label_entity(name, entity.clone())?;
        Ok(entity)
    }

    /// Remove an entity label.
    pub fn unlabel(&mut self, name: &str) -> Result<EntityName, ScriptError> {
        Ok(self.tree.unlabel_entity(name)?)
    }

    /// Rename an entity label (references never change).
    pub fn rename_label(&mut self, old: &str, new: &str) -> Result<(), ScriptError> {
        Ok(self.tree.rename_entity_label(old, new)?)
    }

    /// All entity labels.
    pub fn labels(&self) -> Vec<(String, EntityName)> {
        self.tree
            .labels()
            .entities
            .iter()
            .map(|(l, e)| (l.clone(), e.clone()))
            .collect()
    }

    /// Resolve an entity label.
    pub fn entity_by_label(&self, label: &str) -> Result<EntityName, ScriptError> {
        self.tree
            .entity_by_label(label)
            .cloned()
            .ok_or_else(|| popkern_history::HistoryError::LabelNotFound(label.to_string()).into())
    }

    /// The current BRep of the active body.
    pub fn to_brep(&mut self) -> Result<BRep, ScriptError> {
        let body = self.active()?;
        self.body_brep(body.id)
    }

    pub(crate) fn body_brep(&mut self, body: BodyId) -> Result<BRep, ScriptError> {
        if self.tree.body_brep(body).is_none() {
            self.tree.replay(self.tolerance)?;
        }
        self.tree
            .body_brep(body)
            .cloned()
            .ok_or(ScriptError::BodyNotAlive(body.0 as u64))
    }

    /// The active body as a validated `Part`.
    pub fn to_part(&mut self) -> Result<Part, ScriptError> {
        let brep = self.to_brep()?;
        Ok(Part::from_brep(brep, self.tolerance, "popkern-script")?)
    }

    /// Export the active body as STEP.
    pub fn to_step(&mut self, path: impl AsRef<std::path::Path>) -> Result<(), ScriptError> {
        Ok(self.to_part()?.export_step(path)?)
    }

    /// The active body as STEP Part 21 bytes with `name` as the
    /// product name — what [`Model::to_step`] writes for a path whose
    /// file stem is `name`. Hosts without a filesystem (the browser
    /// playground) hand these to a download.
    pub fn to_step_bytes(&mut self, name: &str) -> Result<Vec<u8>, ScriptError> {
        Ok(self.to_part()?.step_bytes(name)?)
    }

    /// The current BRep of `body` (replaying first when the tree
    /// holds no snapshot for it).
    pub fn brep_of(&mut self, body: Body) -> Result<BRep, ScriptError> {
        self.body_brep(body.id)
    }

    /// Gate: the active body tessellates watertight.
    pub fn assert_watertight(&mut self) -> Result<(), ScriptError> {
        let open_edges = self.to_part()?.open_edge_count()?;
        if open_edges > 0 {
            return Err(ScriptError::NotWatertight { open_edges });
        }
        Ok(())
    }

    /// Gate: STEP export re-imports with the same face count.
    pub fn assert_step_roundtrip(
        &mut self,
        path: impl AsRef<std::path::Path>,
    ) -> Result<(), ScriptError> {
        let part = self.to_part()?;
        part.export_step(&path)?;
        part.assert_step_roundtrip(&path)?;
        Ok(())
    }
}

impl Default for Model {
    fn default() -> Self {
        Model::new()
    }
}
