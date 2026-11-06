//! A body as it is read back (PK-1102): every serialized field of
//! [`BRep`], so a body written before regions were stored (faces with
//! one `shell`, shells owned by a `solid`, a `solids` arena) loads with
//! its regions derived ([`BRep::derive_regions`]). Unknown fields are
//! refused, so a field added to `BRep` without one here fails the
//! round-trip test rather than being dropped on load; retired fields
//! are named and ignored.

use std::collections::{BTreeMap, HashMap};

use serde::Deserialize;
use serde::de::IgnoredAny;
use slotmap::SlotMap;

use crate::geometry::edge_triple::BoundarySide;
use crate::geometry::pcurve::Pcurve;
use crate::naming::{NamePool, OperationLabels};
use crate::tolerance::Tolerance;
use crate::topology::brep_face_use_trim_storage::FaceUseTrimKey;
use crate::topology::*;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct StoredBRep {
    vertices: SlotMap<VertexHandle, Vertex>,
    half_edges: SlotMap<HalfEdgeHandle, HalfEdge>,
    edges: SlotMap<EdgeHandle, Edge>,
    loops: SlotMap<LoopHandle, Loop>,
    faces: SlotMap<FaceHandle, Face>,
    shells: SlotMap<ShellHandle, Shell>,
    #[serde(alias = "solids")]
    regions: SlotMap<RegionHandle, Region>,
    #[serde(default)]
    outside: RegionHandle,
    tolerance: Tolerance,
    #[serde(default)]
    name_pool: NamePool,
    #[serde(default)]
    operation_labels: OperationLabels,
    #[serde(default = "BRep::default_max_imprecision")]
    max_imprecision: f64,
    #[serde(default, with = "super::serde_handle_map")]
    surface_boundary_hints: HashMap<HalfEdgeHandle, BoundarySide>,
    #[serde(default, with = "super::serde_handle_map")]
    surface_boundary_declines: HashMap<HalfEdgeHandle, SurfaceBoundaryHintDecline>,
    #[serde(default)]
    surface_boundary_audit_events: Vec<SurfaceBoundaryAuditEvent>,
    #[serde(default)]
    face_use_pcurves: BTreeMap<FaceUseTrimKey, Pcurve>,
    #[serde(default)]
    topology_contacts: Vec<TopologyContact>,
    /// Retired: the recorded free-edge set, read off the fin rings now.
    #[serde(default, rename = "boundary_edges")]
    _boundary_edges: IgnoredAny,
    /// Retired: the side map of non-manifold edges, stored on the fins'
    /// rings now.
    #[serde(default, rename = "nonmanifold_radial")]
    _nonmanifold_radial: IgnoredAny,
}

impl From<StoredBRep> for BRep {
    fn from(s: StoredBRep) -> Self {
        let mut b = BRep::new(s.tolerance);
        b.vertices = s.vertices;
        b.half_edges = s.half_edges;
        b.edges = s.edges;
        b.loops = s.loops;
        b.faces = s.faces;
        b.shells = s.shells;
        b.regions = s.regions;
        b.outside = s.outside;
        b.name_pool = s.name_pool;
        b.operation_labels = s.operation_labels;
        b.max_imprecision = s.max_imprecision;
        b.surface_boundary_hints = s.surface_boundary_hints;
        b.surface_boundary_declines = s.surface_boundary_declines;
        b.surface_boundary_audit_events = s.surface_boundary_audit_events;
        b.face_use_pcurves = s.face_use_pcurves;
        b.topology_contacts = s.topology_contacts;
        b.derive_regions();
        b
    }
}
