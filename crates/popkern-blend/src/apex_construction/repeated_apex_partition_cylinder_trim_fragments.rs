use super::*;

#[derive(Clone)]
pub(super) struct RetiredCylinderTrimFragment {
    pub(super) path: Vec<HalfEdgeHandle>,
    pub(super) vertices: Vec<VertexHandle>,
    pub(super) adjacent_face_runs: Vec<RetiredCylinderTrimRun>,
    pub(super) adjacent_faces: Vec<FaceHandle>,
    pub(super) off_surface_vertices: Vec<VertexHandle>,
    pub(super) materialized_anchor_replacement: bool,
}

#[derive(Clone)]
pub(super) struct RetiredCylinderTrimRun {
    pub(super) adjacent_face: FaceHandle,
    pub(super) path: Vec<HalfEdgeHandle>,
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) vertices: Vec<VertexHandle>,
    pub(super) off_surface_vertices: Vec<VertexHandle>,
}

pub(super) struct ReplacementRunReownershipPlan {
    pub(super) spans: Vec<ReplacementRunSpanOwnership>,
}

pub(super) struct ReplacementRunSpanOwnership {
    pub(super) start: VertexHandle,
    pub(super) end: VertexHandle,
    pub(super) candidate_faces: Vec<ReplacementRunSupportCandidate>,
}

pub(super) struct ReplacementRunSupportCandidate {
    pub(super) face: FaceHandle,
    pub(super) curve_kind: &'static str,
}

pub(super) struct FragmentAuthoringStatus {
    pub(super) authorable: bool,
    pub(super) diagnostic: String,
}

impl RetiredCylinderTrimFragment {
    pub(super) fn from_path(brep: &BRep, face: FaceHandle, path: Vec<HalfEdgeHandle>) -> Self {
        let vertices = path
            .iter()
            .filter_map(|half_edge| brep.half_edges.get(*half_edge).map(|data| data.origin))
            .collect();
        let adjacent_face_runs = retired_trim_runs_by_adjacent_face(brep, face, &path);
        let adjacent_faces = path
            .iter()
            .filter_map(|half_edge| {
                let twin = brep.half_edges.get(*half_edge)?.twin;
                let twin_loop = brep.half_edges.get(twin)?.loop_;
                let twin_face = brep.loops.get(twin_loop)?.face;
                (twin_face != face).then_some(twin_face)
            })
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        let off_surface_vertices = path_off_cylinder_vertices(brep, face, &path);
        Self {
            path,
            vertices,
            adjacent_face_runs,
            adjacent_faces,
            off_surface_vertices,
            materialized_anchor_replacement: false,
        }
    }

    pub(super) fn materialized_anchor_replacement(start: VertexHandle, end: VertexHandle) -> Self {
        Self {
            path: Vec::new(),
            vertices: vec![start, end],
            adjacent_face_runs: Vec::new(),
            adjacent_faces: Vec::new(),
            off_surface_vertices: Vec::new(),
            materialized_anchor_replacement: true,
        }
    }

    pub(super) fn diagnostic_label(&self) -> String {
        if self.materialized_anchor_replacement {
            return format!(
                "materialized-anchor replacement {:?}->{:?}",
                self.vertices.first().copied().unwrap_or_default(),
                self.vertices.get(1).copied().unwrap_or_default()
            );
        }
        format!(
            "{} half-edge(s), vertices {:?}, adjacent faces {:?}, runs [{}], \
             off-cylinder vertices {:?}",
            self.path.len(),
            self.vertices,
            self.adjacent_faces,
            self.adjacent_face_runs
                .iter()
                .map(RetiredCylinderTrimRun::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.off_surface_vertices
        )
    }
}

impl RetiredCylinderTrimRun {
    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "{:?}: {} half-edge(s), {:?}->{:?}, vertices {:?}, off-cylinder {:?}, action {}",
            self.adjacent_face,
            self.path.len(),
            self.start,
            self.end,
            self.vertices,
            self.off_surface_vertices,
            if self.requires_replacement() {
                "replace"
            } else {
                "reuse"
            }
        )
    }

    pub(super) fn requires_replacement(&self) -> bool {
        !self.off_surface_vertices.is_empty()
    }
}

impl ReplacementRunReownershipPlan {
    pub(super) fn diagnostic_label(&self) -> String {
        self.spans
            .iter()
            .map(ReplacementRunSpanOwnership::diagnostic_label)
            .collect::<Vec<_>>()
            .join(", ")
    }
}

impl ReplacementRunSpanOwnership {
    fn diagnostic_label(&self) -> String {
        let candidates = self
            .candidate_faces
            .iter()
            .map(|candidate| format!("{:?}/{}", candidate.face, candidate.curve_kind))
            .collect::<Vec<_>>()
            .join("|");
        format!(
            "{:?}->{:?} candidates [{}]",
            self.start, self.end, candidates
        )
    }
}
