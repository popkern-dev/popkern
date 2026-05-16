//! Read-only loop and shell diagnostics for BRep repair decisions.

use std::collections::{HashMap, HashSet};

use crate::topology::brep_shell_diagnostics_facts::{
    attach_loop_shell_fact_ids, loop_shell_knowledge,
};
use crate::topology::brep_shell_diagnostics_helpers::push_finding;
use crate::topology::{
    BRep, BrepRepairFactId, EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, ShellHandle,
    VertexHandle,
};

/// Severity of a loop/shell diagnostic finding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BrepTopologyDiagnosticSeverity {
    /// Informational topology state useful for repair planning.
    Info,
    /// Suspicious topology that may still be representable.
    Warning,
    /// Invalid topology that should block repair mutation.
    Error,
}

/// Machine-readable loop/shell diagnostic kind.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum BrepTopologyDiagnosticKind {
    /// A stored handle points at no entity.
    MissingHandle,
    /// A half-edge loop walk did not close cleanly.
    BrokenLoopCycle,
    /// A loop references a face that does not own it.
    LoopFaceOwnershipMismatch,
    /// A face references a loop whose stored parent differs.
    FaceLoopOwnershipMismatch,
    /// A half-edge appears in more than one loop walk.
    DuplicateHalfEdgeMembership,
    /// A shell has no valid anchor face.
    MissingShellAnchor,
    /// A face references a missing shell.
    MissingFaceShell,
    /// A shell contains disconnected face components.
    FragmentedShell,
    /// Boundary edges do not form closed boundary cycles.
    DanglingBoundaryChain,
    /// Non-manifold radial data is malformed or present.
    RadialConflict,
}

/// One structured loop/shell diagnostic finding.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrepTopologyDiagnosticFinding {
    /// Finding severity.
    pub severity: BrepTopologyDiagnosticSeverity,
    /// Machine-readable finding kind.
    pub kind: BrepTopologyDiagnosticKind,
    /// Optional shell involved in the finding.
    pub shell: Option<ShellHandle>,
    /// Optional face involved in the finding.
    pub face: Option<FaceHandle>,
    /// Optional loop involved in the finding.
    pub loop_: Option<LoopHandle>,
    /// Optional edge involved in the finding.
    pub edge: Option<EdgeHandle>,
    /// Optional half-edge involved in the finding.
    pub half_edge: Option<HalfEdgeHandle>,
    /// Human-readable repair/debug message.
    pub message: String,
    /// Shared repair knowledge facts that explain this finding.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

/// Aggregate loop/shell diagnostic report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct BrepLoopShellDiagnosticReport {
    /// Number of loops in the BRep.
    pub loops: usize,
    /// Number of shells in the BRep.
    pub shells: usize,
    /// Number of faces in the BRep.
    pub faces: usize,
    /// Number of explicit boundary edges.
    pub boundary_edges: usize,
    /// Number of closed boundary cycles.
    pub boundary_cycles: usize,
    /// Number of open boundary chains/components.
    pub boundary_chains: usize,
    /// Number of connected face components across all shells.
    pub face_components: usize,
    /// Number of non-manifold radial entries.
    pub radial_entries: usize,
    /// Structured findings.
    pub findings: Vec<BrepTopologyDiagnosticFinding>,
    /// Shared repair knowledge facts referenced by this report.
    pub shared_fact_ids: Vec<BrepRepairFactId>,
}

impl BrepLoopShellDiagnosticReport {
    /// True when no error-severity findings were produced.
    pub fn has_no_errors(&self) -> bool {
        self.findings
            .iter()
            .all(|finding| finding.severity != BrepTopologyDiagnosticSeverity::Error)
    }
}

/// Analyze loop and shell topology without mutating the BRep.
pub fn analyze_brep_loop_shell_validity(brep: &BRep) -> BrepLoopShellDiagnosticReport {
    let knowledge = loop_shell_knowledge(brep);
    let mut findings = Vec::new();
    let loop_walks = analyze_loops(brep, &mut findings);
    analyze_face_loop_ownership(brep, &mut findings);
    analyze_shell_ownership(brep, &mut findings);
    let face_components = analyze_shell_components(brep, &mut findings);
    let boundary = analyze_boundary_graph(brep, &mut findings);
    analyze_radial_conflicts(brep, &mut findings);

    let report = BrepLoopShellDiagnosticReport {
        loops: brep.loops.len(),
        shells: brep.shells.len(),
        faces: brep.faces.len(),
        boundary_edges: brep.boundary_edges.len(),
        boundary_cycles: boundary.cycles,
        boundary_chains: boundary.chains,
        face_components,
        radial_entries: brep.nonmanifold_radial.len(),
        findings,
        shared_fact_ids: Vec::new(),
    }
    .with_duplicate_half_edge_findings(loop_walks);
    attach_loop_shell_fact_ids(report, brep, &knowledge)
}

#[derive(Debug, Default)]
struct LoopWalks {
    half_edge_loops: HashMap<HalfEdgeHandle, Vec<LoopHandle>>,
}

fn analyze_loops(brep: &BRep, findings: &mut Vec<BrepTopologyDiagnosticFinding>) -> LoopWalks {
    let mut walks = LoopWalks::default();
    let bound = brep.half_edges.len().max(1);
    for (loop_handle, loop_data) in &brep.loops {
        let mut cursor = loop_data.half_edge;
        let mut closed = false;
        for _ in 0..=bound {
            let Some(half_edge) = brep.half_edges.get(cursor) else {
                push_finding(
                    findings,
                    BrepTopologyDiagnosticSeverity::Error,
                    BrepTopologyDiagnosticKind::MissingHandle,
                    format!("loop {loop_handle:?} references missing half-edge {cursor:?}"),
                )
                .loop_(loop_handle)
                .half_edge(cursor);
                break;
            };
            walks
                .half_edge_loops
                .entry(cursor)
                .or_default()
                .push(loop_handle);
            if half_edge.loop_ != loop_handle {
                push_finding(
                    findings,
                    BrepTopologyDiagnosticSeverity::Error,
                    BrepTopologyDiagnosticKind::LoopFaceOwnershipMismatch,
                    format!("half-edge {cursor:?} belongs to loop {:?}", half_edge.loop_),
                )
                .loop_(loop_handle)
                .half_edge(cursor);
            }
            cursor = half_edge.next;
            if cursor == loop_data.half_edge {
                closed = true;
                break;
            }
        }
        if !closed {
            push_finding(
                findings,
                BrepTopologyDiagnosticSeverity::Error,
                BrepTopologyDiagnosticKind::BrokenLoopCycle,
                format!("loop {loop_handle:?} did not close in bounded walk"),
            )
            .loop_(loop_handle);
        }
        if !brep.faces.contains_key(loop_data.face) {
            push_finding(
                findings,
                BrepTopologyDiagnosticSeverity::Error,
                BrepTopologyDiagnosticKind::MissingHandle,
                format!(
                    "loop {loop_handle:?} references missing face {:?}",
                    loop_data.face
                ),
            )
            .loop_(loop_handle)
            .face(loop_data.face);
        }
    }
    walks
}

fn analyze_face_loop_ownership(brep: &BRep, findings: &mut Vec<BrepTopologyDiagnosticFinding>) {
    for (face_handle, face) in &brep.faces {
        check_face_loop(brep, findings, face_handle, face.outer_loop);
        for &loop_handle in &face.inner_loops {
            check_face_loop(brep, findings, face_handle, loop_handle);
        }
    }
}

fn check_face_loop(
    brep: &BRep,
    findings: &mut Vec<BrepTopologyDiagnosticFinding>,
    face: FaceHandle,
    loop_: LoopHandle,
) {
    let Some(loop_data) = brep.loops.get(loop_) else {
        push_finding(
            findings,
            BrepTopologyDiagnosticSeverity::Error,
            BrepTopologyDiagnosticKind::MissingHandle,
            format!("face {face:?} references missing loop {loop_:?}"),
        )
        .face(face)
        .loop_(loop_);
        return;
    };
    if loop_data.face != face {
        push_finding(
            findings,
            BrepTopologyDiagnosticSeverity::Error,
            BrepTopologyDiagnosticKind::FaceLoopOwnershipMismatch,
            format!(
                "face {face:?} references loop {loop_:?}, but loop parent is {:?}",
                loop_data.face
            ),
        )
        .face(face)
        .loop_(loop_);
    }
}

fn analyze_shell_ownership(brep: &BRep, findings: &mut Vec<BrepTopologyDiagnosticFinding>) {
    for (shell_handle, shell) in &brep.shells {
        if !brep.faces.contains_key(shell.face) {
            push_finding(
                findings,
                BrepTopologyDiagnosticSeverity::Error,
                BrepTopologyDiagnosticKind::MissingShellAnchor,
                format!(
                    "shell {shell_handle:?} references missing face {:?}",
                    shell.face
                ),
            )
            .shell(shell_handle)
            .face(shell.face);
        }
    }
    for (face_handle, face) in &brep.faces {
        if !brep.shells.contains_key(face.shell) {
            push_finding(
                findings,
                BrepTopologyDiagnosticSeverity::Error,
                BrepTopologyDiagnosticKind::MissingFaceShell,
                format!(
                    "face {face_handle:?} references missing shell {:?}",
                    face.shell
                ),
            )
            .face(face_handle)
            .shell(face.shell);
        }
    }
}

fn analyze_shell_components(
    brep: &BRep,
    findings: &mut Vec<BrepTopologyDiagnosticFinding>,
) -> usize {
    let mut adjacency: HashMap<FaceHandle, Vec<FaceHandle>> = HashMap::new();
    for face in brep.faces.keys() {
        adjacency.entry(face).or_default();
    }
    for (edge_handle, edge) in &brep.edges {
        if brep.boundary_edges.contains(&edge_handle) {
            continue;
        }
        let Some(a) = face_for_half_edge(brep, edge.half_edge) else {
            continue;
        };
        let Some(half_edge) = brep.half_edges.get(edge.half_edge) else {
            continue;
        };
        let Some(b) = face_for_half_edge(brep, half_edge.twin) else {
            continue;
        };
        if brep.faces[a].shell == brep.faces[b].shell {
            adjacency.entry(a).or_default().push(b);
            adjacency.entry(b).or_default().push(a);
        }
    }

    let mut total_components = 0usize;
    for shell in brep.shells.keys() {
        let faces: Vec<_> = brep
            .faces
            .iter()
            .filter_map(|(face, data)| (data.shell == shell).then_some(face))
            .collect();
        let components = count_components(&faces, &adjacency);
        total_components += components;
        if components > 1 {
            push_finding(
                findings,
                BrepTopologyDiagnosticSeverity::Warning,
                BrepTopologyDiagnosticKind::FragmentedShell,
                format!("shell {shell:?} has {components} disconnected face components"),
            )
            .shell(shell);
        }
    }
    total_components
}

#[derive(Debug, Default)]
struct BoundarySummary {
    cycles: usize,
    chains: usize,
}

fn analyze_boundary_graph(
    brep: &BRep,
    findings: &mut Vec<BrepTopologyDiagnosticFinding>,
) -> BoundarySummary {
    let mut graph: HashMap<VertexHandle, Vec<VertexHandle>> = HashMap::new();
    for &edge in &brep.boundary_edges {
        let Some((start, end)) = edge_vertices(brep, edge) else {
            push_finding(
                findings,
                BrepTopologyDiagnosticSeverity::Error,
                BrepTopologyDiagnosticKind::MissingHandle,
                format!("boundary edge {edge:?} has incomplete endpoint topology"),
            )
            .edge(edge);
            continue;
        };
        graph.entry(start).or_default().push(end);
        graph.entry(end).or_default().push(start);
    }
    let mut visited = HashSet::new();
    let mut summary = BoundarySummary::default();
    for &vertex in graph.keys() {
        if !visited.insert(vertex) {
            continue;
        }
        let mut stack = vec![vertex];
        let mut closed = true;
        while let Some(current) = stack.pop() {
            let neighbors = graph.get(&current).map(Vec::as_slice).unwrap_or(&[]);
            if neighbors.len() != 2 {
                closed = false;
            }
            for &next in neighbors {
                if visited.insert(next) {
                    stack.push(next);
                }
            }
        }
        if closed {
            summary.cycles += 1;
        } else {
            summary.chains += 1;
        }
    }
    if summary.chains > 0 {
        push_finding(
            findings,
            BrepTopologyDiagnosticSeverity::Warning,
            BrepTopologyDiagnosticKind::DanglingBoundaryChain,
            format!(
                "boundary graph has {} open chain components",
                summary.chains
            ),
        );
    }
    summary
}

fn analyze_radial_conflicts(brep: &BRep, findings: &mut Vec<BrepTopologyDiagnosticFinding>) {
    for (&from, &to) in &brep.nonmanifold_radial {
        if !brep.half_edges.contains_key(from) || !brep.half_edges.contains_key(to) {
            push_finding(
                findings,
                BrepTopologyDiagnosticSeverity::Error,
                BrepTopologyDiagnosticKind::RadialConflict,
                format!("radial entry {from:?} -> {to:?} references missing half-edge"),
            )
            .half_edge(from);
        } else {
            push_finding(
                findings,
                BrepTopologyDiagnosticSeverity::Warning,
                BrepTopologyDiagnosticKind::RadialConflict,
                format!("non-manifold radial entry {from:?} -> {to:?} is present"),
            )
            .half_edge(from);
        }
    }
}

impl BrepLoopShellDiagnosticReport {
    fn with_duplicate_half_edge_findings(mut self, loop_walks: LoopWalks) -> Self {
        for (half_edge, loops) in loop_walks.half_edge_loops {
            if loops.len() > 1 {
                self.findings.push(BrepTopologyDiagnosticFinding {
                    severity: BrepTopologyDiagnosticSeverity::Error,
                    kind: BrepTopologyDiagnosticKind::DuplicateHalfEdgeMembership,
                    shell: None,
                    face: None,
                    loop_: loops.first().copied(),
                    edge: None,
                    half_edge: Some(half_edge),
                    message: format!(
                        "half-edge {half_edge:?} appears in {} loop walks",
                        loops.len()
                    ),
                    shared_fact_ids: Vec::new(),
                });
            }
        }
        self
    }
}

fn face_for_half_edge(brep: &BRep, half_edge: HalfEdgeHandle) -> Option<FaceHandle> {
    let half_edge = brep.half_edges.get(half_edge)?;
    Some(brep.loops.get(half_edge.loop_)?.face)
}

fn edge_vertices(brep: &BRep, edge: EdgeHandle) -> Option<(VertexHandle, VertexHandle)> {
    let half_edge = brep.half_edges.get(brep.edges.get(edge)?.half_edge)?;
    let next = brep.half_edges.get(half_edge.next)?;
    Some((half_edge.origin, next.origin))
}

fn count_components(
    faces: &[FaceHandle],
    adjacency: &HashMap<FaceHandle, Vec<FaceHandle>>,
) -> usize {
    let face_set: HashSet<_> = faces.iter().copied().collect();
    let mut visited = HashSet::new();
    let mut components = 0usize;
    for &face in faces {
        if !visited.insert(face) {
            continue;
        }
        components += 1;
        let mut stack = vec![face];
        while let Some(current) = stack.pop() {
            for &next in adjacency.get(&current).map(Vec::as_slice).unwrap_or(&[]) {
                if face_set.contains(&next) && visited.insert(next) {
                    stack.push(next);
                }
            }
        }
    }
    components
}

#[cfg(test)]
#[path = "brep_shell_diagnostics_tests.rs"]
mod tests;
