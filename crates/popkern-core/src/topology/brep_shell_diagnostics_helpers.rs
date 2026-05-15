//! Helper builder for BRep shell diagnostic findings.

use crate::topology::{
    BrepTopologyDiagnosticFinding, BrepTopologyDiagnosticKind, BrepTopologyDiagnosticSeverity,
    EdgeHandle, FaceHandle, HalfEdgeHandle, LoopHandle, ShellHandle,
};

pub(super) fn push_finding(
    findings: &mut Vec<BrepTopologyDiagnosticFinding>,
    severity: BrepTopologyDiagnosticSeverity,
    kind: BrepTopologyDiagnosticKind,
    message: String,
) -> FindingBuilder<'_> {
    FindingBuilder {
        findings,
        finding: BrepTopologyDiagnosticFinding {
            severity,
            kind,
            shell: None,
            face: None,
            loop_: None,
            edge: None,
            half_edge: None,
            message,
        },
    }
}

pub(super) struct FindingBuilder<'a> {
    findings: &'a mut Vec<BrepTopologyDiagnosticFinding>,
    finding: BrepTopologyDiagnosticFinding,
}

impl FindingBuilder<'_> {
    pub(super) fn shell(mut self, shell: ShellHandle) -> Self {
        self.finding.shell = Some(shell);
        self
    }

    pub(super) fn face(mut self, face: FaceHandle) -> Self {
        self.finding.face = Some(face);
        self
    }

    pub(super) fn loop_(mut self, loop_: LoopHandle) -> Self {
        self.finding.loop_ = Some(loop_);
        self
    }

    pub(super) fn edge(mut self, edge: EdgeHandle) -> Self {
        self.finding.edge = Some(edge);
        self
    }

    pub(super) fn half_edge(mut self, half_edge: HalfEdgeHandle) -> Self {
        self.finding.half_edge = Some(half_edge);
        self
    }
}

impl Drop for FindingBuilder<'_> {
    fn drop(&mut self) {
        self.findings.push(self.finding.clone());
    }
}
