use super::*;

impl FaceFamilyBridgeRegionPrimitivePlan {
    pub(in crate::apex_construction) fn support_loop(&self) -> LoopHandle {
        self.support_loop
    }

    pub(in crate::apex_construction) fn outside_loop(&self) -> LoopHandle {
        self.outside_loop
    }

    pub(in crate::apex_construction) fn support_trims(&self) -> &[BridgeSupportTrimPlan; 2] {
        &self.trims
    }

    pub(in crate::apex_construction) fn support_paths(&self) -> &FaceFamilyBridgeLoopIntervals {
        &self.support_paths
    }

    pub(in crate::apex_construction) fn outside_paths(&self) -> &FaceFamilyBridgeLoopIntervals {
        &self.outside_paths
    }

    pub(in crate::apex_construction) fn boundary_cycles(&self) -> &[FaceFamilyBridgeBoundaryCycle] {
        &self.boundary_cycles
    }

    pub(in crate::apex_construction) fn boundary_exposure(
        &self,
    ) -> &FaceFamilyBridgeBoundaryExposurePlan {
        &self.boundary_exposure
    }

    pub(in crate::apex_construction) fn boundary_readiness(
        &self,
    ) -> &FaceFamilyBridgeBoundaryReadiness {
        &self.boundary_readiness
    }

    pub(in crate::apex_construction) fn expected_edges(&self) -> usize {
        self.expected_edges
    }

    pub(in crate::apex_construction) fn expected_faces(&self) -> usize {
        self.expected_faces
    }

    pub(in crate::apex_construction) fn expected_loops(&self) -> usize {
        self.expected_loops
    }

    pub(in crate::apex_construction) fn with_boundary_cycles(
        mut self,
        brep: &BRep,
        boundary_cycles: Vec<FaceFamilyBridgeBoundaryCycle>,
    ) -> Self {
        let boundary_exposure = FaceFamilyBridgeBoundaryExposurePlan::build(brep, &boundary_cycles);
        self.boundary_readiness = boundary_exposure.readiness().clone();
        self.expected_edges = boundary_cycles.len() * 2
            + boundary_exposure.fragment_retire_count()
            + boundary_exposure.self_owned_exposure_count();
        self.expected_faces = boundary_cycles.len();
        self.expected_loops = boundary_cycles.len();
        self.boundary_cycles = boundary_cycles;
        self.boundary_exposure = boundary_exposure;
        self
    }

    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
        format!(
            "support={:?} support_loop={:?} outside={:?}/{:?} trims={} support_paths={} outside_paths={} boundary_readiness={} boundary_exposure={} cycles=[{}] delta=dV=0 dE={} dF={} dL={}",
            self.support_face,
            self.support_loop,
            self.outside_face,
            self.outside_loop,
            self.trims.len(),
            self.support_paths.diagnostic_label(),
            self.outside_paths.diagnostic_label(),
            self.boundary_readiness.diagnostic_label(),
            self.boundary_exposure.diagnostic_label(),
            self.boundary_cycles
                .iter()
                .map(FaceFamilyBridgeBoundaryCycle::diagnostic_label)
                .collect::<Vec<_>>()
                .join("; "),
            self.expected_edges,
            self.expected_faces,
            self.expected_loops
        )
    }
}

impl FaceFamilyBridgeBoundaryCycle {
    pub(in crate::apex_construction) fn new(
        support_interval: Vec<HalfEdgeHandle>,
        outside_interval: Vec<HalfEdgeHandle>,
    ) -> Self {
        Self {
            support_interval,
            outside_interval,
        }
    }

    pub(in crate::apex_construction) fn support_interval(&self) -> &[HalfEdgeHandle] {
        &self.support_interval
    }

    pub(in crate::apex_construction) fn outside_interval(&self) -> &[HalfEdgeHandle] {
        &self.outside_interval
    }

    pub(in crate::apex_construction) fn is_empty(&self) -> bool {
        self.support_interval.len() < 2 || self.outside_interval.len() < 2
    }

    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
        format!(
            "support_interval={} outside_interval={}",
            self.support_interval.len(),
            self.outside_interval.len()
        )
    }
}

impl FaceFamilyBridgeLoopIntervals {
    pub(in crate::apex_construction) fn first_to_second(&self) -> &[HalfEdgeHandle] {
        &self.first_to_second
    }

    pub(in crate::apex_construction) fn second_to_first(&self) -> &[HalfEdgeHandle] {
        &self.second_to_first
    }

    pub(in crate::apex_construction) fn diagnostic_label(&self) -> String {
        format!(
            "loop_len={} first_to_second={} second_to_first={}",
            self.loop_len,
            self.first_to_second.len(),
            self.second_to_first.len()
        )
    }

    pub(in crate::apex_construction) fn are_complementary(&self) -> bool {
        self.first_to_second.len() + self.second_to_first.len() == self.loop_len + 2
    }
}

impl FaceFamilyBridgeOperandGap {
    pub(in crate::apex_construction) fn diagnostic_label(&self) -> &'static str {
        match self {
            FaceFamilyBridgeOperandGap::UnsupportedTrimCount => "unsupported-trim-count",
            FaceFamilyBridgeOperandGap::MixedSupportLoops => "mixed-support-loops",
            FaceFamilyBridgeOperandGap::MixedOutsideLoops => "mixed-outside-loops",
            FaceFamilyBridgeOperandGap::MissingSupportPath => "missing-support-path",
            FaceFamilyBridgeOperandGap::MissingOutsidePath => "missing-outside-path",
            FaceFamilyBridgeOperandGap::NonComplementarySupportIntervals => {
                "non-complementary-support-intervals"
            }
            FaceFamilyBridgeOperandGap::NonComplementaryOutsideIntervals => {
                "non-complementary-outside-intervals"
            }
            FaceFamilyBridgeOperandGap::SharedBridgeOperandEdges => "shared-bridge-operand-edges",
        }
    }
}

pub(super) fn two_bridge_trims(
    trims: &[BridgeSupportTrimPlan],
) -> Option<[&BridgeSupportTrimPlan; 2]> {
    match trims {
        [first, second] => Some([first, second]),
        _ => None,
    }
}

pub(super) fn full_loop_interval_from(
    brep: &BRep,
    loop_handle: LoopHandle,
    start: HalfEdgeHandle,
) -> Option<Vec<HalfEdgeHandle>> {
    let mut interval = Vec::new();
    let mut cursor = start;
    for _ in 0..=brep.half_edges.len() {
        let half_edge = brep.half_edges.get(cursor)?;
        if half_edge.loop_ != loop_handle {
            return None;
        }
        interval.push(cursor);
        cursor = half_edge.next;
        if cursor == start {
            interval.push(start);
            return Some(interval);
        }
    }
    None
}

pub(super) fn shorter_interval(
    first: &[HalfEdgeHandle],
    second: &[HalfEdgeHandle],
) -> Vec<HalfEdgeHandle> {
    if first.len() <= second.len() {
        first.to_vec()
    } else {
        second.to_vec()
    }
}
