//! Result payloads from executable side-support continuation.

use popkern_core::topology::{FaceHandle, LoopHandle};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct SideSupportRegionContinuationResult {
    bridge_faces: Vec<FaceHandle>,
    bridge_loops: Vec<LoopHandle>,
}

impl SideSupportRegionContinuationResult {
    pub(super) fn not_applied() -> Self {
        Self::default()
    }

    pub(super) fn from_bridge_payload(
        bridge_faces: Vec<FaceHandle>,
        bridge_loops: Vec<LoopHandle>,
    ) -> Self {
        Self {
            bridge_faces,
            bridge_loops,
        }
    }

    pub(super) fn was_applied(&self) -> bool {
        !self.bridge_faces.is_empty() || !self.bridge_loops.is_empty()
    }

    pub(super) fn bridge_loops(&self) -> &[LoopHandle] {
        &self.bridge_loops
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "side_support_result=applied={} bridge_faces={:?} bridge_loops={:?}",
            self.was_applied(),
            self.bridge_faces,
            self.bridge_loops
        )
    }
}
