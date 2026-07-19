//! Result payloads from executable side-support continuation.

use popkern_core::topology::{FaceHandle, LoopHandle};

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) struct SideSupportRegionContinuationResult {
    primitive: SideSupportContinuationPrimitive,
    bridge_faces: Vec<FaceHandle>,
    bridge_loops: Vec<LoopHandle>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub(super) enum SideSupportContinuationPrimitive {
    #[default]
    None,
    FaceFamilyBridge,
    OutsideOutsideCylinderMediated,
    /// Mixed plan: face-family bridge regions and cylinder-mediated
    /// continuations applied from one selection (the cube_extreme
    /// CrossOutsideOwners corner).
    MixedBridgeAndCylinderMediated,
}

impl SideSupportRegionContinuationResult {
    pub(super) fn not_applied() -> Self {
        Self::default()
    }

    pub(super) fn from_face_family_bridge_payload(
        bridge_faces: Vec<FaceHandle>,
        bridge_loops: Vec<LoopHandle>,
    ) -> Self {
        Self::from_bridge_payload(
            SideSupportContinuationPrimitive::FaceFamilyBridge,
            bridge_faces,
            bridge_loops,
        )
    }

    pub(super) fn from_outside_outside_cylinder_payload(
        bridge_faces: Vec<FaceHandle>,
        bridge_loops: Vec<LoopHandle>,
    ) -> Self {
        Self::from_bridge_payload(
            SideSupportContinuationPrimitive::OutsideOutsideCylinderMediated,
            bridge_faces,
            bridge_loops,
        )
    }

    fn from_bridge_payload(
        primitive: SideSupportContinuationPrimitive,
        bridge_faces: Vec<FaceHandle>,
        bridge_loops: Vec<LoopHandle>,
    ) -> Self {
        Self {
            primitive,
            bridge_faces,
            bridge_loops,
        }
    }

    /// Combine two applied payloads from one mixed selection.
    pub(super) fn merged_with(self, other: Self) -> Self {
        let mut bridge_faces = self.bridge_faces;
        let mut bridge_loops = self.bridge_loops;
        bridge_faces.extend(other.bridge_faces);
        bridge_loops.extend(other.bridge_loops);
        Self {
            primitive: SideSupportContinuationPrimitive::MixedBridgeAndCylinderMediated,
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

    pub(super) fn bridge_face_count(&self) -> usize {
        self.bridge_faces.len()
    }

    pub(super) fn bridge_loop_count(&self) -> usize {
        self.bridge_loops.len()
    }

    pub(super) fn primitive_code(&self) -> &'static str {
        self.primitive.code()
    }

    pub(super) fn diagnostic_label(&self) -> String {
        format!(
            "side_support_result=primitive={} applied={} bridge_faces={:?} bridge_loops={:?}",
            self.primitive_code(),
            self.was_applied(),
            self.bridge_faces,
            self.bridge_loops
        )
    }
}

impl SideSupportContinuationPrimitive {
    fn code(&self) -> &'static str {
        match self {
            SideSupportContinuationPrimitive::None => "none",
            SideSupportContinuationPrimitive::FaceFamilyBridge => "face_family_bridge",
            SideSupportContinuationPrimitive::OutsideOutsideCylinderMediated => {
                "outside_outside_cylinder_mediated"
            }
            SideSupportContinuationPrimitive::MixedBridgeAndCylinderMediated => {
                "face_family_bridge+outside_outside_cylinder_mediated"
            }
        }
    }
}
