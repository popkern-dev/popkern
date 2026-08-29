//! Return type for Gregory cap surface assembly.

use crate::types::{ApexConstructionDiagnostic, CornerPatchInfo};

pub(super) struct GregoryCapSurfaceAssembly {
    pub corner_patches: Vec<CornerPatchInfo>,
    pub diagnostics: Vec<ApexConstructionDiagnostic>,
}
