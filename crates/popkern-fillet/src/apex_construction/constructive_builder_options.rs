//! Execution options for constructive apex topology application.

use super::cap_loop_construction_options::CapLoopConstructionOptions;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct ConstructiveApexExecutionOptions {
    cap_loop: CapLoopConstructionOptions,
}

impl ConstructiveApexExecutionOptions {
    pub(crate) fn from_environment_adapter() -> Self {
        Self {
            cap_loop: CapLoopConstructionOptions::from_environment_adapter(),
        }
    }

    pub(super) fn cap_loop(self) -> CapLoopConstructionOptions {
        self.cap_loop
    }
}
