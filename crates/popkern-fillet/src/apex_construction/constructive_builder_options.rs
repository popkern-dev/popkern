//! Execution options for constructive apex topology application.

use super::cap_loop_construction_options::CapLoopConstructionOptions;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct ConstructiveApexExecutionOptions {
    cap_loop: CapLoopConstructionOptions,
    require_v_free_repeated_apex: bool,
}

impl ConstructiveApexExecutionOptions {
    pub(crate) fn from_environment_adapter() -> Self {
        Self {
            cap_loop: CapLoopConstructionOptions::from_environment_adapter(),
            require_v_free_repeated_apex: std::env::var_os(
                "POPKERN_A2_REQUIRE_V_FREE_REPEATED_APEX",
            )
            .is_some(),
        }
    }

    pub(super) fn cap_loop(self) -> CapLoopConstructionOptions {
        self.cap_loop
    }

    pub(super) fn require_v_free_repeated_apex(self) -> bool {
        self.require_v_free_repeated_apex
    }
}
