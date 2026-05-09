//! Execution options for cap-loop construction.

use super::repeated_apex_partition::RepeatedApexPartitionExecutionOptions;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct CapLoopConstructionOptions {
    pub(super) repeated_apex_partition: RepeatedApexPartitionExecutionOptions,
}

impl Default for CapLoopConstructionOptions {
    fn default() -> Self {
        Self {
            repeated_apex_partition: RepeatedApexPartitionExecutionOptions::new(
                true, false, false, false, false, false,
            ),
        }
    }
}

impl CapLoopConstructionOptions {
    pub(crate) fn from_environment_adapter() -> Self {
        Self {
            repeated_apex_partition: RepeatedApexPartitionExecutionOptions::new(
                true,
                std::env::var_os("POPKERN_A2_REQUIRE_CYLINDER_TRIM_REWRITE_READY").is_some(),
                std::env::var_os("POPKERN_A2_REQUIRE_CYLINDER_TRIM_REBUILD_PLAN").is_some(),
                std::env::var_os("POPKERN_A2_ADOPT_CYLINDER_FILLET_END_ENDPOINTS").is_some(),
                std::env::var_os("POPKERN_A2_REQUIRE_SHARED_CYLINDER_FILLET_END").is_some(),
                std::env::var_os("POPKERN_A2_REQUIRE_V_FREE_REPEATED_APEX").is_some(),
            ),
        }
    }
}
