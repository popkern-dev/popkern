//! Execution options for cap-loop construction.

use super::repeated_apex_partition::RepeatedApexPartitionExecutionOptions;

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct CapLoopConstructionOptions {
    pub(super) allow_repeated_apex_mef_fallback: bool,
    pub(super) repeated_apex_partition: RepeatedApexPartitionExecutionOptions,
}

impl CapLoopConstructionOptions {
    pub(crate) fn from_environment_adapter() -> Self {
        Self {
            allow_repeated_apex_mef_fallback: std::env::var_os(
                "POPKERN_A2_ALLOW_REPEATED_APEX_LOOP",
            )
            .is_some(),
            repeated_apex_partition: RepeatedApexPartitionExecutionOptions::new(
                std::env::var_os("POPKERN_A2_APPLY_CYLINDER_TRIM_REWRITE").is_some(),
                std::env::var_os("POPKERN_A2_REQUIRE_CYLINDER_TRIM_REWRITE_READY").is_some(),
                std::env::var_os("POPKERN_A2_REQUIRE_CYLINDER_TRIM_REBUILD_PLAN").is_some(),
                std::env::var_os("POPKERN_A2_ADOPT_CYLINDER_FILLET_END_ENDPOINTS").is_some(),
                std::env::var_os("POPKERN_A2_REQUIRE_SHARED_CYLINDER_FILLET_END").is_some(),
            ),
        }
    }
}
