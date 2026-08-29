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
    /// Copy of these options with the V-free repeated-apex contract
    /// resolved (the mode's valence scoping lives in
    /// `ConstructiveApexExecutionOptions`).
    pub(crate) fn with_v_free_repeated_apex(mut self, enabled: bool) -> Self {
        self.repeated_apex_partition = self
            .repeated_apex_partition
            .with_v_free_repeated_apex(enabled);
        self
    }
}
