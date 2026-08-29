//! Execution options for the repeated-apex partition operator.

#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub(crate) struct RepeatedApexPartitionExecutionOptions {
    pub(in crate::apex_construction) apply_cylinder_trim_rewrite: bool,
    pub(in crate::apex_construction) require_cylinder_trim_rewrite_ready: bool,
    pub(in crate::apex_construction) require_cylinder_trim_rebuild_plan: bool,
    pub(in crate::apex_construction) adopt_cylinder_fillet_end_endpoints: bool,
    pub(in crate::apex_construction) require_shared_cylinder_fillet_end: bool,
    pub(in crate::apex_construction) require_v_free_repeated_apex: bool,
}

impl RepeatedApexPartitionExecutionOptions {
    pub(crate) fn new(
        apply_cylinder_trim_rewrite: bool,
        require_cylinder_trim_rewrite_ready: bool,
        require_cylinder_trim_rebuild_plan: bool,
        adopt_cylinder_fillet_end_endpoints: bool,
        require_shared_cylinder_fillet_end: bool,
        require_v_free_repeated_apex: bool,
    ) -> Self {
        Self {
            apply_cylinder_trim_rewrite,
            require_cylinder_trim_rewrite_ready,
            require_cylinder_trim_rebuild_plan,
            adopt_cylinder_fillet_end_endpoints,
            require_shared_cylinder_fillet_end,
            require_v_free_repeated_apex,
        }
    }

    /// Copy of these options with the V-free repeated-apex contract
    /// resolved for a specific apex.
    pub(crate) fn with_v_free_repeated_apex(mut self, enabled: bool) -> Self {
        self.require_v_free_repeated_apex = enabled;
        self
    }
}
