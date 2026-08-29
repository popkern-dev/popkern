//! Structured timing helpers for KP20 polar dispatch diagnostics.

#[derive(Debug, Clone, Copy, Default)]
pub(super) struct Kp20DispatchStageTiming {
    pub(super) support_domain_ms: u128,
    pub(super) support_transition_ms: u128,
    pub(super) support_normal_lift_ms: u128,
    pub(super) origin_inner_ms: u128,
    pub(super) paired_inner_ms: u128,
    pub(super) support_inner_ms: u128,
    pub(super) support_reparameterization_ms: u128,
}

impl Kp20DispatchStageTiming {
    pub(super) fn total_ms(self) -> u128 {
        self.support_domain_ms
            + self.support_transition_ms
            + self.support_normal_lift_ms
            + self.origin_inner_ms
            + self.paired_inner_ms
            + self.support_inner_ms
            + self.support_reparameterization_ms
    }

    pub(super) fn label(self) -> String {
        format!(
            "timing_total_ms={};timing_support_domain_ms={};timing_support_transition_ms={};timing_support_normal_lift_ms={};timing_origin_inner_ms={};timing_paired_inner_ms={};timing_support_inner_ms={};timing_support_reparameterization_ms={}",
            self.total_ms(),
            self.support_domain_ms,
            self.support_transition_ms,
            self.support_normal_lift_ms,
            self.origin_inner_ms,
            self.paired_inner_ms,
            self.support_inner_ms,
            self.support_reparameterization_ms,
        )
    }
}

pub(super) fn timed<T>(f: impl FnOnce() -> T) -> (T, u128) {
    let start = popkern_core::time::Stopwatch::start();
    let value = f();
    (value, start.elapsed_or_zero_ms())
}
