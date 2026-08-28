//! Centralized direct-apex diagnostic trace switches.

fn trace_switch_enabled(primary: &str, alias: &str) -> bool {
    // Compatibility alias for pre-DirectApex local workflows. Delete once
    // internal docs and task PK-038 no longer reference the retired A2 switch names.
    std::env::var_os(primary).is_some() || std::env::var_os(alias).is_some()
}

pub(crate) fn enabled() -> bool {
    trace_switch_enabled("POPKERN_DIRECT_APEX_TRACE", "POPKERN_A2_TRACE")
}

pub(crate) fn cap_g0_enabled() -> bool {
    trace_switch_enabled(
        "POPKERN_DIRECT_APEX_TRACE_CAP_G0",
        "POPKERN_A2_TRACE_CAP_G0",
    )
}

pub(crate) fn prepass_enabled() -> bool {
    trace_switch_enabled(
        "POPKERN_DIRECT_APEX_PREPASS_TRACE",
        "POPKERN_A2_PREPASS_TRACE",
    )
}

pub(crate) fn support_frames_enabled() -> bool {
    trace_switch_enabled(
        "POPKERN_DIRECT_APEX_TRACE_SUPPORT_FRAMES",
        "POPKERN_A2_TRACE_SUPPORT_FRAMES",
    )
}
