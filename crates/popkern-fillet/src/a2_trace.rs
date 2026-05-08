//! Centralized A2 diagnostic trace switches.

pub(crate) fn enabled() -> bool {
    std::env::var_os("POPKERN_A2_TRACE").is_some()
}

pub(crate) fn cap_g0_enabled() -> bool {
    std::env::var_os("POPKERN_A2_TRACE_CAP_G0").is_some()
}

pub(crate) fn prepass_enabled() -> bool {
    std::env::var_os("POPKERN_A2_PREPASS_TRACE").is_some()
}

pub(crate) fn support_frames_enabled() -> bool {
    std::env::var_os("POPKERN_A2_TRACE_SUPPORT_FRAMES").is_some()
}
