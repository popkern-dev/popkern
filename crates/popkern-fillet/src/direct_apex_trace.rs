//! Centralized direct-apex diagnostic trace switches.
//!
//! Thin topic wrappers over [`popkern_core::trace::trace_enabled`] so
//! call sites keep one name per channel.

pub(crate) fn enabled() -> bool {
    popkern_core::trace::trace_enabled("direct-apex")
}

pub(crate) fn cap_g0_enabled() -> bool {
    popkern_core::trace::trace_enabled("direct-apex-cap-g0")
}

pub(crate) fn prepass_enabled() -> bool {
    popkern_core::trace::trace_enabled("direct-apex-prepass")
}

pub(crate) fn support_frames_enabled() -> bool {
    popkern_core::trace::trace_enabled("direct-apex-support-frames")
}
