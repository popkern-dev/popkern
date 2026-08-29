//! Monotonic stopwatch that cannot trap on clockless platforms.
//!
//! `std::time::Instant::now()` panics at runtime on
//! `wasm32-unknown-unknown` (the browser target has no ambient
//! monotonic clock), which made every timed kernel pipeline stage a
//! runtime trap in the browser (measured under PK-369: the first
//! `extrude` aborted with `unreachable`). Kernel code therefore never
//! calls `Instant::now()` directly; it uses [`Stopwatch`], which
//! degrades to *absent* timing instead of aborting:
//!
//! - Native and WASI targets: wraps `std::time::Instant`; timing is
//!   always available.
//! - `wasm32-unknown-unknown`: reads a host-registered millisecond
//!   clock (the browser playground registers `performance.now` at
//!   startup via [`register_monotonic_clock`]). With no registered
//!   clock, [`Stopwatch::elapsed`] returns `None` — telemetry is
//!   honestly absent, never a fabricated zero.
//!
//! Timing here is observability (stage events, telemetry records,
//! trace lines), never control flow: kernel algorithms must not make
//! decisions on elapsed wall-clock time, so an absent clock degrades
//! reports only.

use std::time::Duration;

/// Host-registered millisecond clock for targets without an ambient
/// monotonic clock. Set once per process.
#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
static HOST_CLOCK_MS: std::sync::OnceLock<fn() -> f64> = std::sync::OnceLock::new();

/// Register the host's monotonic millisecond clock (browser target
/// only; e.g. a wasm-bindgen shim over `performance.now`). Returns
/// `false` if a clock was already registered (first registration
/// wins).
#[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
pub fn register_monotonic_clock(clock_ms: fn() -> f64) -> bool {
    HOST_CLOCK_MS.set(clock_ms).is_ok()
}

/// Whether elapsed-time measurement is available in this process.
/// Always `true` off the browser target; on
/// `wasm32-unknown-unknown`, `true` once a host clock is registered.
pub fn monotonic_clock_available() -> bool {
    #[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
    {
        HOST_CLOCK_MS.get().is_some()
    }
    #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
    {
        true
    }
}

/// A started stopwatch. See the module docs for the platform
/// contract.
#[derive(Debug, Clone, Copy)]
pub struct Stopwatch {
    #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
    started: std::time::Instant,
    /// Start timestamp in host-clock milliseconds; `None` when no
    /// host clock was registered at start time.
    #[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
    started_ms: Option<f64>,
}

impl Stopwatch {
    /// Start measuring now.
    pub fn start() -> Stopwatch {
        #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
        {
            Stopwatch {
                started: std::time::Instant::now(),
            }
        }
        #[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
        {
            Stopwatch {
                started_ms: HOST_CLOCK_MS.get().map(|clock| clock()),
            }
        }
    }

    /// Elapsed time since [`Stopwatch::start`], or `None` when no
    /// monotonic clock exists on this platform (clockless browser
    /// target). Telemetry consumers record the absence, not a zero.
    pub fn elapsed(&self) -> Option<Duration> {
        #[cfg(not(all(target_arch = "wasm32", target_os = "unknown")))]
        {
            Some(self.started.elapsed())
        }
        #[cfg(all(target_arch = "wasm32", target_os = "unknown"))]
        {
            let started = self.started_ms?;
            let clock = HOST_CLOCK_MS.get()?;
            Some(Duration::from_secs_f64(
                (clock() - started).max(0.0) / 1.0e3,
            ))
        }
    }

    /// Elapsed milliseconds, or `None` when unmeasurable (see
    /// [`Stopwatch::elapsed`]).
    pub fn elapsed_ms(&self) -> Option<f64> {
        self.elapsed().map(|d| d.as_secs_f64() * 1.0e3)
    }

    /// Elapsed whole milliseconds for human-readable diagnostic
    /// detail strings ONLY; reads 0 on clockless platforms.
    /// Structured numeric telemetry must use [`Stopwatch::elapsed`]
    /// and record absence instead of a zero.
    pub fn elapsed_or_zero_ms(&self) -> u128 {
        self.elapsed().map_or(0, |d| d.as_millis())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stopwatch_measures_on_clocked_platforms() {
        let sw = Stopwatch::start();
        assert!(monotonic_clock_available());
        let elapsed = sw.elapsed().expect("native clock is always present");
        assert!(elapsed >= Duration::ZERO);
        assert!(sw.elapsed_ms().expect("ms present") >= 0.0);
    }
}
