//! Peak resident memory of the six-axis arm build in the test profile.
//!
//! On 2026-10-14 every arm test held 13.5-16.5 GB: debug-assertion
//! builds ran `validate_all()` after every boolean purely to print
//! stderr, on rayon's long-lived workers, and each worker kept ~1 GB.
//! Nine arm tests in one nextest sweep (~125 GB) exhausted the host's
//! RAM and swap. Without that diagnostic the build peaks at 0.86 GB
//! (16 rayon workers). No gate measured memory, so a 17x regression
//! hid behind passing tests that were merely slow.
//!
//! This binary holds ONE test so the process peak is the build's
//! alone (`cargo test` runs a binary's tests as threads in one
//! process). The budget is ~3x the measured peak: allocator variance
//! is small, and the old defect retained ~1 GB per worker, so it
//! fails here on any machine with two or more cores. Linux only
//! (`/proc/self/status`); elsewhere it skips loudly.

mod arm_support;

/// Peak resident set budget for one arm build (bytes).
const PEAK_BUDGET: u64 = 2_500 * 1024 * 1024;

/// `VmHWM` (peak resident set) of this process, in bytes.
fn peak_resident() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    let line = status.lines().find(|l| l.starts_with("VmHWM:"))?;
    let kb: u64 = line.split_whitespace().nth(1)?.parse().ok()?;
    Some(kb * 1024)
}

#[test]
fn the_arm_builds_within_its_peak_memory_budget() {
    if peak_resident().is_none() {
        eprintln!("[arm-memory] SKIP: no /proc/self/status VmHWM on this platform");
        return;
    }
    let asm = arm_support::build();
    let peak = peak_resident().unwrap_or(u64::MAX);
    let workers = std::thread::available_parallelism().map_or(1, |n| n.get());
    eprintln!(
        "[arm-memory] {} tree rows, peak resident {:.2} GB with {workers} hardware threads (budget {:.2} GB)",
        asm.tree().len(),
        peak as f64 / 1073741824.0,
        PEAK_BUDGET as f64 / 1073741824.0
    );
    assert!(
        peak <= PEAK_BUDGET,
        "the arm build peaked at {:.2} GB resident, over its {:.2} GB budget",
        peak as f64 / 1073741824.0,
        PEAK_BUDGET as f64 / 1073741824.0
    );
}
