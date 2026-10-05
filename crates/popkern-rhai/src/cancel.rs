//! Host cancellation of a running script (PK-665): each thread that
//! runs scripts owns a flag; a host holding the flag (a server's socket
//! thread, while the script runs on its session thread) sets it, and
//! the run stops at the next interpreter progress tick with the
//! `Cancelled` class. A kernel operation already in progress finishes
//! first: the flag is checked between script steps, not inside kernel
//! algorithms.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};

thread_local! {
    static FLAG: Arc<AtomicBool> = Arc::new(AtomicBool::new(false));
}

/// The calling thread's cancellation flag (share it with whoever may
/// cancel this thread's runs).
pub fn cancel_handle() -> Arc<AtomicBool> {
    FLAG.with(Arc::clone)
}

/// Clear the calling thread's flag (a run is starting).
pub(crate) fn clear() {
    FLAG.with(|f| f.store(false, Ordering::SeqCst));
}

/// The message a cancelled run terminates with.
pub(crate) const CANCELLED: &str = "cancelled by the host";
