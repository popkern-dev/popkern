//! Cooperative cancellation shared by every long-running kernel
//! operation, plus the detach-on-cancel worker layer that makes the
//! "clean termination < 1 s" claim structural.
//!
//! Design: `docs/architecture/interactive-operation-control.md`
//! (PK-276), extending the shipped fillet mechanism of
//! `docs/architecture/fillet-cancellation-latency.md` (PK-190).
//!
//! # One token for the whole kernel
//!
//! [`CancellationToken`] is the single cancellation type for every
//! operation crate. `popkern_bool::cancellation`,
//! `popkern_blend::cancellation`, and `popkern_core::healing`
//! re-export it, so one token held by a UI Cancel button can fan into
//! a boolean, a fillet, a healing pass, and a re-tessellation at once.
//! (History: three structurally identical tokens existed before
//! PK-276, with divergent memory orderings; this definition settles on
//! `SeqCst` — the strongest ordering, unmeasurable at poll frequency.)
//!
//! # Detach-on-cancel
//!
//! Cooperative polls alone bound cancellation latency by the longest
//! poll-free stretch, which regresses silently as pipelines grow. The
//! detach layer makes the latency guarantee structural instead: the
//! operation runs on a worker thread over private clones of its
//! inputs, and the calling thread returns as soon as the token fires —
//! regardless of any phase's internal granularity. The abandoned
//! worker observes the same token at its next cooperative poll, drops
//! its clones, and exits; cooperative polls thus bound *wasted CPU*,
//! not user-facing latency.
//!
//! [`run_detach_on_cancel`] is that layer as a helper.
//! [`run_detach_on_cancel_with_events`] additionally lets the detached
//! worker deliver progress events to a caller-borrowed reporter: the
//! worker cannot hold the caller's `&mut` reporter across the thread
//! boundary, so it sends values through an [`EventSink`] and the
//! caller's poll loop forwards them in FIFO order.

use std::sync::Arc;
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc;
use std::time::Duration;

/// A cooperative cancellation flag shared across threads.
///
/// Operation pipelines poll [`Self::is_cancelled`] between phases (and
/// inside dominant loops). Calling [`Self::cancel`] from any thread
/// causes the next poll to return `true`, triggering a clean teardown —
/// no inconsistent state survives a cancellation. Cheap to clone
/// (`Arc<AtomicBool>`); all clones observe the same flag.
#[derive(Debug, Clone, Default)]
pub struct CancellationToken {
    flag: Arc<AtomicBool>,
}

impl CancellationToken {
    /// Create a fresh, non-cancelled token.
    pub fn new() -> Self {
        Self {
            flag: Arc::new(AtomicBool::new(false)),
        }
    }

    /// Signal cancellation. All clones of this token observe the
    /// request. Idempotent.
    pub fn cancel(&self) {
        self.flag.store(true, Ordering::SeqCst);
    }

    /// Returns `true` if any holder of this token has called
    /// [`Self::cancel`].
    pub fn is_cancelled(&self) -> bool {
        self.flag.load(Ordering::SeqCst)
    }
}

/// Outcome of a detach-on-cancel run.
#[derive(Debug)]
pub enum DetachOutcome<T> {
    /// The worker finished; its result was never discarded (completion
    /// wins a tie against a simultaneously fired token).
    Completed(T),
    /// The token fired first. The worker was abandoned: it keeps
    /// running until its next cooperative poll, then drops its private
    /// clones and exits. Callers map this to their `Cancelled` error;
    /// the caller-visible model is unchanged.
    Detached,
}

/// Sending half of the detached worker's progress channel.
///
/// Handed to the worker closure by
/// [`run_detach_on_cancel_with_events`]; the worker calls
/// [`Self::emit`] wherever it would call a reporter. Cheap to clone so
/// pipeline stages can carry their own handles.
#[derive(Debug, Clone)]
pub struct EventSink<E> {
    sender: mpsc::Sender<E>,
}

impl<E> EventSink<E> {
    /// Send one event to the caller's reporter. Never blocks; if the
    /// caller has detached (receiver dropped), the event is silently
    /// discarded — progress after cancellation is meaningless.
    pub fn emit(&self, event: E) {
        let _ = self.sender.send(event);
    }
}

/// Interval at which the calling thread polls the worker's completion
/// and the token. 2 ms keeps worst-case added latency far below the
/// 1 s product claim while costing nothing measurable.
const DETACH_POLL: Duration = Duration::from_millis(2);

/// Run `work` on a worker thread, returning early if `cancel` fires.
///
/// Contract (identical to the shipped PK-190 fillet mechanism):
///
/// - A pre-cancelled token returns [`DetachOutcome::Detached`] without
///   spawning a thread (and without running `work`).
/// - The worker owns pre-cloned inputs captured in `work`; the caller
///   polls completion/token every 2 ms.
/// - **Completion wins a tie** — a finished result is returned even if
///   the token also fired.
/// - A worker panic is re-raised on the calling thread
///   (`resume_unwind`): panic parity with a non-threaded call.
/// - On detach the worker is *not* joined; it must observe the same
///   token at its own cooperative polls to bound wasted CPU.
pub fn run_detach_on_cancel<T: Send + 'static>(
    cancel: &CancellationToken,
    work: impl FnOnce() -> T + Send + 'static,
) -> DetachOutcome<T> {
    run_detach_on_cancel_with_events(
        cancel,
        move |_sink: EventSink<Silent>| work(),
        |event| match event {},
    )
}

/// Uninhabited event type for [`run_detach_on_cancel`]: the no-events
/// case is the with-events case over a channel that cannot carry
/// anything.
#[derive(Debug, Clone, Copy)]
pub enum Silent {}

/// Like [`run_detach_on_cancel`], but the worker receives an
/// [`EventSink`] and every event it emits is forwarded to `on_event`
/// **on the calling thread**, in FIFO order. This is how a detached
/// operation reports progress to a caller-borrowed
/// `&mut dyn ProgressReporter` (see [`crate::progress`]).
///
/// On [`DetachOutcome::Completed`] the queue is fully drained before
/// returning, so no event is lost. On [`DetachOutcome::Detached`]
/// undelivered events are dropped with the receiver.
pub fn run_detach_on_cancel_with_events<T, E>(
    cancel: &CancellationToken,
    work: impl FnOnce(EventSink<E>) -> T + Send + 'static,
    mut on_event: impl FnMut(E),
) -> DetachOutcome<T>
where
    T: Send + 'static,
    E: Send + 'static,
{
    if cancel.is_cancelled() {
        return DetachOutcome::Detached;
    }
    let (sender, receiver) = mpsc::channel();
    let worker = std::thread::spawn(move || work(EventSink { sender }));
    loop {
        for event in receiver.try_iter() {
            on_event(event);
        }
        // Completion is checked before the token so a finished result
        // is never discarded on a tie.
        if worker.is_finished() {
            let result = match worker.join() {
                Ok(value) => value,
                Err(panic) => std::panic::resume_unwind(panic),
            };
            // The worker's sink is dropped; drain whatever it sent
            // between the last poll and completion.
            for event in receiver.try_iter() {
                on_event(event);
            }
            return DetachOutcome::Completed(result);
        }
        if cancel.is_cancelled() {
            return DetachOutcome::Detached;
        }
        std::thread::sleep(DETACH_POLL);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicUsize;

    #[test]
    fn fresh_token_is_not_cancelled() {
        assert!(!CancellationToken::new().is_cancelled());
        assert!(!CancellationToken::default().is_cancelled());
    }

    #[test]
    fn cancel_propagates_to_clones() {
        let token = CancellationToken::new();
        let clone = token.clone();
        token.cancel();
        assert!(clone.is_cancelled());
        assert!(token.is_cancelled());
    }

    #[test]
    fn detach_returns_completed_result() {
        let token = CancellationToken::new();
        match run_detach_on_cancel(&token, || 41 + 1) {
            DetachOutcome::Completed(v) => assert_eq!(v, 42),
            DetachOutcome::Detached => panic!("uncancelled work must complete"),
        }
    }

    #[test]
    fn pre_cancelled_token_detaches_without_running_work() {
        let token = CancellationToken::new();
        token.cancel();
        let ran = Arc::new(AtomicBool::new(false));
        let ran_flag = ran.clone();
        let outcome = run_detach_on_cancel(&token, move || {
            ran_flag.store(true, Ordering::SeqCst);
        });
        assert!(matches!(outcome, DetachOutcome::Detached));
        assert!(
            !ran.load(Ordering::SeqCst),
            "pre-cancelled detach must not spawn the worker"
        );
    }

    #[test]
    fn mid_run_cancel_detaches_promptly_and_worker_unwinds_on_its_own_poll() {
        let token = CancellationToken::new();
        let cancel_from_worker = token.clone();
        // The worker fires the token, then stays alive until the test
        // releases it — guaranteeing the caller's poll sees
        // cancelled-but-not-finished (a genuine mid-run detach, not a
        // completion race).
        let release = Arc::new(AtomicBool::new(false));
        let release_gate = release.clone();
        let worker_exited = Arc::new(AtomicBool::new(false));
        let exited_flag = worker_exited.clone();
        let outcome = run_detach_on_cancel(&token, move || {
            cancel_from_worker.cancel();
            while !release_gate.load(Ordering::SeqCst) {
                std::thread::sleep(Duration::from_millis(1));
            }
            exited_flag.store(true, Ordering::SeqCst);
        });
        assert!(matches!(outcome, DetachOutcome::Detached));
        // The abandoned worker keeps running, then exits on its own.
        release.store(true, Ordering::SeqCst);
        for _ in 0..500 {
            if worker_exited.load(Ordering::SeqCst) {
                return;
            }
            std::thread::sleep(Duration::from_millis(2));
        }
        panic!("abandoned worker never exited after release");
    }

    #[test]
    #[should_panic(expected = "worker panic parity")]
    fn worker_panic_is_reraised_on_the_caller_thread() {
        let token = CancellationToken::new();
        let _ = run_detach_on_cancel(&token, || panic!("worker panic parity"));
    }

    #[test]
    fn events_are_delivered_in_fifo_order_and_fully_drained() {
        let token = CancellationToken::new();
        let mut received = Vec::new();
        let outcome = run_detach_on_cancel_with_events(
            &token,
            |sink| {
                for i in 0..100usize {
                    sink.emit(i);
                }
                "done"
            },
            |event| received.push(event),
        );
        assert!(matches!(outcome, DetachOutcome::Completed("done")));
        assert_eq!(received, (0..100).collect::<Vec<_>>());
    }

    #[test]
    fn sink_emit_after_detach_is_silently_discarded() {
        let token = CancellationToken::new();
        let cancel_from_worker = token.clone();
        let release = Arc::new(AtomicBool::new(false));
        let release_gate = release.clone();
        let emitted_late = Arc::new(AtomicUsize::new(0));
        let late = emitted_late.clone();
        let outcome = run_detach_on_cancel_with_events(
            &token,
            move |sink| {
                cancel_from_worker.cancel();
                while !release_gate.load(Ordering::SeqCst) {
                    std::thread::sleep(Duration::from_millis(1));
                }
                // The caller has provably detached; emitting into the
                // dropped receiver must neither panic nor block.
                sink.emit(7usize);
                late.store(1, Ordering::SeqCst);
            },
            |_event: usize| {},
        );
        assert!(matches!(outcome, DetachOutcome::Detached));
        release.store(true, Ordering::SeqCst);
        for _ in 0..500 {
            if emitted_late.load(Ordering::SeqCst) == 1 {
                return;
            }
            std::thread::sleep(Duration::from_millis(2));
        }
        panic!("abandoned worker never reached its post-detach emit");
    }
}
