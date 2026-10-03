//! Wrapped-ring pairing, stashing, and pre-dissolution for the
//! boundary-first band slit (split from `boundary_first_wrap_slit.rs`
//! under the 500-line structural budget).

use popkern_core::geometry::surface::SurfaceGeometry;

use super::BoundaryRing;
use super::boundary_first_wrap_slit::{collapse_stubs, net_winding, slit_band};

/// Pair and join a face's wrapped rings (moved out of
/// `boundary_first.rs` under its size ratchet): exactly two
/// |winding| = 1 rings slit into one chartable polygon; anything else
/// keeps the honest decline with the ring dump.
pub(super) fn join_wrapped_rings(
    face_handle: popkern_core::topology::FaceHandle,
    wrapped_rings: &[BoundaryRing],
    p: f64,
    surface: &SurfaceGeometry,
    orientation: bool,
    chord_height: f64,
) -> Option<Vec<BoundaryRing>> {
    // One ring winding once: a cap around a chart pole (PK-636).
    if let [ring] = wrapped_rings
        && let Some(cap) = super::boundary_first_pole_cap::close_through_pole(
            ring,
            p,
            surface,
            orientation,
            chord_height,
        )
    {
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!(
                "    [bf] face {face_handle:?} pole cap: {} pts",
                cap.uv.len()
            );
        }
        return Some(vec![cap]);
    }
    // Stored walk direction decides the SIGN; the slit orients each
    // piece explicitly, so only |winding| matters.
    let in_scope = !wrapped_rings.is_empty()
        && wrapped_rings.len() <= 4
        && wrapped_rings
            .iter()
            .all(|r| (1..=8).contains(&net_winding(r, p).abs()));
    if in_scope && let Some(joined) = slit_band(wrapped_rings, p, surface) {
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!(
                "    [bf] face {face_handle:?} annular-band slit: {} polygon(s), {} pts total",
                joined.len(),
                joined.iter().map(|r| r.uv.len()).sum::<usize>()
            );
        }
        return Some(joined);
    }
    if popkern_core::trace::trace_enabled("boundary-first") {
        let winds: Vec<i64> = wrapped_rings.iter().map(|r| net_winding(r, p)).collect();
        eprintln!(
            "    [bf] face {face_handle:?} u-wrap decline ({} wrapped ring(s), windings {winds:?})",
            wrapped_rings.len()
        );
    }
    for r in wrapped_rings {
        super::boundary_first_unwrap::dump_ring(r, "u");
    }
    None
}

/// Dissolve a wrap-failed ring's slits, pinches, and stubs, then
/// stash the surviving lobes for pairing. The ring arrives
/// consecutively unwrapped, so slit dissolution is UNAMBIGUOUS (the
/// periodic scoping in `dissolve_slits` guards raw-branch seam
/// coincidence, which the unwrap already resolved), and pinches must
/// dissolve here because the original decline returned before the
/// pinch stage — a dissolved-seam stub left in place collapses the
/// slit polygon's coverage with zero-width T-points (measured: an
/// out-and-back spike at u = 5π/2 on the m5 tip band's bottom ring).
pub(super) fn stash_wrapped_ring(ring: BoundaryRing, wrapped_rings: &mut Vec<BoundaryRing>) {
    for piece in super::boundary_first_slit::dissolve_slits(ring, (None, None)) {
        for mut lobe in super::boundary_first_slit::dissolve_pinches(piece) {
            collapse_stubs(&mut lobe);
            if lobe.uv.len() >= 3 {
                wrapped_rings.push(lobe);
            }
        }
    }
}
