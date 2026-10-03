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
    transposed: bool,
) -> Option<Vec<BoundaryRing>> {
    // One ring winding once: a cap around a chart pole (PK-636).
    if !transposed
        && let [ring] = wrapped_rings
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
    let sided = band_side(wrapped_rings, p, surface, orientation, transposed);
    let wrapped_rings: &[BoundaryRing] = sided.as_deref().unwrap_or(wrapped_rings);
    if in_scope && let Some(joined) = slit_band(wrapped_rings, p, surface, transposed, chord_height)
    {
        // Rings that wrapped in v were slit as u bands; back to (u, v).
        let joined: Vec<BoundaryRing> = if transposed {
            joined.into_iter().map(transposed_ring).collect()
        } else {
            joined
        };
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

/// Note a face's wrapped ring as `(period, transposed)`; a face whose
/// rings wrap in BOTH directions (u and v) declines (`None`).
pub(super) fn record_wrap(slot: &mut Option<(f64, bool)>, wrap: (f64, bool)) -> Option<()> {
    match slot.replace(wrap) {
        Some((_, t)) if t != wrap.1 => None,
        _ => Some(()),
    }
}

/// Stash a ring that wraps in v with u and v exchanged (positions and
/// spans unchanged): a band whose rings wrap in v (a torus elbow's tube
/// circles) is slit by the u-band machinery in the transposed chart.
pub(super) fn stash_transposed(ring: BoundaryRing, wrapped_rings: &mut Vec<BoundaryRing>) {
    stash_wrapped_ring(transposed_ring(ring), wrapped_rings);
}

fn transposed_ring(mut ring: BoundaryRing) -> BoundaryRing {
    for q in &mut ring.uv {
        *q = (q.1, q.0);
    }
    ring
}

/// On a surface periodic ACROSS the band too (a torus), two wrapped
/// rings bound two bands, both chartable; the face is the one its flag
/// names. The kernel flag says the loop runs counter-clockwise in the
/// chart, so a ring walking toward +first coordinate has the face on
/// its +second side (handedness flips in a transposed chart). When the
/// lower ring does not, the face is the complementary band: the lower
/// ring moves one period up. `None` keeps the rings as they are
/// (measured: jack_in_the_box's S-bend tubes meshed the outer half of
/// one torus, enclosing 0.14 of 58.2).
fn band_side(
    rings: &[BoundaryRing],
    p: f64,
    surface: &SurfaceGeometry,
    orientation: bool,
    transposed: bool,
) -> Option<Vec<BoundaryRing>> {
    let [a, b] = rings else { return None };
    let (period_u, period_v) = surface.periodicity();
    let cross = if transposed { period_u } else { period_v }?;
    let mean = |r: &BoundaryRing| r.uv.iter().map(|q| q.1).sum::<f64>() / r.uv.len().max(1) as f64;
    let (low, high) = if mean(a) <= mean(b) { (a, b) } else { (b, a) };
    let face_above_low = (orientation != transposed) == (net_winding(low, p) > 0);
    if face_above_low {
        return None;
    }
    let mut lifted = BoundaryRing {
        uv: low.uv.clone(),
        pos: low.pos.clone(),
        sources: low.sources.clone(),
    };
    for q in &mut lifted.uv {
        q.1 += cross;
    }
    let other = BoundaryRing {
        uv: high.uv.clone(),
        pos: high.pos.clone(),
        sources: high.sources.clone(),
    };
    Some(vec![other, lifted])
}
