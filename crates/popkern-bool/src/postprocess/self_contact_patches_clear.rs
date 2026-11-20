//! Whether a NURBS face's window lies provably farther than a pad from
//! an analytic surface: its Bézier patches, halved where near, each
//! proved clear by its control box (the patch lies in it, weights
//! positive) against the surface in closed form
//! (`ff_parallel_reach::box_clear_of_surface`). The face-pair pass's
//! contact search skips a pair so proved (tangency design §4.8): a
//! threaded part's helical walls beside a bore took 0.74 s a pair to
//! search for contacts that cannot be there.

use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;

use super::patches::{Patch, window_patches};

/// Halvings of a patch before it is left undecided.
const DEPTH: usize = 12;

/// Patches tested before the proof gives up.
const BUDGET: usize = 20_000;

/// Whether `s` over `window` lies farther than `pad` from the whole of
/// `other`: true only when proved.
pub(crate) fn window_clear_of_surface(
    s: &NurbsSurface,
    window: &((f64, f64), (f64, f64)),
    other: &SurfaceGeometry,
    pad: f64,
) -> bool {
    let Some(patches) = window_patches(s, window) else {
        return false;
    };
    let mut stack: Vec<(Patch, usize)> = patches.into_iter().map(|p| (p, 0)).collect();
    let mut tested = 0usize;
    while let Some((patch, depth)) = stack.pop() {
        tested += 1;
        if tested > BUDGET {
            return false;
        }
        if crate::parallel::ff_parallel::ff_parallel_reach::box_clear_of_surface(
            &patch.bx, other, pad,
        ) {
            continue;
        }
        if depth == DEPTH {
            return false;
        }
        let (along_u, along_v) = patch.extents();
        let Some((lo, hi)) = patch.halves(along_u >= along_v) else {
            return false;
        };
        stack.push((lo, depth + 1));
        stack.push((hi, depth + 1));
    }
    true
}
