//! Contact discovery in the face-pair pass (W3a, tangency design §4.8):
//! each unique surface pair the pass solves also searches its windows for
//! contacts where it found no curve. Read by nothing until `SourceKind::Contact` lands (W3b); with
//! the shadow log set, one `ContactDiscovery` line per boolean counts the
//! pairs searched, the contacts and touches within the reach, and the
//! search's time beside the solve's.

use popkern_core::geometry::intersection::contact_classify::ContactScale;
use popkern_core::geometry::intersection::contact_discovery::find_surface_contacts;
use popkern_core::tolerance::Tolerance;

use super::PairInputs;

/// With the shadow log set, search every unique pair for contacts and log
/// the counts and the search's cost (`solve_ms`: the pass's solve time).
pub(super) fn shadow_discover(
    unique: &[&PairInputs<'_>],
    solved: &[Result<
        popkern_core::geometry::intersection::SurfaceSurfaceIntersectionResult,
        String,
    >],
    tolerance: &Tolerance,
    solve_ms: f64,
) {
    if std::env::var_os("POPKERN_COMPLEX_SHADOW_LOG").is_none() || unique.is_empty() {
        return;
    }
    let clock = popkern_core::time::Stopwatch::start();
    let (mut contacts, mut within) = (0usize, 0usize);
    // Only the pairs whose solve found nothing: a pair meeting in a curve
    // has its meetings as records, and a touch is a pair whose boxes
    // overlap with no curve between them (the design's skip, §4.8).
    let mut searched = 0usize;
    for (inp, solve) in unique.iter().zip(solved) {
        if solve.as_ref().map_or(true, |r| !r.intersections.is_empty()) {
            continue;
        }
        searched += 1;
        // The faces' size: the larger window's diagonal in space.
        let diagonal = |s: &popkern_core::geometry::surface::SurfaceGeometry,
                        ((u0, u1), (v0, v1)): ((f64, f64), (f64, f64))| {
            use popkern_core::geometry::surface::SurfaceEval as _;
            s.point_at(u0, v0).distance_to(&s.point_at(u1, v1))
        };
        let extent = diagonal(inp.surf_a, inp.r1)
            .max(diagonal(inp.surf_b, inp.r2))
            .max(tolerance.linear());
        let search = find_surface_contacts(
            inp.surf_a,
            inp.r1,
            inp.surf_b,
            inp.r2,
            ContactScale {
                linear: tolerance.linear(),
                extent,
            },
        );
        contacts += search.contacts.len();
        within += search
            .contacts
            .iter()
            .filter(|c| c.gap.abs() <= tolerance.linear())
            .count();
    }
    crate::complex_shadow::append_log_line(&format!(
        "ContactDiscovery pairs {} searched {searched} contacts {contacts} within_tolerance {within} search_ms {:.2} solve_ms {solve_ms:.2}",
        unique.len(),
        clock.elapsed_ms().unwrap_or(0.0)
    ));
}
