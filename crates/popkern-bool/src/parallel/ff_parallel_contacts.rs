//! Contact discovery in the face-pair pass (W3a, tangency design §4.8):
//! each unique surface pair the pass solves also searches its windows for
//! contacts. Read by nothing until `SourceKind::Contact` lands (W3b); with
//! the shadow log set, one `ContactDiscovery` line per boolean counts the
//! pairs searched, the contacts and touches within the reach, and the
//! search's time beside the solve's.

use popkern_core::geometry::intersection::contact_classify::ContactScale;
use popkern_core::geometry::intersection::contact_discovery::find_surface_contacts;
use popkern_core::tolerance::Tolerance;

use super::PairInputs;

/// With the shadow log set, search every unique pair for contacts and log
/// the counts and the search's cost (`solve_ms`: the pass's solve time).
pub(super) fn shadow_discover(unique: &[&PairInputs<'_>], tolerance: &Tolerance, solve_ms: f64) {
    if std::env::var_os("POPKERN_COMPLEX_SHADOW_LOG").is_none() || unique.is_empty() {
        return;
    }
    let clock = popkern_core::time::Stopwatch::start();
    let (mut contacts, mut within) = (0usize, 0usize);
    for inp in unique {
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
        "ContactDiscovery pairs {} contacts {contacts} within_tolerance {within} search_ms {:.2} solve_ms {solve_ms:.2}",
        unique.len(),
        clock.elapsed_ms().unwrap_or(0.0)
    ));
}
