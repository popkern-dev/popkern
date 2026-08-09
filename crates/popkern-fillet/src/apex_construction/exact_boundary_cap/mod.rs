//! Exact-boundary apex cap primitives.

pub(crate) mod boundary;
pub(crate) mod cap_quality_postcondition;
pub(crate) mod cap_spoke_commit;
pub(crate) mod commit;
pub(crate) mod commit_diagnostics;
pub(crate) mod corner;
mod corner_pole_search_hints;
mod cross_spoke_g1;
pub(crate) mod facade;
mod facade_pole_search;
pub(crate) mod facts;
#[cfg(test)]
#[path = "informational_channel_tests.rs"]
mod informational_channel_tests;
pub(crate) mod multipatch;
pub(crate) mod pcurves;
pub(crate) mod performance;
mod pns_candidate;
pub(crate) mod ribbon;
pub(crate) mod ribbon_projection;
pub(crate) mod sector;
pub(crate) mod solve;
pub(crate) mod support;
pub(crate) mod support_subdivision;
pub(crate) mod topology;
pub(crate) mod validation;

#[cfg(test)]
#[path = "adversarial_tests.rs"]
mod adversarial_tests;

#[cfg(test)]
#[path = "boundary_tests.rs"]
mod boundary_tests;

#[cfg(test)]
#[path = "cap_degeneracy_matrix_tests.rs"]
mod cap_degeneracy_matrix_tests;

#[cfg(test)]
#[path = "corner_tests.rs"]
mod corner_tests;

#[cfg(test)]
#[path = "facade_tests.rs"]
mod facade_tests;

#[cfg(test)]
#[path = "commit_tests.rs"]
mod commit_tests;

#[cfg(test)]
#[path = "ribbon_tests.rs"]
mod ribbon_tests;

#[cfg(test)]
#[path = "n_sector_tests.rs"]
mod n_sector_tests;

#[cfg(test)]
#[path = "performance_determinism_tests.rs"]
mod performance_determinism_tests;

#[cfg(test)]
#[path = "pcurves_tests.rs"]
mod pcurves_tests;

#[cfg(test)]
#[path = "sector_tests.rs"]
mod sector_tests;

#[cfg(test)]
#[path = "solve_tests.rs"]
mod solve_tests;

#[cfg(test)]
#[path = "support_candidate_tests.rs"]
mod support_candidate_tests;

#[cfg(test)]
#[path = "topology_tests.rs"]
mod topology_tests;

#[cfg(test)]
#[path = "validation_tests.rs"]
mod validation_tests;
