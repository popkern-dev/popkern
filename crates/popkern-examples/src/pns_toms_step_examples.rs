//! STEP export showcase for the ACM TOMS Algorithm 1032 PnS fixtures.
//!
//! The fixture OBJs live in the optional, deliberately untracked ToMS
//! package (license — PK-315); the whole showcase gates on
//! `cfg(toms_package)` so fresh clones still build (PK-314). The real
//! implementation is `pns_toms_step_showcase.rs`.

#[cfg(toms_package)]
#[path = "pns_toms_step_showcase.rs"]
mod pns_toms_step_showcase;

#[cfg(toms_package)]
fn main() {
    pns_toms_step_showcase::run()
}

#[cfg(not(toms_package))]
fn main() {
    eprintln!(
        "pns_toms_step_examples: ACM ToMS 1032 package absent \
         (docs/code/ACM-ToMS-PnS-1032) — showcase unavailable on this \
         checkout; see PK-314/PK-315"
    );
}
