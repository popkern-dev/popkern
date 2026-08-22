//! Visible skip notice for the optional ACM ToMS 1032 package (see
//! PK-314/PK-315): the PnS fixture suites gate on `cfg(toms_package)`
//! and vanish when the package is absent — this test keeps the skip
//! honest in the run output instead of silently shrinking the suite.

#[cfg(not(toms_package))]
#[test]
fn toms_package_absent_notice() {
    eprintln!(
        "[toms] ACM ToMS 1032 package absent — pns_toms_step_roundtrip and \
         pns_toms_fixture_matrix suites skipped (see PK-314/PK-315)"
    );
}
