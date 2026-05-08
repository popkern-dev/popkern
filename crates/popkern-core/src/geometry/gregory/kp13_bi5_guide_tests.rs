//! Tests for the bi-5 guide construction (KP13 § 3).

#![allow(clippy::needless_range_loop)]

use super::super::bernstein_tri5::{
    TRI5_CP_COUNT, bernstein_tri5_all, bernstein_tri5_all_extended,
};
use super::super::kp13_bi5_free_params::{Bi5FreeParams, apply_free_params};
use super::super::kp13_bi5_gtri5::{IDX_005, IDX_050, IDX_500};
use super::*;

#[path = "kp13_bi5_guide_tests_corner_identities.rs"]
mod tests_corner_identities;
#[path = "kp13_bi5_guide_tests_degree_raise.rs"]
mod tests_degree_raise;
#[path = "kp13_bi5_guide_tests_determinism.rs"]
mod tests_determinism;
#[path = "kp13_bi5_guide_tests_eq8_bilinear.rs"]
mod tests_eq8_bilinear;
#[path = "kp13_bi5_guide_tests_eq8_cc.rs"]
mod tests_eq8_cc;
#[path = "kp13_bi5_guide_tests_from_gtri5.rs"]
mod tests_from_gtri5;
#[path = "kp13_bi5_guide_tests_support.rs"]
mod tests_support;
#[path = "kp13_bi5_guide_tests_validation_eval.rs"]
mod tests_validation_eval;
