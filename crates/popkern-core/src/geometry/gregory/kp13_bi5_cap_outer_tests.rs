//! Tests for KP13 I.14.b outer-corner-jet transfer.

#![allow(clippy::needless_range_loop)]

use super::super::kp13_bi5_boundary::TensorBorder;
use super::super::kp13_bi5_cap::{Bi5Cap, Bi5CapPiece, CapPieceId};
use super::super::kp13_bi5_chi::ChiCorner;
use super::super::kp13_bi5_jet::TaylorJet2;
use super::*;

#[path = "kp13_bi5_cap_outer_tests_corner_jets.rs"]
mod tests_corner_jets;
#[path = "kp13_bi5_cap_outer_tests_corner_transfer.rs"]
mod tests_corner_transfer;
#[path = "kp13_bi5_cap_outer_tests_outer_layers.rs"]
mod tests_outer_layers;
#[path = "kp13_bi5_cap_outer_tests_support.rs"]
mod tests_support;
