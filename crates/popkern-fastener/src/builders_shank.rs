//! The shank-plus-thread stage of headed screws, with the modelled or
//! cosmetic thread (PK-847).
//!
//! A cosmetic thread is the common CAD convention for threads in
//! assemblies: the threaded length is a plain cylinder at the external
//! thread's root diameter d3 = d − 2·(17/24)H (ISO 724, the
//! drawing convention), and the feature keeps its full thread spec (it
//! still names the M5 thread; it just does not cut the helix). A screw
//! then clears a tap-drill bore (M5: 4.019 in 4.2, M4: 3.141 in 3.3)
//! and a real tapped hole's crests (M6: 4.773 under D1 = 4.917)
//! instead of overlapping the bore by the thread depth or riding the
//! crests tangentially (the basic profile's root, D1, is the internal
//! crest itself, and the boolean declined that contact on the arm's
//! gripper screws) — and it builds without the helical cut: the
//! six-axis arm's 126 screws reported as interfering by design.

use popkern_core::{tolerance::Tolerance, topology::BRep};
use popkern_sketch::ProfileShape;

use crate::FastenerError;
use crate::builders::tip_chamfer;
use crate::engine::{cut_thread, extrude_shape, union};
use crate::spec::ThreadSpec;
use crate::thread_form::external_groove;

/// Shared shank-plus-thread stage for headed screws: shank from
/// `overlap_from` (buried in the head) to `length`, tip chamfer,
/// external thread over the last `thread_length` (clamped clear of
/// the head at `head_clear_z`). A cosmetic thread makes that length a
/// root-diameter cylinder (no tip chamfer: the step to the root
/// diameter is the thread's lead-in).
#[allow(clippy::too_many_arguments)]
pub(crate) fn headed_screw(
    head: BRep,
    spec: &ThreadSpec,
    length: f64,
    thread_length: Option<f64>,
    overlap_from: f64,
    head_clear_z: f64,
    chamfer_tip: bool,
    cosmetic: bool,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let requested = thread_length.unwrap_or(length);
    let z_start = (length - requested).max(head_clear_z);
    if cosmetic {
        return cosmetic_shank(head, spec, length, overlap_from, z_start, tolerance);
    }
    let groove = external_groove(spec)?;
    let shank = extrude_shape(
        ProfileShape::Circle {
            radius: spec.major_diameter() / 2.0,
        },
        overlap_from,
        length - overlap_from,
        "screw shank",
        tolerance,
    )?;
    let mut body = union(&head, &shank, "screw head ∪ shank", tolerance)?;
    if chamfer_tip {
        body = tip_chamfer(&body, spec, length, tolerance)?;
    }
    cut_thread(
        &body,
        spec,
        &groove,
        z_start,
        length,
        Some(head_clear_z + crate::thread_form::head_side_full_major(spec)?),
        None,
        "external thread",
        tolerance,
    )
}

/// The external root diameter d3 = d − 2·(17/24)H a cosmetic thread
/// is drawn at.
fn cosmetic_diameter(spec: &ThreadSpec) -> f64 {
    spec.major_diameter() - 2.0 * (17.0 / 24.0) * spec.fundamental_height()
}

/// The cosmetic shank: major diameter from the head to `z_start`, the
/// root diameter from there to the tip.
fn cosmetic_shank(
    head: BRep,
    spec: &ThreadSpec,
    length: f64,
    overlap_from: f64,
    z_start: f64,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let root = extrude_shape(
        ProfileShape::Circle {
            radius: cosmetic_diameter(spec) / 2.0,
        },
        overlap_from,
        length - overlap_from,
        "cosmetic thread",
        tolerance,
    )?;
    let mut body = union(&head, &root, "screw head ∪ cosmetic thread", tolerance)?;
    // An unthreaded run under the head keeps the major diameter.
    if z_start > overlap_from + 1e-9 * length.abs().max(1.0) {
        let plain = extrude_shape(
            ProfileShape::Circle {
                radius: spec.major_diameter() / 2.0,
            },
            overlap_from,
            z_start - overlap_from,
            "screw shank",
            tolerance,
        )?;
        body = union(&body, &plain, "screw ∪ plain shank", tolerance)?;
    }
    Ok(body)
}

/// A cosmetic threaded rod: the root-diameter cylinder from 0 to
/// `length`.
pub(crate) fn cosmetic_rod(
    spec: &ThreadSpec,
    length: f64,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    extrude_shape(
        ProfileShape::Circle {
            radius: cosmetic_diameter(spec) / 2.0,
        },
        0.0,
        length,
        "cosmetic threaded rod",
        tolerance,
    )
}
