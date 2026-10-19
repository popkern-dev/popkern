//! Receiving-hole cuts on an existing body: tapped holes (real
//! modeled internal thread), clearance holes (with countersink /
//! counterbore), and external threads cut onto an existing shank.
//!
//! Every cut takes an arbitrary axis — `origin` on (or just above)
//! the entry surface, `direction` pointing INTO the material — so
//! angled fittings need no extra machinery, and a fastener built on
//! the same axis mates by construction: a tapped hole's entry cone
//! clears the head-side runout of a screw seated flush on it.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::ProfileShape;

use crate::FastenerError;
use crate::engine::{cut_thread, extrude_shape, revolved_tool, subtract};
use crate::entry_ring::{entry_ring_profile, probe_entry};
use crate::spec::ThreadSpec;
use crate::standards::ClearanceFit;
use crate::thread_form::{self, external_groove, internal_bore_radius, internal_ridge};

/// How deep a hole goes, measured along its axis from the origin.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum HoleDepth {
    /// A blind hole of the given depth.
    Blind(f64),
    /// Cut through the whole body (the extent is measured from the
    /// body's geometry).
    Through,
}

/// Options for [`tapped_hole`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TapOptions {
    /// Radial/flank fit clearance to the nominal external thread
    /// (0 = exact basic profile; 0.15–0.3 typical for FDM printing).
    pub clearance: f64,
    /// Cut a 45° lead-in chamfer at the entry (default on).
    pub entry_chamfer: bool,
}

impl Default for TapOptions {
    fn default() -> TapOptions {
        TapOptions {
            clearance: 0.0,
            entry_chamfer: true,
        }
    }
}

/// Options for [`clearance_hole`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ClearanceHoleOptions {
    /// Hole diameter series (default [`ClearanceFit::Normal`]).
    pub fit: ClearanceFit,
    /// Sink a 90° countersink for a flush countersunk head, sized
    /// from the head standard plus this extra radial room
    /// (`Some(0.0)` = nominal head diameter).
    pub countersink: Option<f64>,
    /// Sink a flat-bottomed counterbore for a socket head:
    /// `Some(extra)` sizes diameter and depth from the socket-head
    /// standard plus `extra` radial/axial room.
    pub counterbore: Option<f64>,
    /// Cut a 45° entry chamfer (default on; skipped automatically
    /// under a countersink or counterbore).
    pub entry_chamfer: bool,
}

impl Default for ClearanceHoleOptions {
    fn default() -> ClearanceHoleOptions {
        ClearanceHoleOptions {
            fit: ClearanceFit::Normal,
            countersink: None,
            counterbore: None,
            entry_chamfer: true,
        }
    }
}

/// How far above the origin the BORE extends, so sloped entry faces
/// are fully drilled (the region is outside the material on a flat
/// entry; a plain cylinder there costs nothing).
fn entry_headroom(spec: &ThreadSpec) -> f64 {
    0.5 * spec.major_diameter() + spec.lead()
}

/// How far above the origin the THREAD extends. Narrower than the
/// bore headroom deliberately: every extra lead of thread lengthens
/// the helical tool and its marched intersection curves (the march
/// budget bounds a cut to ~6 leads reliably today — the long-thread
/// campaign is the recorded follow-up). One lead covers entry faces
/// sloped up to ~lead/major.
fn thread_headroom(spec: &ThreadSpec) -> f64 {
    spec.lead()
}

/// Body vertex-hull extent along `direction` past `origin` (the far
/// side of a through cut, before margins).
fn through_reach(body: &BRep, origin: Point3, direction: Vec3) -> f64 {
    body.vertices
        .iter()
        .map(|(_, v)| (v.position - origin).dot(&direction))
        .fold(0.0_f64, f64::max)
}

/// The rigid placement carrying canonical +Z-axis tools onto the
/// hole axis.
/// The x axis seed of a hole's insertion frame along `direction`: the
/// world axis most orthogonal to it. A hole cut along an axis carries
/// its thread phase in this frame (z down the axis, x from this seed),
/// so a fastener placed in the same frame mates; placing one on a hole
/// needs exactly this rule (PK-897).
pub fn insertion_x_seed(direction: Vec3) -> Vec3 {
    let d = direction;
    if d.x.abs() <= d.y.abs() && d.x.abs() <= d.z.abs() {
        Vec3::new(1.0, 0.0, 0.0)
    } else if d.y.abs() <= d.z.abs() {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    }
}

fn axis_placement(
    origin: Point3,
    direction: Vec3,
    stage: &str,
) -> Result<Transform3, FastenerError> {
    let x_seed = insertion_x_seed(direction);
    Transform3::from_axis_placement(origin, direction, x_seed).map_err(|source| {
        FastenerError::Construction {
            stage: stage.to_string(),
            detail: format!("axis placement failed: {source}"),
        }
    })
}

/// Place a canonical-frame (+Z axis) tool onto the hole axis and
/// subtract it.
fn subtract_placed(
    body: &BRep,
    tool: &BRep,
    placement: &Transform3,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let placed = tool
        .transform(placement)
        .map_err(|source| FastenerError::Construction {
            stage: stage.to_string(),
            detail: format!("tool placement failed: {source}"),
        })?;
    subtract(body, &placed, stage, tolerance)
}

/// Cut a tapped (threaded) hole: a bore at the internal minor
/// diameter plus a real modeled internal thread with the requested
/// fit clearance, entered at `origin` along `direction` (pointing
/// into the material).
///
/// A blind tap bores to `depth` and threads to `depth − ¾ lead` (the
/// tap-relief shoulder); a through tap threads all the way out the
/// far side.
///
/// # Errors
///
/// [`FastenerError`] for spec/clearance problems or kernel declines
/// (the input body is unchanged on failure).
pub fn tapped_hole(
    body: &BRep,
    spec: &ThreadSpec,
    origin: Point3,
    direction: Vec3,
    depth: HoleDepth,
    options: &TapOptions,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let stage = "tapped hole";
    let ridge = internal_ridge(spec, options.clearance)?;
    let r_bore = internal_bore_radius(spec, options.clearance);
    let up = entry_headroom(spec);
    let (bore_depth, thread_end) = match depth {
        HoleDepth::Blind(d) => {
            if !(d.is_finite() && d > spec.lead()) {
                return Err(FastenerError::InvalidSpec {
                    field: "depth",
                    value: d,
                    requirement: "finite and > one thread lead",
                });
            }
            // Keep the whole tool inside the bore: its helix runs half a
            // lead past the nominal end and rounds up to a quarter turn
            // (`engine::plan_turns`), and the ridge reaches its axial
            // half width beyond that; a quarter pitch more keeps it off
            // the floor. Ending ¾ lead up, the tip met the floor 0.07
            // past it at every whole-pitch M3 depth, right on the floor
            // rim's seam vertex, and the subtraction declined (PK-895).
            let reach = 0.75 * spec.lead() + ridge.axial_half_width() + 0.25 * spec.pitch();
            (d, d - reach)
        }
        HoleDepth::Through => {
            let reach = through_reach(body, origin, direction);
            // Bore well past the far side; the thread only needs to
            // run off it by half a lead (its own runout adds the
            // rest).
            (reach + spec.major_diameter(), reach + 0.5 * spec.lead())
        }
    };
    let placement = axis_placement(origin, direction, stage)?;
    // Probe the entry face on the intact body: the bore and thread
    // cuts don't move its plane, and the pre-cut hull keeps the ring
    // lid conservative.
    let entry = probe_entry(body, origin, direction, up);
    let bore = extrude_shape(
        ProfileShape::Circle { radius: r_bore },
        -up,
        up + bore_depth,
        "tap bore",
        tolerance,
    )?;
    let mut result = subtract_placed(body, &bore, &placement, stage, tolerance)?;
    result = cut_thread(
        &result,
        spec,
        &ridge,
        thread_form::internal_phase_start(-thread_headroom(spec), spec),
        thread_end,
        None,
        Some(&placement),
        stage,
        tolerance,
    )?;
    if options.entry_chamfer {
        // A seated bearing-head screw's groove starts at its runout
        // (`head_side_full_major` below the head): every tap ridge
        // whose mating groove that clip removed meets the screw's
        // full-major crest, up to the runout plus a ridge half width.
        // The 45° cone reaches the major radius there, plus P/8 and the
        // clearance. Ending at P/8 + clearance, the tap's first ridge
        // met a flush M6 socket screw's full-major top turn (PK-849).
        let r_mouth = spec.major_diameter() / 2.0
            + thread_form::head_side_full_major(spec)?
            + external_groove(spec)?.axial_half_width()
            + spec.pitch() / 8.0
            + options.clearance;
        // Strictly inside the bore void — see the hex-nut mouth
        // chamfer: r_bore - embed is exactly the groove-base radius
        // and a coincident cutter cylinder poisons classification.
        let r_in = r_bore - 2.0 * thread_form::embed_margin(spec);
        let eps = 0.25 * spec.pitch();
        let profile = entry_ring_profile(r_in, r_mouth, eps, &entry, stage)?;
        let ring = revolved_tool(&profile, "tap entry chamfer", tolerance)?;
        result = subtract_placed(&result, &ring, &placement, stage, tolerance)?;
    }
    Ok(result)
}

/// Cut a plain clearance hole for a screw of `spec`, with optional
/// countersink (flush 90° head) or counterbore (sunk socket head).
///
/// # Errors
///
/// As [`tapped_hole`].
pub fn clearance_hole(
    body: &BRep,
    spec: &ThreadSpec,
    origin: Point3,
    direction: Vec3,
    depth: HoleDepth,
    options: &ClearanceHoleOptions,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let stage = "clearance hole";
    let hole_d = spec.iso_clearance_diameter(options.fit);
    let up = entry_headroom(spec);
    let bore_depth = match depth {
        HoleDepth::Blind(d) => {
            if !(d.is_finite() && d > 0.0) {
                return Err(FastenerError::InvalidSpec {
                    field: "depth",
                    value: d,
                    requirement: "finite and > 0",
                });
            }
            d
        }
        HoleDepth::Through => through_reach(body, origin, direction) + spec.major_diameter(),
    };
    let placement = axis_placement(origin, direction, stage)?;
    // Probe the entry face on the intact body (see `tapped_hole`).
    let entry = probe_entry(body, origin, direction, up);
    let bore = extrude_shape(
        ProfileShape::Circle {
            radius: hole_d / 2.0,
        },
        -up,
        up + bore_depth,
        "clearance bore",
        tolerance,
    )?;
    let mut result = subtract_placed(body, &bore, &placement, stage, tolerance)?;
    let r_hole = hole_d / 2.0;
    if let Some(extra) = options.countersink {
        let r_head = spec.iso_countersunk_head_diameter() / 2.0 + extra;
        let eps = 0.25 * spec.pitch();
        // 90° cone from the head radius at the surface down to the
        // hole wall; hypotenuse on the line r + z = r_head.
        let r_in = r_hole - eps;
        let profile = entry_ring_profile(r_in, r_head, eps, &entry, stage)?;
        let ring = revolved_tool(&profile, "countersink", tolerance)?;
        result = subtract_placed(&result, &ring, &placement, stage, tolerance)?;
    } else if let Some(extra) = options.counterbore {
        let dims = spec.iso_socket_head();
        let r_cb = dims.diameter / 2.0 + extra;
        let cb = extrude_shape(
            ProfileShape::Circle { radius: r_cb },
            -up,
            up + dims.height + extra,
            "counterbore",
            tolerance,
        )?;
        result = subtract_placed(&result, &cb, &placement, stage, tolerance)?;
    } else if options.entry_chamfer {
        let c = 0.1 * spec.major_diameter();
        let eps = 0.25 * spec.pitch();
        let profile = entry_ring_profile(r_hole - eps, r_hole + c, eps, &entry, stage)?;
        let ring = revolved_tool(&profile, "hole entry chamfer", tolerance)?;
        result = subtract_placed(&result, &ring, &placement, stage, tolerance)?;
    }
    Ok(result)
}

/// Cut an external thread onto an existing (cylindrical) shank
/// region: the groove is cut between `z_from` and `z_to` measured
/// along `direction` from `origin` (the shank axis).
///
/// # Errors
///
/// As [`tapped_hole`]. The cut does not verify the region is a
/// cylinder of the major diameter; it cuts the thread space wherever
/// material intersects the tool.
pub fn external_thread(
    body: &BRep,
    spec: &ThreadSpec,
    origin: Point3,
    direction: Vec3,
    z_from: f64,
    z_to: f64,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let stage = "external thread";
    let groove = external_groove(spec)?;
    if !(z_to.is_finite() && z_from.is_finite() && z_to > z_from) {
        return Err(FastenerError::InvalidSpec {
            field: "z_to - z_from",
            value: z_to - z_from,
            requirement: "finite and > 0",
        });
    }
    let placement = axis_placement(origin, direction, stage)?;
    cut_thread(
        body,
        spec,
        &groove,
        z_from,
        z_to,
        None,
        Some(&placement),
        stage,
        tolerance,
    )
}

#[cfg(test)]
#[path = "seated_tap_tests.rs"]
mod seated_tap_tests;
