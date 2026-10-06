//! Fastener body builders: threaded rod, hex bolt, socket-cap screw,
//! countersunk screw, hex nut — parameter structs plus one build
//! function each, returning validated breps in the insertion frame
//! (axis +Z, origin at the bearing plane, threaded end pointing +Z;
//! heads sit at z ≤ 0, a countersunk head sinks into z ≥ 0).
//!
//! Head proportions default from [`standards`](crate::standards) and
//! can be overridden per parameter struct.

use popkern_core::{tolerance::Tolerance, topology::BRep};
use popkern_math::PortableMath;
use popkern_sketch::ProfileShape;

use crate::FastenerError;
use crate::engine::{
    chamfer_circular_rims_at, cut_thread, extrude_shape, revolved_tool, subtract, union,
};
use crate::spec::ThreadSpec;
use crate::thread_form::{self, external_groove, internal_bore_radius, internal_ridge};

/// Hexagon corner points for a width across flats, one flat facing +Y.
pub(crate) fn hexagon(across_flats: f64) -> ProfileShape<f64> {
    let circumradius = across_flats / 3.0_f64.sqrt();
    ProfileShape::Polygon {
        points: (0..6)
            .map(|i| {
                let a = std::f64::consts::FRAC_PI_6 + f64::from(i) * std::f64::consts::FRAC_PI_3;
                [circumradius * a.pcos(), circumradius * a.psin()]
            })
            .collect(),
    }
}

fn positive(stage: &str, field: &'static str, value: f64) -> Result<(), FastenerError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(FastenerError::InvalidSpec {
            field,
            value,
            requirement: "finite and > 0",
        })
        .map_err(|e| stage_context(stage, e))
    }
}

/// Spec-level errors already carry their own structure; the stage
/// context only matters for construction reporting, so pass them
/// through unchanged.
fn stage_context(_stage: &str, e: FastenerError) -> FastenerError {
    e
}

/// Chamfer the tip of a shank at `z` to just below the minor
/// diameter (thread lead-in; 45°).
pub(crate) fn tip_chamfer(
    body: &BRep,
    spec: &ThreadSpec,
    z: f64,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let distance = spec.thread_depth() + spec.pitch() / 8.0;
    chamfer_circular_rims_at(body, z, distance, "tip chamfer", tolerance)
}

/// Diameter of the hex corner-chamfer circle on the exposed face, as
/// a fraction of the width across flats.
///
/// ISO 4014/4017/4032 draw the 15°–30° top chamfer with its circle
/// at roughly the across-flats width. Exactly `1.0` would make the
/// cone TANGENT to every flat along the face edge (a graze contact
/// the boolean must classify at six points); `0.95` keeps the cone
/// transversal to the face plane (a circle), to each flat (an arc
/// crowning the flat) and to each corner edge, which is also the
/// conventional modelled look of a chamfered hex.
const HEX_CHAMFER_CIRCLE_RATIO: f64 = 0.95;

/// Angle between the chamfer cone and the exposed face (30°).
const HEX_CHAMFER_ANGLE_DEG: f64 = 30.0;

/// The 30° corner-chamfer ring for a hex prism: a revolved cutter
/// whose cone meets the exposed face on the chamfer circle
/// ([`HEX_CHAMFER_CIRCLE_RATIO`]) and rises at 30° INTO the material,
/// trimming the hexagon corners the way a facing cut does.
/// `outward = -1.0` for a face looking down −Z, `+1.0` for +Z.
///
/// The cutter's base sits `eps` outside the face so no tool face is
/// coplanar with the prism's. The base point is slid down the SAME
/// cone line (radially inward by `eps / tan 30°`), so the cone still
/// passes through the chamfer circle at the face: offsetting the
/// base only axially (the PK-669 defect) tilted the cone to 34.1°
/// and moved its face circle to 0.5738·s against the 0.5774·s
/// circumradius, leaving six near-graze slivers 0.0024·s high in
/// place of the chamfer (0.013 mm on M3, where the cellular boolean
/// declined on them).
pub(crate) fn hex_corner_ring(
    across_flats: f64,
    z_face: f64,
    outward: f64,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let slope = HEX_CHAMFER_ANGLE_DEG.to_radians().ptan();
    let r_face = HEX_CHAMFER_CIRCLE_RATIO * across_flats / 2.0;
    // Just beyond the corners (circumradius s/√3): the cutter reaches
    // only 0.076·s into the material, so thin heads/nuts keep their
    // opposite face clear of it.
    let r_out = 1.05 * across_flats / 3.0_f64.sqrt();
    let eps = across_flats * 0.05;
    let r_base = r_face - eps / slope;
    let rise = (r_out - r_face) * slope;
    // Base edge just OUTSIDE the exposed face, cone rising INTO the
    // material through the chamfer circle.
    let points = [
        (r_base, z_face + outward * eps),
        (r_out, z_face + outward * eps),
        (r_out, z_face - outward * rise),
    ];
    revolved_tool(&points, "hex corner chamfer", tolerance)
}

/// A straight threaded rod: full-length external thread, both ends
/// free-running.
#[derive(Debug, Clone, PartialEq)]
pub struct RodParams {
    /// The thread.
    pub spec: ThreadSpec,
    /// Rod length (z = 0..length).
    pub length: f64,
    /// 45° lead-in chamfer on the +Z tip (default on).
    pub tip_chamfer: bool,
}

/// Build a threaded rod.
///
/// # Errors
///
/// [`FastenerError`] for spec problems or kernel sweep / boolean /
/// blend declines, with the failing stage named.
pub fn threaded_rod(params: &RodParams, tolerance: Tolerance) -> Result<BRep, FastenerError> {
    let stage = "threaded rod";
    params.spec.validate()?;
    positive(stage, "length", params.length)?;
    let groove = external_groove(&params.spec)?;
    let mut body = extrude_shape(
        ProfileShape::Circle {
            radius: params.spec.major_diameter() / 2.0,
        },
        0.0,
        params.length,
        stage,
        tolerance,
    )?;
    if params.tip_chamfer {
        body = tip_chamfer(&body, &params.spec, params.length, tolerance)?;
    }
    cut_thread(
        &body,
        &params.spec,
        &groove,
        0.0,
        params.length,
        None,
        None,
        "external thread",
        tolerance,
    )
}

/// Shared shank-plus-thread stage for headed screws: shank from
/// `overlap_from` (buried in the head) to `length`, tip chamfer,
/// external thread over the last `thread_length` (clamped clear of
/// the head at `head_clear_z`).
#[allow(clippy::too_many_arguments)]
fn headed_screw(
    head: BRep,
    spec: &ThreadSpec,
    length: f64,
    thread_length: Option<f64>,
    overlap_from: f64,
    head_clear_z: f64,
    chamfer_tip: bool,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
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
    let lead = spec.lead();
    let requested = thread_length.unwrap_or(length);
    let z_start = (length - requested).max(head_clear_z);
    cut_thread(
        &body,
        spec,
        &groove,
        z_start,
        length,
        Some(head_clear_z + 0.1 * lead),
        None,
        "external thread",
        tolerance,
    )
}

/// A hex-head bolt (ISO 4017 proportions by default).
#[derive(Debug, Clone, PartialEq)]
pub struct HexBoltParams {
    /// The thread.
    pub spec: ThreadSpec,
    /// Under-head length (head at z ≤ 0, thread tip at z = length).
    pub length: f64,
    /// Thread only the last `thread_length` of the shank (None =
    /// full thread up to the head runout).
    pub thread_length: Option<f64>,
    /// Override the hex width across flats.
    pub across_flats: Option<f64>,
    /// Override the head height.
    pub head_height: Option<f64>,
    /// 30° hex corner chamfer on the exposed face (default on).
    pub head_chamfer: bool,
    /// 45° tip lead-in chamfer (default on).
    pub tip_chamfer: bool,
}

/// Build a hex-head bolt.
///
/// # Errors
///
/// As [`threaded_rod`].
pub fn hex_bolt(params: &HexBoltParams, tolerance: Tolerance) -> Result<BRep, FastenerError> {
    let stage = "hex bolt";
    params.spec.validate()?;
    positive(stage, "length", params.length)?;
    let dims = params.spec.iso_hex_head();
    let s = params.across_flats.unwrap_or(dims.across_flats);
    let k = params.head_height.unwrap_or(dims.height);
    positive(stage, "across_flats", s)?;
    positive(stage, "head_height", k)?;
    let mut head = extrude_shape(hexagon(s), -k, k, "hex head", tolerance)?;
    if params.head_chamfer {
        let ring = hex_corner_ring(s, -k, -1.0, tolerance)?;
        head = subtract(&head, &ring, "hex head chamfer", tolerance)?;
    }
    headed_screw(
        head,
        &params.spec,
        params.length,
        params.thread_length,
        -k / 2.0,
        0.0,
        params.tip_chamfer,
        tolerance,
    )
}

/// A socket-head cap screw (ISO 4762 proportions by default).
#[derive(Debug, Clone, PartialEq)]
pub struct SocketCapParams {
    /// The thread.
    pub spec: ThreadSpec,
    /// Under-head length.
    pub length: f64,
    /// Thread only the last `thread_length` of the shank.
    pub thread_length: Option<f64>,
    /// 45° tip lead-in chamfer (default on).
    pub tip_chamfer: bool,
}

/// Build a socket-head cap screw.
///
/// # Errors
///
/// As [`threaded_rod`].
pub fn socket_cap_screw(
    params: &SocketCapParams,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let stage = "socket cap screw";
    params.spec.validate()?;
    positive(stage, "length", params.length)?;
    let dims = params.spec.iso_socket_head();
    let head_cyl = extrude_shape(
        ProfileShape::Circle {
            radius: dims.diameter / 2.0,
        },
        -dims.height,
        dims.height,
        "socket head",
        tolerance,
    )?;
    // Hex socket recess from the exposed face, 0.55·k deep.
    let socket_depth = 0.55 * dims.height;
    let recess = extrude_shape(
        hexagon(dims.socket_across_flats),
        -dims.height - 0.1 * dims.height,
        socket_depth + 0.1 * dims.height,
        "hex socket recess",
        tolerance,
    )?;
    let head = subtract(&head_cyl, &recess, "hex socket recess", tolerance)?;
    // Shank end midway through the web under the recess floor: at k/2
    // it refilled the recess to k/2 deep and, past M1.6's socket
    // corners (0.866 vs shank r 0.8), left six floor slivers (PK-671).
    let web = dims.height - socket_depth;
    headed_screw(
        head,
        &params.spec,
        params.length,
        params.thread_length,
        -web / 2.0,
        0.0,
        params.tip_chamfer,
        tolerance,
    )
}

/// A countersunk (flat-head) screw: 90° head sinking into the
/// material, top face flush at z = 0; `length` is overall (flush
/// plane to tip).
#[derive(Debug, Clone, PartialEq)]
pub struct CountersunkParams {
    /// The thread.
    pub spec: ThreadSpec,
    /// Overall length (flush plane to tip).
    pub length: f64,
    /// Thread only the last `thread_length` of the shank.
    pub thread_length: Option<f64>,
    /// 45° tip lead-in chamfer (default on).
    pub tip_chamfer: bool,
}

/// Build a countersunk screw.
///
/// # Errors
///
/// As [`threaded_rod`].
pub fn countersunk_screw(
    params: &CountersunkParams,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let stage = "countersunk screw";
    params.spec.validate()?;
    positive(stage, "length", params.length)?;
    let r_head = params.spec.iso_countersunk_head_diameter() / 2.0;
    let r_shank = params.spec.major_diameter() / 2.0;
    let head_depth = r_head - r_shank; // 90° cone: 1:1
    if params.length <= head_depth + params.spec.lead() {
        return Err(FastenerError::InvalidSpec {
            field: "length",
            value: params.length,
            requirement: "longer than the head depth plus one thread lead",
        });
    }
    // Head: flat top disc at z = 0 tapering at 90° to the shank.
    let head = revolved_tool(
        &[
            (0.0, 0.0),
            (r_head, 0.0),
            (r_shank, head_depth),
            (0.0, head_depth),
        ],
        "countersunk head",
        tolerance,
    )?;
    headed_screw(
        head,
        &params.spec,
        params.length,
        params.thread_length,
        head_depth / 2.0,
        head_depth,
        params.tip_chamfer,
        tolerance,
    )
}

/// A hex nut (ISO 4032 proportions by default) with a real modeled
/// internal thread, fit clearance, and lead-in chamfers.
#[derive(Debug, Clone, PartialEq)]
pub struct NutParams {
    /// The thread.
    pub spec: ThreadSpec,
    /// Override the nut thickness (None = ISO 4032 nominal).
    pub thickness: Option<f64>,
    /// Override the hex width across flats.
    pub across_flats: Option<f64>,
    /// Radial/flank fit clearance to the nominal external thread
    /// (0 = exact basic profile; 0.15–0.3 typical for FDM printing).
    pub clearance: f64,
    /// 30° hex corner chamfers on both faces (default on).
    pub corner_chamfer: bool,
    /// 45° thread-mouth lead-in chamfers on both faces (default on).
    pub mouth_chamfer: bool,
}

/// Build a hex nut.
///
/// # Errors
///
/// As [`threaded_rod`]; additionally
/// [`FastenerError::ProfileOverlap`] when the clearance consumes the
/// thread form.
pub fn hex_nut(params: &NutParams, tolerance: Tolerance) -> Result<BRep, FastenerError> {
    let stage = "hex nut";
    params.spec.validate()?;
    let m = params
        .thickness
        .unwrap_or_else(|| params.spec.iso_hex_nut_thickness());
    let s = params
        .across_flats
        .unwrap_or_else(|| params.spec.iso_hex_head().across_flats);
    positive(stage, "thickness", m)?;
    positive(stage, "across_flats", s)?;
    let ridge = internal_ridge(&params.spec, params.clearance)?;
    let r_bore = internal_bore_radius(&params.spec, params.clearance);
    let mut body = extrude_shape(hexagon(s), 0.0, m, "nut blank", tolerance)?;
    if params.corner_chamfer {
        body = subtract(
            &body,
            &hex_corner_ring(s, 0.0, -1.0, tolerance)?,
            "nut corner chamfer",
            tolerance,
        )?;
        body = subtract(
            &body,
            &hex_corner_ring(s, m, 1.0, tolerance)?,
            "nut corner chamfer",
            tolerance,
        )?;
    }
    let bore = extrude_shape(
        ProfileShape::Circle { radius: r_bore },
        -0.1 * m,
        1.2 * m,
        "nut bore",
        tolerance,
    )?;
    body = subtract(&body, &bore, "nut bore", tolerance)?;
    body = cut_thread(
        &body,
        &params.spec,
        &ridge,
        thread_form::internal_phase_start(0.0, &params.spec),
        m,
        None,
        None,
        "internal thread",
        tolerance,
    )?;
    if params.mouth_chamfer {
        // Lead-in mouth chamfers on both faces, cut after the thread
        // so the thread ends are chamfered into the cone.
        let r_mouth =
            params.spec.major_diameter() / 2.0 + params.spec.pitch() / 8.0 + params.clearance;
        // Strictly INSIDE the bore void (2x the embed margin):
        // r_bore - embed is EXACTLY the thread groove-base radius,
        // and a cutter cylinder coincident with the groove-base
        // bands poisons boolean classification samples ON the
        // shared surface (measured, PK-464: the propagate seed at
        // (r_in, 0, 0.24) read Inside by parity luck and the
        // chamfer cone's void-side band was kept — 5 self-twinned
        // duplicate edges in the serve). The deeper cut is void-only.
        let r_in = r_bore - 2.0 * thread_form::embed_margin(&params.spec);
        let eps = 0.05 * m;
        body = subtract(
            &body,
            &revolved_tool(
                &[(r_in, -eps), (r_mouth + eps, -eps), (r_in, r_mouth - r_in)],
                "thread mouth chamfer",
                tolerance,
            )?,
            "thread mouth chamfer",
            tolerance,
        )?;
        body = subtract(
            &body,
            &revolved_tool(
                &[
                    (r_in, m + eps),
                    (r_mouth + eps, m + eps),
                    (r_in, m - (r_mouth - r_in)),
                ],
                "thread mouth chamfer",
                tolerance,
            )?,
            "thread mouth chamfer",
            tolerance,
        )?;
    }
    Ok(body)
}
