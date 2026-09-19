//! The thread-cut engine: builds helical tool bodies from
//! [`ToolProfile`](crate::thread_form::ToolProfile) cross-sections
//! and subtracts them, one boolean per thread start — plus the shared
//! brep-level construction helpers (validated extrudes, revolved
//! tools, rim chamfers) the builders compose.
//!
//! Everything here works in the canonical fastener frame (axis = +Z);
//! callers place holes by transforming the tools, fasteners by
//! transforming the finished body.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{AnalyticProfile, ProfileShape, SketchPlane};
use popkern_sweep::{
    FrameMethod, HelixPitch, HelixSpec, SweepOptions, SweepProfile, helical_sweep, sweep,
};

use crate::FastenerError;
use crate::spec::{Handedness, ThreadSpec};
use crate::thread_form::ToolProfile;

/// Measured frame mapping for helical sweeps (pinned by
/// `frame_convention_*` in `engine_tests.rs`): the profile +x axis
/// maps to the frame normal, which points INWARD (toward the axis),
/// so an outward radial coordinate authors as `−r`.
///
/// The engine sweeps with an **axis-referenced auxiliary frame**
/// (`FrameMethod::Auxiliary` with the helix axis as the up-curve),
/// not Frenet: the kernel's stitched rational helix is exact in
/// radius but its axial parameterization ripples ≈ 0.5% of the lead
/// within each quarter arc (z(t) linear in the true azimuth is not
/// rationally representable), and the Frenet normal — built from the
/// curve's second derivative — inherits that ripple as a ≈ 0.7°
/// rotation of the thread form about the spine (measured on the P1
/// probe). The axis-referenced normal is exactly radial at every
/// frame sample regardless of the parameterization.
const X_PER_RADIAL: f64 = -1.0;

/// Measured binormal direction along the axis (same pin): +axis for
/// a right-handed helix, −axis for a left-handed one — the profile
/// +y axis maps to the binormal (tilted off the axis by the true
/// lead angle: the thread cross-section lies perpendicular to the
/// helix, the standard modeled-thread geometry).
fn y_per_axial(handedness: Handedness) -> f64 {
    match handedness {
        Handedness::Right => 1.0,
        Handedness::Left => -1.0,
    }
}

/// Where a thread cut lands along +Z: the groove centreline crosses
/// `z_start` on its first nominal turn and `z_end` on its last; the
/// tool helix overshoots each end by half a lead (runout) and the
/// realized turn count rounds UP to the next quarter turn, where the
/// stitched rational helix is exact (`helix.rs` ceils to 90° arcs —
/// fractional remainders would distort the radius).
pub(crate) fn plan_turns(z_start: f64, z_end: f64, lead: f64) -> (f64, f64) {
    let z0 = z_start - 0.5 * lead;
    let span = z_end + 0.5 * lead - z0;
    let turns = (span / lead * 4.0).ceil() / 4.0;
    (z0, turns)
}

/// Validate a constructed body the way a served product must hold:
/// structural + geometric validation plus the serving-boundary
/// curve-surface gate, and never empty.
pub(crate) fn validated(
    brep: BRep,
    stage: &str,
    _tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    if brep.faces.is_empty() {
        return Err(FastenerError::Construction {
            stage: stage.to_string(),
            detail: "produced an empty body (0 faces)".into(),
        });
    }
    brep.validate()
        .map_err(|source| FastenerError::Validation {
            stage: stage.to_string(),
            source,
        })?;
    brep.validate_curve_surface()
        .map_err(|source| FastenerError::Validation {
            stage: stage.to_string(),
            source,
        })?;
    Ok(brep)
}

/// Boolean subtraction with the fastener stage named.
pub(crate) fn subtract(
    body: &BRep,
    tool: &BRep,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    popkern_bool::boolean(
        popkern_bool::BoolOperation::Subtraction,
        body.clone(),
        tool.clone(),
        tolerance,
        None,
    )
    .map_err(|source| FastenerError::Boolean {
        stage: stage.to_string(),
        source,
    })
    .and_then(|r| validated(r, stage, tolerance))
}

/// Boolean union with the fastener stage named.
pub(crate) fn union(
    body: &BRep,
    other: &BRep,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    popkern_bool::boolean(
        popkern_bool::BoolOperation::Union,
        body.clone(),
        other.clone(),
        tolerance,
        None,
    )
    .map_err(|source| FastenerError::Boolean {
        stage: stage.to_string(),
        source,
    })
    .and_then(|r| validated(r, stage, tolerance))
}

/// Subtract one helical tool per thread start from `body`.
///
/// `stage` names the fastener stage for diagnostics (e.g.
/// `"external thread"`). `tool_floor` clamps how far below `z_start`
/// the tool helix may begin: a free-running thread end runs half a
/// lead past it (rod ends, through holes), but a thread stopping
/// against material (a bolt head, a countersunk cone) must keep the
/// tool clear of it. `placement` moves the canonical +Z tools onto
/// an arbitrary hole axis. The subtraction is transactional per
/// start; a failure leaves the input untouched and carries the
/// kernel boolean diagnostics.
///
/// One tool per start covers the WHOLE thread: cutting long threads
/// in abutting runs was tried and declines — a later run's tool
/// surfaces are exact tangent-continuations of the groove walls the
/// earlier run just cut, meeting at the shared cap station (a
/// coincident/tangent boolean class beyond the flush planar caps the
/// kernel has closed). The flank×shank intersections are traced by
/// the radial level-curve arm (PK-461,
/// `popkern_core::geometry::intersection::radial_level_curve`),
/// which has no step budget — the historical ~35-turn bound came
/// from the marching solver's 200k-step budget. The serve frontier
/// is currently the arm's 40-wrap containment (PK-494: the
/// boolean/mesh substrate past 40 wraps): beyond it, and whenever
/// the arm's structure gates decline a pair, marching runs, its
/// budget truncates, and the merge's strand gate declines honestly.
#[allow(clippy::too_many_arguments)]
pub(crate) fn cut_thread(
    body: &BRep,
    spec: &ThreadSpec,
    profile: &ToolProfile,
    z_start: f64,
    z_end: f64,
    tool_floor: Option<f64>,
    placement: Option<&Transform3>,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let (mut z0, turns) = plan_turns(z_start, z_end, spec.lead());
    if let Some(floor) = tool_floor {
        z0 = z0.max(floor);
    }
    // Serve containment (PK-494): the boolean/mesh substrate is
    // verified through 46-tool-turn thread serves (per-turn probe
    // ladder; PK-494 layers 1-2 moved the bound from 40). Beyond it
    // the chart-corner triangle's inverted winding cascades into
    // either an honest cellular decline or a graph-merge loop-gap —
    // and historically past-frontier serves have TORN instead of
    // declining, which is worse. Declining HERE keeps the failure
    // structured and instant; the radial level-curve arm's own wrap
    // containment is defense-in-depth for other callers. PK-494
    // layer 3 (region-adjacency seam-side assignment in the walk)
    // lifts this further.
    const MAX_SERVED_TURNS: f64 = 46.5;
    // POPKERN_MULTIWRAP_SERVE_CONTAINMENT=0 lifts the gate for
    // PK-494's frontier probes (PK-239 kill-switch precedent;
    // per-call read) — never for production serves.
    let containment_lifted = std::env::var_os("POPKERN_MULTIWRAP_SERVE_CONTAINMENT")
        .map(|v| v == "0")
        .unwrap_or(false);
    if turns > MAX_SERVED_TURNS && !containment_lifted {
        return Err(FastenerError::Construction {
            stage: stage.to_string(),
            detail: format!(
                "thread cut spans {turns} tool turns, past the verified {MAX_SERVED_TURNS}-turn \
                 serve frontier (PK-494: multi-wrap boolean/mesh substrate) — split the thread \
                 into shorter segments or reduce the threaded length"
            ),
        });
    }
    // ONE boolean for ALL starts: the starts' grooves are cut from a
    // single multi-shell tool (the disjoint groove coils
    // concatenated). Sequential per-start subtracts were measured
    // silently dropping the SECOND start's groove: the first cut
    // turns the shank wall into multi-wrap helical band faces whose
    // chart embedding spans several periods, the chart substrate
    // declines ("anchored interior spans 19.3 > period 6.28"), the
    // 3D-Newell fallback cannot partition a non-contractible band,
    // and the second subtract served the body unchanged — an 8.3%
    // volume lie on an M8x1 two-start rod (PK-460). Against the
    // INTACT wall, one cut with every groove present partitions with
    // the tuned machinery.
    let mut combined: Option<BRep> = None;
    for start_index in 0..spec.starts() {
        let mut tool = helical_tool(spec, profile, z0, turns, start_index, stage, tolerance)?;
        if let Some(transform) = placement {
            tool = tool
                .transform(transform)
                .map_err(|e| FastenerError::Construction {
                    stage: stage.to_string(),
                    detail: format!("tool placement transform failed: {e}"),
                })?;
        }
        combined = Some(match combined {
            None => tool,
            Some(existing) => popkern_bool::concatenate_solids(existing, tool).map_err(|e| {
                FastenerError::Construction {
                    stage: stage.to_string(),
                    detail: format!("thread tool concatenation failed: {e}"),
                }
            })?,
        });
    }
    let Some(tool) = combined else {
        return Ok(body.clone());
    };
    // The PK-475 tool-sense containment (reverse the tool to match a
    // folded-inward body) is retired: the boolean normalizes operand
    // material sense at entry (PK-478 stage 4), so mixed-sense pairs
    // are handled by the kernel, not per caller.
    subtract(body, &tool, stage, tolerance)
}

/// Build the helical tool body for one thread start.
pub(crate) fn helical_tool(
    spec: &ThreadSpec,
    profile: &ToolProfile,
    z0: f64,
    turns: f64,
    start_index: u32,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let axial = y_per_axial(spec.handedness());
    let points: Vec<[f64; 2]> = profile
        .points
        .iter()
        .map(|&(r, y)| [X_PER_RADIAL * r, axial * y])
        .collect();
    let wire = polygon_wire(&points, stage)?;
    // The sweep recentres the wire on its own sampled centroid and
    // places that centroid ON the helix; replicate the sampling so
    // the helix radius and axial origin land the authored profile at
    // its absolute coordinates.
    let (cx, cy) = sampled_centroid(&wire);
    let radius = cx * X_PER_RADIAL;
    let phase = 2.0 * std::f64::consts::PI * f64::from(start_index) / f64::from(spec.starts());
    let z_base = z0 + cy * axial;
    let helix = HelixSpec {
        axis_origin: Point3::new(0.0, 0.0, z_base),
        axis_dir: Vec3::new(0.0, 0.0, 1.0),
        start_radius: radius,
        end_radius: radius,
        pitch: HelixPitch::Constant(spec.lead()),
        turns,
        start_angle: phase,
        left_handed: spec.handedness() == Handedness::Left,
    };
    let options = SweepOptions {
        frame: FrameMethod::Auxiliary(axis_line(z_base, z_base + turns * spec.lead())),
        ..SweepOptions::default()
    };
    let result = helical_sweep(&SweepProfile::Wire(wire), &helix, &options, tolerance, None)
        .map_err(|source| FastenerError::Sweep {
            stage: format!("{stage} tool (start {start_index}, {turns} turns)"),
            source,
        })?;
    validated(
        result.brep,
        &format!("{stage} tool (start {start_index}, {turns} turns)"),
        tolerance,
    )
}

/// The helix axis as an auxiliary up-curve: a straight segment from
/// `z_from` to `z_to` on the z-axis, matched to the tool helix so
/// the proportional evaluation stays alongside the path.
pub(crate) fn axis_line(z_from: f64, z_to: f64) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(0.0, 0.0, z_from), Point3::new(0.0, 0.0, z_to)],
        weights: vec![1.0, 1.0],
    }
}

/// Build an axisymmetric tool body by revolving a closed polygon in
/// the meridian half-plane — points are `(radial, z)` with
/// `radial ≥ 0` — a full turn about +Z. Cone mouths, countersinks,
/// hex corner chamfers, and countersunk heads are all built with
/// these.
pub(crate) fn revolved_tool(
    points: &[(f64, f64)],
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let wire = AnalyticProfile {
        shape: ProfileShape::Polygon {
            points: points.iter().map(|&(r, z)| [r, z]).collect(),
        },
        center: [0.0, 0.0],
    }
    .placed_wire(&SketchPlane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, -1.0, 0.0),
    })
    .map_err(|e| FastenerError::Construction {
        stage: stage.to_string(),
        detail: format!("meridian profile compile failed: {e}"),
    })?;
    let result = popkern_sweep::revolve(
        &SweepProfile::Wire(wire),
        &Point3::new(0.0, 0.0, 0.0),
        &Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::TAU,
        tolerance,
        None,
    )
    .map_err(|source| FastenerError::Sweep {
        stage: format!("{stage} (revolved tool)"),
        source,
    })?;
    validated(result.brep, stage, tolerance)
}

/// Extrude a closed profile shape from `z0` by `depth` along +Z (the
/// same exact compile-and-sweep path `popkern-part`'s `Profile` uses,
/// including the sampled-centroid compensation that keeps authored
/// coordinates absolute).
pub(crate) fn extrude_shape(
    shape: ProfileShape<f64>,
    z0: f64,
    depth: f64,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    let wire = AnalyticProfile {
        shape,
        center: [0.0, 0.0],
    }
    .local_wire(&SketchPlane {
        origin: Point3::new(0.0, 0.0, z0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    })
    .map_err(|e| FastenerError::Construction {
        stage: stage.to_string(),
        detail: format!("profile compile failed: {e}"),
    })?;
    // The sweep recentres the wire on its sampled centroid; fold the
    // same centroid into the path so authored coordinates hold.
    let (cx, cy) = sampled_centroid(&wire);
    let path = NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![Point3::new(cx, cy, z0), Point3::new(cx, cy, z0 + depth)],
        weights: vec![1.0, 1.0],
    };
    let result = sweep(
        &SweepProfile::Wire(wire),
        &path,
        &SweepOptions::default(),
        tolerance,
        None,
    )
    .map_err(|source| FastenerError::Sweep {
        stage: stage.to_string(),
        source,
    })?;
    validated(result.brep, stage, tolerance)
}

/// Equal-distance chamfer on every circular edge at height `z` (the
/// shank tip rims and nut bore rims — closed circles whose centres
/// sit on the fastener axis).
pub(crate) fn chamfer_circular_rims_at(
    body: &BRep,
    z: f64,
    distance: f64,
    stage: &str,
    tolerance: Tolerance,
) -> Result<BRep, FastenerError> {
    use popkern_blend::{ChamferKind, ChamferOptions, ChamferSpec, chamfer_edges};
    use popkern_core::geometry::curve::CurveGeometry;
    let band = tolerance.linear() * 10.0;
    let specs: Vec<ChamferSpec> = body
        .edges
        .iter()
        .filter_map(|(eh, edge)| {
            let curve = edge.curve.as_ref()?;
            let center_z = match &curve.geometry {
                CurveGeometry::Circle { center, .. } | CurveGeometry::Arc { center, .. } => {
                    center.z
                }
                _ => return None,
            };
            ((center_z - z).abs() < band).then_some(ChamferSpec {
                edge: eh,
                kind: ChamferKind::Equal(distance),
            })
        })
        .collect();
    if specs.is_empty() {
        return Err(FastenerError::Construction {
            stage: stage.to_string(),
            detail: format!("no circular rim found at z = {z} to chamfer"),
        });
    }
    let result = chamfer_edges(body, &specs, &ChamferOptions::default(), tolerance, None).map_err(
        |source| FastenerError::Blend {
            stage: stage.to_string(),
            source,
        },
    )?;
    validated(result.brep, stage, tolerance)
}

/// Compile a closed polygon to an exact NURBS wire in the sweep
/// profile's local frame (the same compiler every sketch uses).
pub(crate) fn polygon_wire(points: &[[f64; 2]], stage: &str) -> Result<NurbsCurve, FastenerError> {
    AnalyticProfile {
        shape: ProfileShape::Polygon {
            points: points.to_vec(),
        },
        center: [0.0, 0.0],
    }
    .local_wire(&SketchPlane {
        origin: Point3::new(0.0, 0.0, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    })
    .map_err(|e| FastenerError::Construction {
        stage: stage.to_string(),
        detail: format!("profile compile failed: {e}"),
    })
}

/// The sweep kernel's wire recentring: the centroid of 64 uniform
/// parameter samples (mirrors `Profile::extrude`; the absolute
/// thread-surface pins in the engine tests catch any drift between
/// the two).
pub(crate) fn sampled_centroid(wire: &NurbsCurve) -> (f64, f64) {
    use popkern_core::geometry::curve::CurveEval;
    let (t0, t1) = wire.parameter_range();
    let n = 64;
    let mut sx = 0.0;
    let mut sy = 0.0;
    for i in 0..n {
        let t = t0 + (t1 - t0) * (i as f64) / (n as f64);
        let p = wire.point_at(t);
        sx += p.x;
        sy += p.y;
    }
    (sx / n as f64, sy / n as f64)
}
