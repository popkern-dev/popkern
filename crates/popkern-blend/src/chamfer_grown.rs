//! Chamfers whose bevels consume a face, grown from a smaller size
//! (PK-1001).
//!
//! **The class.** Planar chamfers wide enough that the bevels meet
//! beyond a face between them, so the exact result has no piece of it
//! (two adjacent slant-edge bevels of a narrow-faced cone). The
//! one-pass pipeline trims every face by its rails and cannot delete a
//! face, so it declines; the same request at a smaller size serves.
//!
//! **Construction.** When the one-pass pipeline declines and every
//! requested edge is a line between planar faces whose bevel is a
//! plane, the request runs at `1/2`, then `1/4`, then `1/8` of every
//! distance until one serves; that body's bevels are then grown to
//! their requested planes by [`crate::construct::bevel_growth`],
//! which deletes consumed faces and flips vanishing edges at the exact
//! instants they happen. Each bevel face is matched to its spec by
//! its plane at the smaller size (exactly one face of the served body
//! lies on it), and the bevel's offset from its edge must scale with
//! the distance (the growth is a pure translation).
//!
//! **Gates.** The grown body passes the pipeline's serve gate
//! (structural validation and check 9, curve-surface consistency)
//! after the same pcurve maintenance. **Naming:** the stages run
//! untracked and the result is named once against the original input.
//! **Fallback:** an ineligible request, no serving fraction, a failed
//! growth or a failed gate returns the one-pass pipeline's own decline
//! unchanged (the growth's reason under the `bevel-growth` trace
//! topic). **Tolerance:** the caller's; positions gate at `linear × 10`.

use crate::analysis::analyze_edges;
use crate::cancellation::CancellationToken;
use crate::construct::bevel_growth::grow_bevels;
use crate::error::BlendError;
use crate::geometry_compute::compute_chamfer_geometry;
use crate::types::{BlendResult, ChamferKind, ChamferSpec};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

/// The fractions of the requested size tried, largest first.
const GROWTH_FRACTIONS: [f64; 3] = [0.5, 0.25, 0.125];

/// A requested bevel: a point on its edge, its normal, and its plane's
/// offset from that point at the requested size.
struct Bevel {
    edge_point: Point3,
    normal: Vec3,
    offset: f64,
}

fn scaled(kind: &ChamferKind, s: f64) -> ChamferKind {
    match *kind {
        ChamferKind::Equal(d) => ChamferKind::Equal(d * s),
        ChamferKind::TwoDistance {
            face_1,
            dist_1,
            dist_2,
        } => ChamferKind::TwoDistance {
            face_1,
            dist_1: dist_1 * s,
            dist_2: dist_2 * s,
        },
        ChamferKind::DistanceAngle {
            face,
            distance,
            angle,
        } => ChamferKind::DistanceAngle {
            face,
            distance: distance * s,
            angle,
        },
    }
}

/// The bevel planes of `specs` at fraction `s`, when every spec is a
/// planar bevel of a line edge.
fn bevel_planes(
    brep: &BRep,
    specs: &[ChamferSpec],
    s: f64,
    tolerance: &Tolerance,
) -> Option<Vec<(Point3, Vec3)>> {
    let edges: Vec<_> = specs.iter().map(|spec| spec.edge).collect();
    let analyses = analyze_edges(brep, &edges, false).ok()?;
    if analyses.len() != specs.len() {
        return None;
    }
    specs
        .iter()
        .zip(&analyses)
        .map(|(spec, analysis)| {
            let geometry =
                compute_chamfer_geometry(analysis, &scaled(&spec.kind, s), brep, tolerance).ok()?;
            match geometry.surface {
                SurfaceGeometry::Plane { origin, normal, .. } => {
                    Some((origin, normal.normalized()?))
                }
                _ => None,
            }
        })
        .collect()
}

/// The requested bevels, when every spec is a planar bevel of an edge
/// whose ends are vertices.
fn requested_bevels(
    brep: &BRep,
    specs: &[ChamferSpec],
    tolerance: &Tolerance,
) -> Option<Vec<Bevel>> {
    let planes = bevel_planes(brep, specs, 1.0, tolerance)?;
    specs
        .iter()
        .zip(planes)
        .map(|(spec, (origin, normal))| {
            let he = brep.half_edges.get(brep.edges.get(spec.edge)?.half_edge)?;
            let twin = brep.half_edges.get(he.twin)?;
            let (a, b) = (
                brep.vertices.get(he.origin)?.position,
                brep.vertices.get(twin.origin)?.position,
            );
            let edge_point = a + (b - a) * 0.5;
            Some(Bevel {
                edge_point,
                normal,
                offset: (origin - edge_point).dot(&normal),
            })
        })
        .collect()
}

/// Serve `specs` by growth after the one-pass pipeline returned
/// `declined` (see the module docs). `pipeline` is the untracked
/// single pass.
pub(crate) fn chamfer_grown(
    brep: &BRep,
    specs: &[ChamferSpec],
    tolerance: Tolerance,
    tracking: Option<&mut TrackingContext>,
    cancel: Option<&CancellationToken>,
    declined: BlendError,
    pipeline: impl Fn(
        &BRep,
        &[ChamferSpec],
        Option<&CancellationToken>,
    ) -> Result<BlendResult, BlendError>,
) -> Result<BlendResult, BlendError> {
    if !matches!(declined, BlendError::PartialFailure { .. }) {
        return Err(declined);
    }
    let Some(bevels) = requested_bevels(brep, specs, &tolerance) else {
        return Err(declined);
    };
    let trace = popkern_core::trace::trace_enabled("bevel-growth");
    let gate = (tolerance.linear() * 10.0).max(1e-9);
    for s in GROWTH_FRACTIONS {
        let small: Vec<ChamferSpec> = specs
            .iter()
            .map(|spec| ChamferSpec {
                edge: spec.edge,
                kind: scaled(&spec.kind, s),
            })
            .collect();
        let served = match pipeline(brep, &small, cancel) {
            Ok(served) => served,
            Err(BlendError::Cancelled) => return Err(BlendError::Cancelled),
            Err(_) => continue,
        };
        match grow(brep, &small, &bevels, served, s, gate, &tolerance) {
            Ok(mut result) => {
                if trace {
                    eprintln!(
                        "[bevel-growth] served from {s} of the requested size: {} flips, {} \
                         faces consumed",
                        result.repair_report.bevel_growth_flips,
                        result.repair_report.bevel_consumed_faces
                    );
                }
                if let Some(ctx) = tracking {
                    let pre = crate::naming::PreOperationSnapshot::capture(brep);
                    crate::naming::assign_fillet_names(&mut result.brep, ctx, &pre.faces);
                    crate::naming::assign_output_edge_vertex_names(
                        &mut result.brep,
                        ctx,
                        &pre,
                        &tolerance,
                    );
                }
                return Ok(result);
            }
            Err(why) => {
                if trace {
                    eprintln!("[bevel-growth] from {s} of the requested size: {why}");
                }
            }
        }
    }
    Err(declined)
}

/// Grow `served` (the request at fraction `s`) to the requested bevels
/// and run the serve gate.
fn grow(
    brep: &BRep,
    small: &[ChamferSpec],
    bevels: &[Bevel],
    served: BlendResult,
    s: f64,
    gate: f64,
    tolerance: &Tolerance,
) -> Result<BlendResult, String> {
    let starts =
        bevel_planes(brep, small, 1.0, tolerance).ok_or("the smaller bevels are not planar")?;
    let mut body = served.brep;
    let mut growing: Vec<(FaceHandle, (Point3, Vec3))> = Vec::new();
    for (bevel, (origin, normal)) in bevels.iter().zip(starts) {
        // The growth is a translation: same normal, offset scaled by s.
        if normal.cross(&bevel.normal).length() > 1e-9
            || ((origin - bevel.edge_point).dot(&bevel.normal) - s * bevel.offset).abs() > gate
        {
            return Err("a bevel does not scale with its distance".into());
        }
        let on: Vec<FaceHandle> = body
            .faces
            .keys()
            .filter(|&f| {
                let ring = body.face_outer_loop_positions(f);
                !ring.is_empty()
                    && ring
                        .iter()
                        .all(|p| (*p - origin).dot(&normal).abs() <= gate)
            })
            .collect();
        let [face] = on[..] else {
            return Err(format!("{} served faces lie on a bevel's plane", on.len()));
        };
        let end = bevel.edge_point + bevel.normal * bevel.offset;
        growing.push((face, (end, bevel.normal)));
    }
    let report = grow_bevels(&mut body, &growing, gate)?;
    body.fix_stale_edge_windows();
    body.fix_stale_pcurves();
    body.recompute_all_pcurves_lenient();
    let _ = crate::construct::clamp_bounded_pcurves(&mut body);
    let mut repair_report = served.repair_report;
    repair_report.orientation_flag_reconciliations += body.reconcile_orientation_flags();
    body.validate()
        .and_then(|()| body.validate_curve_surface())
        .map_err(|e| format!("the grown body fails the serve gate: {e}"))?;
    body.repair_canonical_forward();
    repair_report.bevel_growth_flips = report.edge_flips;
    repair_report.bevel_consumed_faces = report.consumed_faces.len();
    Ok(BlendResult {
        brep: body,
        edge_warnings: served.edge_warnings,
        corner_patches: served.corner_patches,
        repair_report,
    })
}
