//! Window-SET clipping of intersection records to face boundaries.
//!
//! A record curve can enter and leave a face's trimmed domain several
//! times — an analytic sphere×plane circle crosses a cube face's
//! bounds four times, a line through a reflex outline twice per bite.
//! The per-face verdict is therefore a set of in-face parameter
//! windows ([`UvClipResult::Windows`]), and the joint verdict is the
//! pairwise intersection of the two faces' window sets (mod the
//! parameter period for closed records). Every surviving window
//! becomes its own clipped `CurveOnSurface` — collapsing them to one
//! span (the pre-PK-111 behavior) manufactured records that ran
//! THROUGH regions outside the face and over foreign junctions.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::curve_on_surface::CurveOnSurface;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::topology::{BRep, FaceHandle};

use super::clip_uv::{UvClipResult, clip_curve_to_crossings, clip_curve_wrapped, uv_clip_one_face};

#[path = "clip_window_sets.rs"]
mod clip_window_sets;
use clip_window_sets::{
    intersect_window_sets, record_is_transversal, subtract_window_sets, trace_record_geometry,
};

/// Trim-residue emission gate (PK-465): residue records partition a
/// face wherever the OTHER operand's surface crosses it beyond the
/// shared trim window. Built end to end (emission here, kept-flank
/// filtering at insertion, single-face registration, slit-aware wedge
/// matching in the walk).
const EMIT_TRIM_RESIDUES: bool = true;

/// Clip an intersection record to both faces' trimmed UV domains.
/// Returns EVERY in-face span as its own `CurveOnSurface` (empty =
/// the record misses at least one face entirely). A record judged
/// in-face over its whole closed range survives as a single closed
/// clone.
pub(crate) fn clip_curve_uv_space(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
    residues: Option<&mut Vec<crate::clip::TrimResidue>>,
    cache: &mut super::clip_pierce_windows::PierceInventoryCache,
) -> Vec<CurveOnSurface> {
    // OFF-PRINCIPAL-BRANCH RETRY (PK-480, PK-564): a record pcurve can
    // sit whole periods from a face whose chart occupies a NEGATIVE
    // branch (the angled tap's entry-cone face at u ∈ [−2π, 0]; the
    // cone×top-face conic was authored at u ∈ [2π, 4π]) or OVERFLOWS
    // the principal cell (a chamfer cone charted on [π/2, 5π/2]). The tuned single-period
    // normalization inside crossing detection FOLDS into [0, period)
    // and provably cannot reach such a chart — and it is internal, so
    // every containment vote, window verdict, and residue judged the
    // raw branch anyway: the closed conic dropped whole and the
    // served body kept a phantom flap above the entry face (21 open
    // mesh edges). The retry is scoped to NOTHING-PRODUCED: the tuned
    // pass runs first, byte-identical, and only a record it clipped
    // away entirely (no pieces, no residues) is re-clipped with its
    // pcurve translated onto the face's branch — a record that
    // produced anything keeps its tuned output (the flush-gusset
    // union has a negative-branch face whose producing records must
    // not move).
    let mut first_residues: Vec<crate::clip::TrimResidue> = Vec::new();
    let out = clip_curve_uv_space_tuned(
        brep_a,
        face_a,
        brep_b,
        face_b,
        curve,
        Some(&mut first_residues),
        cache,
    );
    let adopt = |chosen: Vec<crate::clip::TrimResidue>| {
        if let Some(res) = residues {
            res.extend(chosen);
        }
    };
    if out.is_empty() && first_residues.is_empty() {
        let a = align_off_principal_branch_side(brep_a, face_a, &curve.curve_uv1);
        let b = align_off_principal_branch_side(brep_b, face_b, &curve.curve_uv2);
        if a.is_some() || b.is_some() {
            let mut adjusted = curve.clone();
            if let Some(c) = a {
                adjusted.curve_uv1 = c;
            }
            if let Some(c) = b {
                adjusted.curve_uv2 = c;
            }
            if popkern_core::trace::trace_enabled("clip") {
                eprintln!(
                    "  [uvclip] {face_a:?}x{face_b:?}: off-principal-branch retry (tuned pass empty)"
                );
            }
            let mut retry_residues: Vec<crate::clip::TrimResidue> = Vec::new();
            let retry = clip_curve_uv_space_tuned(
                brep_a,
                face_a,
                brep_b,
                face_b,
                &adjusted,
                Some(&mut retry_residues),
                cache,
            );
            if !retry.is_empty() || !retry_residues.is_empty() {
                adopt(retry_residues);
                return retry;
            }
        }
    }
    adopt(first_residues);
    out
}

/// The tuned clip body — see [`clip_curve_uv_space`] for the
/// off-principal-branch retry wrapper.
fn clip_curve_uv_space_tuned(
    brep_a: &BRep,
    face_a: FaceHandle,
    brep_b: &BRep,
    face_b: FaceHandle,
    curve: &CurveOnSurface,
    residues: Option<&mut Vec<crate::clip::TrimResidue>>,
    cache: &mut super::clip_pierce_windows::PierceInventoryCache,
) -> Vec<CurveOnSurface> {
    let trace = popkern_core::trace::trace_enabled("clip");
    if trace {
        trace_record_geometry(brep_b, face_a, face_b, curve);
    }
    // Open records may take decisive multi-windows only when the
    // record is TRANSVERSAL — the two faces' surface normals disagree
    // along it. Flush/coincident boundary-running traces (parallel
    // normals) keep the tuned crossing-hull contract (PK-114).
    // SCOPE (PK-472): pierce authority and residue emission serve
    // the graze-junction classes of ANALYTIC/swept faces (a handful
    // of boundary edges — tap, nut, rod mouths). A sketch-derived
    // polygon face defeats both: the pierce inventory is
    // O(boundary-edges × records) (the impeller's 406-point star ×
    // its blade records turned a 40 ms subtraction into 935 s), and
    // its residue swarm destabilized the star-cut walk (broken
    // half-edge cycle, 1.4e-2 loop gap — a served baseline). Faces
    // beyond the edge-count scope keep the tuned machinery
    // byte-identical.
    let in_pierce_scope = |brep: &BRep, face: FaceHandle| -> bool {
        crate::clip_endpoint_refine::boundary_edges(brep, face).len() <= 64
    };
    let pierce_scope = in_pierce_scope(brep_a, face_a) && in_pierce_scope(brep_b, face_b);
    let open_windows_allowed =
        !curve.is_closed && record_is_transversal(brep_a, face_a, brep_b, face_b, curve);
    // Chart-cell route (PK-304): doubly-periodic faces take
    // folded-containment windows (see `chart_cell`); everything else
    // keeps the tuned crossing machinery.
    let cell_clip = super::chart_cell::cell_clip;
    let clip_a = cell_clip(brep_a, face_a, &curve.curve_uv1).unwrap_or_else(|| {
        uv_clip_one_face(
            brep_a,
            face_a,
            &curve.curve_uv1,
            curve.is_closed,
            open_windows_allowed,
        )
    });
    let clip_b = cell_clip(brep_b, face_b, &curve.curve_uv2).unwrap_or_else(|| {
        uv_clip_one_face(
            brep_b,
            face_b,
            &curve.curve_uv2,
            curve.is_closed,
            open_windows_allowed,
        )
    });

    let name = |r: &UvClipResult| match r {
        UvClipResult::NoPcurves => "NoPcurves".to_string(),
        UvClipResult::Outside => "Outside".to_string(),
        UvClipResult::Windows(ws) => format!(
            "Windows[{}]",
            ws.iter()
                .map(|(e, x)| format!("({e:.3},{x:.3})"))
                .collect::<Vec<_>>()
                .join(",")
        ),
    };
    if matches!(clip_a, UvClipResult::Outside) || matches!(clip_b, UvClipResult::Outside) {
        if trace {
            eprintln!(
                "  [uvclip] dropped: A={} B={}",
                name(&clip_a),
                name(&clip_b)
            );
        }
        // DROPPED-SIDE RESIDUES (PK-465): a record fully outside ONE
        // face's trim still partitions the OTHER face wherever that
        // side produced real windows — the record's surface keeps
        // crossing it. Measured on the angled tap: the thread tool's
        // ridge-root band runs flush along the drilled bore wall and
        // exits the plate through the mouth opening; its plate-top
        // crossing records clip to A=Outside (they lie inside the
        // mouth ellipse), the band kept walking top-to-bottom as one
        // region, classified whole from a below-plate sample, and
        // the served body carried the band's above-plate flap as a
        // masked ghost weave (55 self-twinned strands). Emit the
        // windowed side's spans as single-face residues instead of
        // discarding the partition outright. Same gating as the
        // joint-window residue block below.
        // Scope: only records whose Outside verdict is the IN-HOLE
        // shape — inside the face's OUTER trim but inside an inner
        // loop (the other operand's surface exits through the
        // opening; the mouth of a drilled hole is the measured
        // case). A record truly outside the outer trim partitions
        // nothing real on the other face's material — emitting those
        // as slits wove phantom boundaries through a flush union
        // (measured: the project-box boss union's serve deviated
        // 7.8e-3 with unscoped emission).
        let outside_is_in_hole = |outside_a: bool| -> bool {
            use popkern_core::geometry::curve::CurveEval as _;
            let (brep, face, uv) = if outside_a {
                (brep_a, face_a, &curve.curve_uv1)
            } else {
                (brep_b, face_b, &curve.curve_uv2)
            };
            let (t0, t1) = uv.parameter_range();
            let q = uv.point_at(0.5 * (t0 + t1));
            crate::split_inner_loop::point_inside_face_uv(brep, face, q, &brep.tolerance)
                && !crate::split_inner_loop::point_inside_face_region_uv(
                    brep,
                    face,
                    q,
                    &brep.tolerance,
                )
        };
        if let Some(res) = residues
            && EMIT_TRIM_RESIDUES
            && open_windows_allowed
            && pierce_scope
            && outside_is_in_hole(matches!(clip_a, UvClipResult::Outside))
        {
            let (t_min, t_max) = curve.curve_3d.parameter_range();
            let span_eps = (t_max - t_min).abs().max(1e-12) * 1e-6;
            let mut push_windows = |ws: &[(f64, f64)], owner_is_a: bool| {
                for &(t_enter, t_exit) in ws {
                    if t_exit - t_enter < span_eps || t_exit > t_max + 1e-12 {
                        continue;
                    }
                    let (Some(clipped_3d), Some(clipped_uv1), Some(clipped_uv2)) = (
                        clip_curve_to_crossings(&curve.curve_3d, t_enter, t_exit),
                        clip_curve_to_crossings(&curve.curve_uv1, t_enter, t_exit),
                        clip_curve_to_crossings(&curve.curve_uv2, t_enter, t_exit),
                    ) else {
                        continue;
                    };
                    res.push((
                        owner_is_a,
                        CurveOnSurface {
                            curve_3d: clipped_3d,
                            curve_uv1: clipped_uv1,
                            curve_uv2: clipped_uv2,
                            is_closed: false,
                        },
                        // No kept sibling on the dropped-side path —
                        // no shared cut bound to keep consistent.
                        (false, false),
                    ));
                }
            };
            match (&clip_a, &clip_b) {
                (UvClipResult::Windows(wa), UvClipResult::Outside) => push_windows(wa, true),
                (UvClipResult::Outside, UvClipResult::Windows(wb)) => push_windows(wb, false),
                _ => {}
            }
        }
        return Vec::new();
    }
    if matches!(clip_a, UvClipResult::NoPcurves) && matches!(clip_b, UvClipResult::NoPcurves) {
        return vec![curve.clone()];
    }
    if trace {
        eprintln!(
            "  [uvclip] {face_a:?}x{face_b:?}: A={} B={}",
            name(&clip_a),
            name(&clip_b)
        );
    }

    let (t_min, t_max) = curve.curve_3d.parameter_range();
    // Residues are only meaningful for a side that produced REAL
    // windows: the NoPcurves full-range default is no trim knowledge,
    // and subtracting from it would fabricate chart-wide residues.
    let a_real = matches!(clip_a, UvClipResult::Windows(_));
    let b_real = matches!(clip_b, UvClipResult::Windows(_));
    let windows_of = |r: UvClipResult| -> Vec<(f64, f64)> {
        match r {
            UvClipResult::Windows(ws) => ws,
            // NoPcurves on one side: that face imposes no constraint.
            _ => vec![(t_min, t_max)],
        }
    };
    let wa = windows_of(clip_a);
    let wb = windows_of(clip_b);
    let mut joint = intersect_window_sets(&wa, &wb, t_min, t_max, curve.is_closed);

    // PIERCE AUTHORITY — mechanism and measured rationale in
    // `clip_pierce_windows::apply_pierce_authority`. Scoped to
    // transversal OPEN records with no wrapped tuned window: the
    // wrapped/closed bookkeeping keeps its tuned path.
    let mut junction_points: Vec<(f64, popkern_core::geometry::Point3)> = Vec::new();
    if open_windows_allowed && pierce_scope && joint.iter().all(|w| w.1 <= t_max + 1e-12) {
        let (j2, jp) = super::clip_pierce_windows::apply_pierce_authority(
            brep_a, face_a, brep_b, face_b, curve, t_min, t_max, joint, trace, cache,
        );
        joint = j2;
        junction_points = jp;
    }
    if trace && joint.len() > 1 {
        eprintln!(
            "  [uvclip] {} in-face windows survive the joint clip",
            joint.len()
        );
    }

    // Pin a piece end whose window bound is a pierce junction onto
    // the junction's exact 3D point (PK-472): the clipped end
    // otherwise evaluates on the record's fitted curve, and two
    // records meeting at one junction mint duplicate nodes a lateral
    // slop apart (~6e-4 measured on the m6 nut's mouth corner — the
    // assembled shell failed validation on the gap).
    let span_all = (t_max - t_min).abs().max(1e-12);
    let surfs = (
        brep_a
            .faces
            .get(face_a)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry),
        brep_b
            .faces
            .get(face_b)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry),
    );
    let pin_gate = brep_a.tolerance.linear().max(brep_b.tolerance.linear()) * 2000.0;
    let pin_junction_ends = |piece: &mut CurveOnSurface, t_enter: f64, t_exit: f64| {
        let (Some(sa), Some(sb)) = surfs else {
            return;
        };
        for &(tj, pj) in &junction_points {
            for (bound, end) in [(t_enter, false), (t_exit, true)] {
                // Window bounds are the pierce params themselves
                // when the override fired — exact match, with a 3D
                // fit-slop gate as the safety.
                if (tj - bound).abs() < span_all * 1.0e-6 {
                    let p_end = if end {
                        let (_, t1) = piece.curve_3d.parameter_range();
                        piece.curve_3d.point_at(t1)
                    } else {
                        let (t0, _) = piece.curve_3d.parameter_range();
                        piece.curve_3d.point_at(t0)
                    };
                    if p_end.distance_to(&pj) < pin_gate {
                        crate::clip_endpoint_refine::set_endpoint_3d_uv(
                            piece,
                            end,
                            pj,
                            sa,
                            sb,
                            &brep_a.tolerance,
                        );
                    }
                }
            }
        }
    };

    // TRIM-RESIDUE extraction (PK-465): the spans inside ONE face's
    // trim but outside the other's still partition that one face —
    // beyond the shared window the other operand's SURFACE keeps
    // crossing this face, and regions on its two sides must not walk
    // as one fragment (the angled tap's flank merged its above-plate
    // and in-material parts around each z-plane record endpoint and
    // classified them with one verdict). Scoped to transversal OPEN
    // records, the same contract as decisive windows.
    // STAGED (PK-465): residue extraction is built end to end
    // (emission here, op-aware kept-flank filtering at insertion,
    // single-face registration, walk slit handling) and measurably
    // completes the angled tap's flank partition — but discard-flank
    // slits still destabilize one tuned flush fixture's face walk
    // (half_wiring self-twins), so emission stays OFF until the slit
    // walk contract is finished. Flip this gate to re-enable.
    if let Some(res) = residues
        && EMIT_TRIM_RESIDUES
        && open_windows_allowed
        && pierce_scope
        && joint.iter().all(|w| w.1 <= t_max + 1e-12)
    {
        let span_eps = (t_max - t_min).abs().max(1e-12) * 1e-6;
        let mut push_residues = |mine: &[(f64, f64)], owner_is_a: bool| {
            for (t_enter, t_exit) in subtract_window_sets(mine, &joint) {
                if t_exit - t_enter < span_eps {
                    continue;
                }
                let (Some(clipped_3d), Some(clipped_uv1), Some(clipped_uv2)) = (
                    clip_curve_to_crossings(&curve.curve_3d, t_enter, t_exit),
                    clip_curve_to_crossings(&curve.curve_uv1, t_enter, t_exit),
                    clip_curve_to_crossings(&curve.curve_uv2, t_enter, t_exit),
                ) else {
                    continue;
                };
                let mut piece = CurveOnSurface {
                    curve_3d: clipped_3d,
                    curve_uv1: clipped_uv1,
                    curve_uv2: clipped_uv2,
                    is_closed: false,
                };
                pin_junction_ends(&mut piece, t_enter, t_exit);
                // Which ends abut a joint window — the shared cut
                // bounds a kept sibling also carries. Those ends
                // (and ONLY those) must take the same endpoint
                // refinement the sibling gets (PK-472).
                let bound_eps = (t_max - t_min).abs().max(1e-12) * 1e-9;
                let shared = (
                    joint.iter().any(|w| {
                        (w.1 - t_enter).abs() < bound_eps || (w.0 - t_enter).abs() < bound_eps
                    }),
                    joint.iter().any(|w| {
                        (w.0 - t_exit).abs() < bound_eps || (w.1 - t_exit).abs() < bound_eps
                    }),
                );
                res.push((owner_is_a, piece, shared));
            }
        };
        if a_real {
            push_residues(&wa, true);
        }
        if b_real {
            push_residues(&wb, false);
        }
    }

    let mut out = Vec::new();
    for (t_enter, t_exit) in joint {
        if t_exit - t_enter < 1e-15 {
            continue;
        }
        if curve.is_closed && (t_enter - t_min).abs() < 1e-12 && (t_exit - t_max).abs() < 1e-12 {
            out.push(curve.clone());
            continue;
        }
        let wrapped = t_exit > t_max + 1e-12;
        let extract = |c: &NurbsCurve| -> Option<NurbsCurve> {
            if wrapped {
                clip_curve_wrapped(c, t_enter, t_exit)
            } else {
                clip_curve_to_crossings(c, t_enter, t_exit)
            }
        };
        let (Some(clipped_3d), Some(clipped_uv1), Some(clipped_uv2)) = (
            extract(&curve.curve_3d),
            extract(&curve.curve_uv1),
            extract(&curve.curve_uv2),
        ) else {
            continue;
        };
        let mut piece = CurveOnSurface {
            curve_3d: clipped_3d,
            curve_uv1: clipped_uv1,
            curve_uv2: clipped_uv2,
            is_closed: false,
        };
        if !wrapped {
            pin_junction_ends(&mut piece, t_enter, t_exit);
        }
        out.push(piece);
    }
    out
}

use super::clip_negative_branch::align_off_principal_branch_side;
