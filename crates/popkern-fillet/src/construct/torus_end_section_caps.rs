//! Torus-strip mid-chain terminations: section-plane end caps
//! (PK-238).
//!
//! **The class.** A fillet strip promoted to an analytic Torus whose
//! spine ends mid-smooth-chain — the filleted edge tangent-continues
//! into un-blended sharp edges (§12.3 blend-on-blend: the strip-end
//! runout arc joins the cube's straight top edges tangentially).
//! The carve leaves each end as an off-surface mitre complex: the
//! support-plane trim is a straight CHORD between overshoot vertices
//! (the edge endpoints pushed one fillet radius along the tangent),
//! joined by stub chords through the surviving sharp corner — every
//! edge of the complex sits ~`r(√2−1)` off the exact torus (measured
//! 8.3e-3..1.1e-2 at r = 0.02, the pinned §12.3 residual class).
//!
//! **Repair, by construction.** A mid-chain fillet terminus is the
//! torus' own axial section at the end station (the plane through
//! the torus axis and the spine end): the strip ends on the exact
//! minor-circle arc, and the material between that arc and the
//! sharp corner is a flat SECTION CAP in that plane — the standard
//! stopped-fillet end wall. Everything is derived from the torus
//! alone:
//!
//! 1. contract the overshoot stubs (`kve` — the sharp corner
//!    vertex survives, restoring the un-blended edges to their true
//!    corner),
//! 2. split the support-plane chord at each end's plane-contact
//!    point `V_t = spine_end + axis·(±r)` and land the split
//!    vertices exactly there,
//! 3. re-author the middle piece as the exact plane-contact circle
//!    arc (radius = major, center on the axis) and the far trim as
//!    the exact support-contact circle through its own endpoints,
//! 4. one `mef` per end from `V_t` to the on-torus trim end `V_b`
//!    across the strip: the new edge is the exact axial-section arc
//!    (minor circle, radius = tube), the new face is the planar
//!    section cap keeping the sharp corner.
//!
//! **Invariants.** Per strip: V +2, E +2 net from the splits/kve
//! pairs and E +2, F +2, L +2 from the two mefs — χ preserved. The
//! strip keeps its handle and Torus; the caps are new faces with
//! the section Plane; corner vertices and the un-blended edges keep
//! their handles. **Fallback:** every strip is repaired on a clone
//! and accepted only when structural validation holds and every
//! touched edge lies on BOTH adjacent surfaces within the
//! acceptance band; a decline keeps the carve topology and reports
//! via the `corner-carve` trace topic. **Tolerances:** on/off-torus
//! classification 1e-6 (class separation is 8e-3 vs 1e-12);
//! authored-curve endpoint and lie-on checks at 1e-9-per-unit
//! scale, both local to this pass.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, HalfEdgeHandle, VertexHandle};

#[path = "torus_end_section_caps_author.rs"]
mod torus_end_section_caps_author;
use torus_end_section_caps_author::*;

/// On/off-torus vertex classification band. The mitre complex sits
/// `r(√2−1)` (8e-3 at r = 0.02) off the torus; true contacts are
/// exact to the contact-solver residual.
const ON_TORUS_BAND: f64 = 1.0e-6;

/// Endpoint hit and lie-on acceptance for authored exact curves.
const EXACT_BAND: f64 = 1.0e-9;

/// Author section caps on every analytic-torus strip whose end
/// boundary carries the off-surface mitre-chord signature. Returns
/// the number of capped ends.
pub(super) fn author_torus_end_section_caps(result: &mut BRep, tolerance: &Tolerance) -> usize {
    let trace = popkern_core::trace::trace_enabled("corner-carve");
    let strips: Vec<(FaceHandle, Point3, Vec3, f64, f64)> = result
        .faces
        .iter()
        .filter_map(|(h, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(&SurfaceGeometry::Torus {
                center,
                axis,
                major_radius,
                minor_radius,
            }) => axis
                .normalized()
                .map(|a| (h, center, a, major_radius, minor_radius)),
            _ => None,
        })
        .collect();

    let mut capped = 0_usize;
    for (strip, c0, a, major, minor) in strips {
        let mut candidate = result.clone();
        match cap_strip_ends(&mut candidate, strip, c0, a, major, minor, tolerance) {
            Ok(0) => {}
            Ok(n) => {
                *result = candidate;
                capped += n;
                if trace {
                    eprintln!("  [torus-end-cap] strip {strip:?}: {n} section caps authored");
                }
            }
            Err(reason) => {
                if trace {
                    eprintln!("  [torus-end-cap] strip {strip:?} declined: {reason}");
                }
            }
        }
    }
    capped
}

/// Distance from `p` to the torus surface.
fn torus_distance(p: &Point3, c0: &Point3, a: &Vec3, major: f64, minor: f64) -> f64 {
    let w = *p - *c0;
    let axial = w.dot(a);
    let radial = (w - *a * axial).length();
    (((radial - major).powi(2) + axial * axial).sqrt() - minor).abs()
}

/// One strip's end reconstruction. Returns the number of capped
/// ends (0 = the strip does not carry the signature — not an
/// error), or a decline reason.
#[allow(clippy::too_many_arguments)]
fn cap_strip_ends(
    brep: &mut BRep,
    strip: FaceHandle,
    c0: Point3,
    a: Vec3,
    major: f64,
    minor: f64,
    tolerance: &Tolerance,
) -> Result<usize, String> {
    // ── Boundary census ─────────────────────────────────────────
    let loop_h = brep
        .faces
        .get(strip)
        .map(|f| f.outer_loop)
        .ok_or("strip face vanished")?;
    let start = brep
        .loops
        .get(loop_h)
        .map(|l| l.half_edge)
        .ok_or("strip loop vanished")?;
    let mut entries: Vec<(HalfEdgeHandle, EdgeHandle, VertexHandle, Point3, bool)> = Vec::new();
    let mut cur = start;
    for _ in 0..brep.half_edges.len() + 1 {
        let he = brep.half_edges.get(cur).ok_or("loop walk broke")?;
        let pos = brep
            .vertices
            .get(he.origin)
            .map(|v| v.position)
            .ok_or("loop vertex vanished")?;
        let on = torus_distance(&pos, &c0, &a, major, minor) <= ON_TORUS_BAND;
        entries.push((cur, he.edge, he.origin, pos, on));
        cur = he.next;
        if cur == start {
            break;
        }
    }
    let n = entries.len();
    let off_count = entries.iter().filter(|e| !e.4).count();
    if off_count == 0 {
        return Ok(0); // exact boundary — nothing to do
    }
    // Exactly one contiguous cyclic off-run of 4 (two corners + two
    // overshoots) is the mitre-chord signature; anything else is a
    // different end class and declines with the census.
    let mut run_start = None;
    for i in 0..n {
        if !entries[i].4 && entries[(i + n - 1) % n].4 {
            if run_start.is_some() {
                return Err(format!(
                    "multiple off-torus runs on the strip boundary ({off_count} off vertices)"
                ));
            }
            run_start = Some(i);
        }
    }
    let run_start = run_start.ok_or("all boundary vertices off the torus")?;
    let run: Vec<usize> = (0..off_count).map(|k| (run_start + k) % n).collect();
    if off_count != 4 || run.iter().any(|&i| entries[i].4) {
        return Err(format!(
            "off-torus run is not the 4-vertex mitre signature ({off_count} off vertices)"
        ));
    }
    let idx_prev = (run_start + n - 1) % n;
    let idx_after = (run_start + 4) % n;
    let v_b_pos = entries[idx_prev].3; // on-torus trim end, side B
    let v_b2_pos = entries[idx_after].3; // on-torus trim end, side X
    let corner_b = entries[run[0]].2;
    let corner_b_pos = entries[run[0]].3;
    let over_b = entries[run[1]].2;
    let over_x = entries[run[2]].2;
    let corner_x = entries[run[3]].2;
    let corner_x_pos = entries[run[3]].3;
    // entries[i].1 is the edge of the half-edge LEAVING entries[i].2.
    let chord_edge = entries[run[1]].1; // over_b → over_x
    let stub_b_edge = entries[run[0]].1; // corner_b → over_b
    let stub_x_edge = entries[run[2]].1; // over_x → corner_x

    // ── End frames (derived from the torus alone) ───────────────
    let end_frame = |v_pos: Point3| -> Result<(Vec3, Point3, Vec3), String> {
        let w = v_pos - c0;
        let axial = w.dot(&a);
        let u = (w - a * axial)
            .normalized()
            .ok_or("trim end on the torus axis")?;
        let spine = c0 + u * major;
        let sec_n = a
            .cross(&u)
            .normalized()
            .ok_or("degenerate section normal")?;
        Ok((u, spine, sec_n))
    };
    let (_u_b, spine_b, w_b) = end_frame(v_b_pos)?;
    let (_u_x, spine_x, w_x) = end_frame(v_b2_pos)?;
    // The trim end must be a tube point of its section circle.
    for (v_pos, spine) in [(v_b_pos, spine_b), (v_b2_pos, spine_x)] {
        if ((v_pos - spine).length() - minor).abs() > ON_TORUS_BAND {
            return Err("trim end is not on its end-section circle".into());
        }
    }
    // Corners must sit on their section planes (they are the
    // surviving sharp-corner vertices of the tangent chain).
    let plane_band = tolerance.linear().max(1.0e-9) * 10.0;
    if (corner_b_pos - v_b_pos).dot(&w_b).abs() > plane_band
        || (corner_x_pos - v_b2_pos).dot(&w_x).abs() > plane_band
    {
        return Err("mitre corner off its end-section plane".into());
    }
    // Plane-contact side: corners sit ON the support-plane contact
    // plane (axial offset exactly ±minor).
    let axial_b = (corner_b_pos - c0).dot(&a);
    if (axial_b.abs() - minor).abs() > plane_band {
        return Err(format!(
            "corner axial offset {axial_b:.3e} is not ±minor — not the plane-side chord class"
        ));
    }
    let s_a = axial_b.signum();
    let vt_b_pos = spine_b + a * (s_a * minor);
    let vt_x_pos = spine_x + a * (s_a * minor);
    let q_c = c0 + a * (s_a * minor); // plane-contact circle center
    let old_chord_mid = Point3::new(
        (corner_b_pos.x + corner_x_pos.x) * 0.5,
        (corner_b_pos.y + corner_x_pos.y) * 0.5,
        (corner_b_pos.z + corner_x_pos.z) * 0.5,
    );

    // ── Surgery ─────────────────────────────────────────────────
    // 1. Contract the overshoot stubs; the corners survive.
    for (stub, corner) in [(stub_b_edge, corner_b), (stub_x_edge, corner_x)] {
        let he0 = brep
            .edges
            .get(stub)
            .map(|e| e.half_edge)
            .ok_or("stub gone")?;
        let he = if brep.half_edges[he0].origin == corner {
            he0
        } else {
            brep.half_edges[he0].twin
        };
        if brep.half_edges[he].origin != corner {
            return Err("stub half-edge does not originate at the corner".into());
        }
        brep.kve(he)
            .map_err(|e| format!("stub contraction failed: {e}"))?;
        // kve's repair contract: curves incident to the merged
        // vertex still reference the retired overshoot position
        // (measured: the un-blended cube edge left a 0.02 open-loop
        // gap). Re-anchor every straight incident curve; a
        // non-line mismatch is outside this signature.
        rebuild_stale_lines_at(brep, corner)?;
    }
    let _ = (over_b, over_x); // handles retired by kve
    // Rebuild the chord's curve on its surviving endpoints so the
    // splits below are well-defined.
    rebuild_line(brep, chord_edge)?;

    // 2. Split the chord at both plane-contact points and land the
    // vertices exactly.
    let vt_b = split_at(brep, chord_edge, vt_b_pos, corner_b)?;
    // The piece NOT touching corner_b carries on to corner_x.
    let middle = edge_between(brep, vt_b, corner_x).ok_or("middle chord piece not found")?;
    let vt_x = split_at(brep, middle, vt_x_pos, corner_x)?;
    let middle = edge_between(brep, vt_b, vt_x).ok_or("contact-arc piece not found")?;
    rebuild_line(
        brep,
        edge_between(brep, corner_b, vt_b).ok_or("end piece B not found")?,
    )?;
    rebuild_line(
        brep,
        edge_between(brep, corner_x, vt_x).ok_or("end piece X not found")?,
    )?;

    // 3. Exact plane-contact arc for the middle piece.
    author_circle_arc(
        brep,
        middle,
        q_c,
        a,
        vt_b_pos,
        vt_x_pos,
        Some(old_chord_mid),
    )?;
    // Exact support-contact circle for the far trim (the on-torus
    // trim through v_b/v_b2): center on the axis at the trim's own
    // axial offset.
    if let Some(far_trim) = on_torus_trim_edge(brep, strip, &entries, idx_prev, idx_after) {
        let axial_far = (v_b_pos - c0).dot(&a);
        let far_center = c0 + a * axial_far;
        let mid_hint = far_trim_mid(brep, far_trim);
        author_circle_arc(brep, far_trim, far_center, a, v_b_pos, v_b2_pos, mid_hint)?;
    }

    // 4. Section caps: one mef per end.
    let span_dir_b = if (v_b2_pos - v_b_pos).dot(&w_b) >= 0.0 {
        w_b
    } else {
        -w_b
    };
    let span_dir_x = if (v_b_pos - v_b2_pos).dot(&w_x) >= 0.0 {
        w_x
    } else {
        -w_x
    };
    let shell_outward = shell_outward_witness(brep, strip, &vt_b_pos, &(a * s_a));
    for (corner, spine, span_dir, vt_pos, vb_pos) in [
        (corner_b, spine_b, span_dir_b, vt_b_pos, v_b_pos),
        (corner_x, spine_x, span_dir_x, vt_x_pos, v_b2_pos),
    ] {
        author_section_cap(
            brep,
            strip,
            corner,
            spine,
            span_dir,
            vt_pos,
            vb_pos,
            minor,
            shell_outward,
        )?;
    }

    // ── Acceptance ──────────────────────────────────────────────
    brep.validate_loop_closure()
        .map_err(|e| format!("post-cap loop closure failed: {e:?}"))?;
    for eh in touched_edges(brep, strip) {
        let dev = edge_max_deviation(brep, eh);
        if dev > 1.0e-6 {
            return Err(format!(
                "rebuilt boundary edge {eh:?} still {dev:.3e} off an adjacent surface"
            ));
        }
    }
    Ok(2)
}
