//! Fillet strip end-boundary reconstruction on planar end faces.
//!
//! A cylinder fillet strip that terminates on a planar end face (e.g.
//! a vertical corner fillet running out on the part's base) must end
//! in the exact cylinder ∩ plane arc between its two wall tangency
//! feet. The carve phase instead leaves the ORIGINAL corner vertex —
//! at r·√2 from the strip axis, r(√2−1) off the cylinder — joined to
//! the feet by straight chords (an "L"), and earlier repair passes
//! either decline on that off-surface anchor (`edge_reconcile`
//! `AnchorOffSurface`) or leave non-analytic replacement curves, one
//! of which was measured marching 270° the long way around the strip
//! (flagship PK-099: bottom-corner render artifacts, off-surface
//! residual 3.3e-2, u-wrap tessellation declines, spur slits in the
//! base loop).
//!
//! This pass is the strip-end analogue of the route B corner-cap
//! surgery (`corner_arc_boundaries`):
//!
//! 1. contract each off-cylinder end corner into an on-cylinder foot
//!    via `kve` (the L becomes a single foot→foot edge),
//! 2. replace every non-analytic strip/plane runout boundary between
//!    on-cylinder feet with the exact verified cylinder ∩ plane Arc
//!    (short sweep; angles in the Arc evaluator's own `circle_basis`
//!    frame, both normal candidates tried, 1e-9 round-trip
//!    acceptance),
//! 3. remove zero-width spur slits left in the end face's loop by
//!    earlier rebuilds (`kev` on half-edges immediately followed by
//!    their twin).
//!
//! Containment boundary (recorded PK-097 Phase C, retired for the
//! plain case by PK-148): this pass only runs when the blend
//! produced corner patches (chain/blend context) — the caller gates
//! on that. Plain open cylinder fillets are settled by the
//! zero-topology `terminal_arc_trim` pass instead (corner vertex
//! repositioned onto the arc, stub curves re-authored), which the
//! in-place modify machinery also runs so in-place stays equivalent
//! to remove-recreate — the equivalence concern that originally
//! froze the legacy end shape.
//!
//! Per-corner and per-slit transactional: each mutation is applied on
//! a candidate clone and accepted only when `validate()` holds;
//! declines keep the existing topology and report via PK_CARVE_TRACE.
//! Tolerances: the 5e-3 anchor-placement gate (same class as route
//! B's snap gate) for on/off-cylinder classification; 1e-9 for arc
//! acceptance. Both are local to this pass.

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle, HalfEdgeHandle};

/// Anchor-placement gate: how far a vertex may sit from the cylinder
/// (radially) and still count as an on-surface foot; and how far off
/// it must be to count as a carve-artifact corner.
const SNAP_GATE: f64 = 5.0e-3;

/// Minimum |axis·n̂| for a closed end section: below this the end
/// plane is close to parallel with the strip axis and the section
/// degenerates toward generator lines — decline, recorded.
const OBLIQUE_SECTION_MIN_COS: f64 = 0.1;

/// Reconstruct strip end boundaries. Returns (contracted corners,
/// arcs authored, spurs removed).
use super::strip_end_surface::{StripSurface, fill_missing_planar_pcurves, remove_spur_slits};

pub(super) fn reconstruct_strip_end_boundaries(
    brep: &mut BRep,
    blend_surfaces: &[SurfaceGeometry],
) -> (usize, usize, usize) {
    let tol = popkern_core::tolerance::Tolerance::default();
    let strips: Vec<(FaceHandle, bool)> = brep
        .faces
        .iter()
        .filter_map(|(h, f)| match f.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Cylinder { .. }) => Some((h, true)),
            // Planar strips: ONLY this operation's own authored flat
            // bevels (the surface_fix trust rule). Arbitrary planar
            // faces must not be contracted — Gregory support
            // sub-faces legitimately carry corner anchors r·(√2−1)
            // off their plane, owned by the corner passes (measured:
            // unrestricted planar contraction broke the Gregory
            // STEP-fidelity classification).
            Some(&SurfaceGeometry::Plane {
                origin: fo,
                normal: fn_,
                ..
            }) => blend_surfaces
                .iter()
                .any(|bs| match bs {
                    SurfaceGeometry::Plane {
                        origin: bo,
                        normal: bn,
                        ..
                    } => {
                        let (Some(na), Some(nb)) = (fn_.normalized(), bn.normalized()) else {
                            return false;
                        };
                        // Parallel within the angular tolerance and
                        // coplanar within the contraction SNAP_GATE
                        // (post-pipeline passes may drift the chart
                        // origin slightly off the authored plane).
                        (na.dot(&nb).abs() - 1.0).abs() < tol.angular()
                            && ((fo - *bo).dot(&na)).abs() <= SNAP_GATE
                    }
                    _ => false,
                })
                .then_some((h, false)),
            _ => None,
        })
        .collect();

    let mut contracted = 0;
    let mut arcs = 0;
    let mut touched_planes: Vec<FaceHandle> = Vec::new();

    let mut contracted_feet: Vec<popkern_core::topology::VertexHandle> = Vec::new();
    for &(strip, is_cylinder) in &strips {
        contracted += contract_end_corners(brep, strip, &mut touched_planes, &mut contracted_feet);
        if is_cylinder {
            arcs += arcify_runout_boundaries(brep, strip, &mut touched_planes);
        }
    }

    touched_planes.sort();
    touched_planes.dedup();
    let mut spurs = 0;
    for plane in &touched_planes {
        spurs += remove_spur_slits(brep, *plane);
    }
    // Contractions merge edges via kve and the merged half-edges
    // carry no pcurves; fill the gaps SURGICALLY — only faces
    // incident to the surviving contraction feet. MISSING pcurves
    // only (authored analytic pcurves are never overwritten), and
    // the scope stays inside this pass's own surgery: a broad fill
    // was measured to heal pcurve gaps that Gregory repair
    // diagnostics deliberately probe.
    if !contracted_feet.is_empty() {
        let mut affected: Vec<FaceHandle> = Vec::new();
        for (_, l) in brep.loops.iter() {
            let start = l.half_edge;
            let mut he = start;
            for _ in 0..brep.half_edges.len() + 1 {
                let Some(d) = brep.half_edges.get(he) else {
                    break;
                };
                if contracted_feet.contains(&d.origin) {
                    affected.push(l.face);
                    break;
                }
                he = d.next;
                if he == start {
                    break;
                }
            }
        }
        affected.sort();
        affected.dedup();
        for face in affected {
            fill_missing_planar_pcurves(brep, face);
        }
    }
    (contracted, arcs, spurs)
}

fn cylinder_of(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3, f64)> {
    match brep.faces.get(face)?.surface.as_ref()?.geometry {
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
            ..
        } => Some((origin, axis, radius)),
        _ => None,
    }
}

fn plane_of(brep: &BRep, face: FaceHandle) -> Option<(Point3, Vec3)> {
    match brep.faces.get(face)?.surface.as_ref()?.geometry {
        SurfaceGeometry::Plane { origin, normal, .. } => Some((origin, normal)),
        _ => None,
    }
}

fn face_across(brep: &BRep, he: HalfEdgeHandle) -> Option<FaceHandle> {
    let twin = brep.half_edges.get(he)?.twin;
    let l = brep.half_edges.get(twin)?.loop_;
    Some(brep.loops.get(l)?.face)
}

pub(super) fn radial_dist(origin: &Point3, axis: &Vec3, p: &Point3) -> f64 {
    let d = *p - *origin;
    let ax = d.dot(axis);
    (d - *axis * ax).length()
}

fn is_chord(curve: &CurveOnEdge) -> bool {
    match &curve.geometry {
        CurveGeometry::Line { .. } => true,
        CurveGeometry::Nurbs(n) => n.degree == 1 && n.control_points.len() == 2,
        _ => false,
    }
}

/// Phase 1: contract each off-cylinder end corner (both flanking
/// chord edges shared with the SAME plane face) into its on-cylinder
/// predecessor foot via `kve`.
fn contract_end_corners(
    brep: &mut BRep,
    strip: FaceHandle,
    touched_planes: &mut Vec<FaceHandle>,
    contracted_feet: &mut Vec<popkern_core::topology::VertexHandle>,
) -> usize {
    let mut total = 0;
    // Restart the walk after each accepted contraction (loop
    // structure changes).
    for _ in 0..16 {
        let Some(surface) = StripSurface::of(brep, strip) else {
            return total;
        };
        let Some(face) = brep.faces.get(strip) else {
            return total;
        };
        let start_he = brep.loops[face.outer_loop].half_edge;
        let mut target: Option<(HalfEdgeHandle, FaceHandle)> = None;
        let mut he = start_he;
        for _ in 0..brep.half_edges.len() + 1 {
            let d = &brep.half_edges[he];
            let next = d.next;
            // Corner = destination of `he`, flanked by `he` (foot→corner)
            // and `next` (corner→foot).
            let corner = brep.half_edges[next].origin;
            let foot_a = d.origin;
            let foot_b = brep.half_edges[brep.half_edges[next].next].origin;
            let p_corner = brep.vertices[corner].position;
            let p_a = brep.vertices[foot_a].position;
            let p_b = brep.vertices[foot_b].position;
            let across_in = face_across(brep, he);
            let across_out = face_across(brep, next);
            let plane_face = match (across_in, across_out) {
                (Some(f1), Some(f2)) if f1 == f2 && f1 != strip => f1,
                _ => {
                    he = next;
                    if he == start_he {
                        break;
                    }
                    continue;
                }
            };
            let both_chords = [he, next].iter().all(|&h| {
                let e = brep.half_edges[h].edge;
                brep.edges[e].curve.as_ref().map(is_chord).unwrap_or(false)
            });
            if both_chords
                && plane_of(brep, plane_face).is_some()
                && surface.excess(&p_corner) > SNAP_GATE
                && surface.deviation(&p_a) <= SNAP_GATE
                && surface.deviation(&p_b) <= SNAP_GATE
            {
                target = Some((he, plane_face));
                break;
            }
            he = next;
            if he == start_he {
                break;
            }
        }
        let Some((he, plane_face)) = target else {
            return total;
        };
        let mut candidate = brep.clone();
        match candidate.kve(he) {
            Ok(r) => {
                contracted_feet.push(r.vertex);
                super::corner_arc_boundaries::rebuild_incident_straight_edges(
                    &mut candidate,
                    r.vertex,
                );
                if candidate.validate().is_ok() {
                    *brep = candidate;
                    touched_planes.push(plane_face);
                    total += 1;
                    if std::env::var_os("PK_CARVE_TRACE").is_some() {
                        eprintln!("  [strip-end] {strip:?}: contracted end corner into foot");
                    }
                } else {
                    if std::env::var_os("PK_CARVE_TRACE").is_some() {
                        eprintln!("  [strip-end] {strip:?}: contraction ROLLED BACK");
                    }
                    return total;
                }
            }
            Err(e) => {
                if std::env::var_os("PK_CARVE_TRACE").is_some() {
                    eprintln!("  [strip-end] {strip:?}: kve declined: {e:?}");
                }
                return total;
            }
        }
    }
    total
}

/// Phase 2: replace non-analytic strip/plane runout boundaries
/// (angle-spanning edges between on-cylinder feet) with exact
/// verified cylinder ∩ plane arcs.
fn arcify_runout_boundaries(
    brep: &mut BRep,
    strip: FaceHandle,
    touched_planes: &mut Vec<FaceHandle>,
) -> usize {
    let Some((origin, axis, radius)) = cylinder_of(brep, strip) else {
        return 0;
    };
    let Some(face) = brep.faces.get(strip) else {
        return 0;
    };
    let start_he = brep.loops[face.outer_loop].half_edge;
    let mut runouts: Vec<(popkern_core::topology::EdgeHandle, FaceHandle)> = Vec::new();
    let mut he = start_he;
    for _ in 0..brep.half_edges.len() + 1 {
        let d = &brep.half_edges[he];
        let next = d.next;
        let eh = d.edge;
        if let Some(across) = face_across(brep, he)
            && across != strip
            && plane_of(brep, across).is_some()
            && brep.edges[eh]
                .curve
                .as_ref()
                .map(|c| !matches!(c.geometry, CurveGeometry::Arc { .. }))
                .unwrap_or(false)
        {
            let v_o = d.origin;
            let v_d = brep.half_edges[d.twin].origin;
            if v_o != v_d {
                let p_o = brep.vertices[v_o].position;
                let p_d = brep.vertices[v_d].position;
                let on_cyl =
                    |p: &Point3| (radial_dist(&origin, &axis, p) - radius).abs() <= SNAP_GATE;
                // Runout edges SPAN angle; wall trim lines don't.
                let ang = |p: &Point3| {
                    let (f1, f2) =
                        popkern_core::geometry::analytic_surfaces::orthonormal_basis(axis);
                    let d = *p - origin;
                    d.dot(&f2).atan2(d.dot(&f1))
                };
                let mut da = (ang(&p_o) - ang(&p_d)).abs();
                if da > std::f64::consts::PI {
                    da = std::f64::consts::TAU - da;
                }
                if on_cyl(&p_o) && on_cyl(&p_d) && da > 0.1 && da < std::f64::consts::PI - 0.1 {
                    runouts.push((eh, across));
                }
            }
        }
        he = next;
        if he == start_he {
            break;
        }
    }

    let mut rewritten = 0;
    for (eh, plane_face) in runouts {
        let Some((p_origin, p_normal)) = plane_of(brep, plane_face) else {
            continue;
        };
        let Some(n_hat) = p_normal.normalized() else {
            continue;
        };
        let denom = axis.dot(&n_hat);
        // A near-parallel end plane has no closed section (the
        // intersection degenerates toward generator lines) — decline,
        // recorded. Perpendicular ends author the exact circle;
        // oblique ends the exact ellipse (PK-262 — the former
        // recorded decline whose kept incoherent runout curve folded
        // the strip's trim polygon).
        if denom.abs() < OBLIQUE_SECTION_MIN_COS {
            if std::env::var_os("PK_CARVE_TRACE").is_some() {
                eprintln!(
                    "  [strip-end] {strip:?} edge {eh:?}: near-parallel end (|axis·n|={:.4}) — declined",
                    denom.abs()
                );
            }
            continue;
        }
        // Section center: strip axis pierced through the end plane.
        let t = (p_origin - origin).dot(&n_hat) / denom;
        let c0 = origin + axis * t;
        let he0 = brep.edges[eh].half_edge;
        let (v1, v2) = {
            let d = &brep.half_edges[he0];
            (d.origin, brep.half_edges[d.twin].origin)
        };
        let p1 = brep.vertices[v1].position;
        let p2 = brep.vertices[v2].position;
        let section = if denom.abs() >= 1.0 - 1.0e-9 {
            super::strip_end_section_curves::verified_circle_section(
                &origin, &axis, radius, c0, p1, p2, &p_origin, &n_hat,
            )
        } else {
            super::strip_end_section_curves::verified_ellipse_section(
                &origin, &axis, radius, c0, p1, p2, &p_origin, &n_hat, denom,
            )
        };
        let Some((new_curve, s1, s2)) = section else {
            if std::env::var_os("PK_CARVE_TRACE").is_some() {
                eprintln!("  [strip-end] {strip:?} edge {eh:?}: no verified section candidate");
            }
            continue;
        };
        let mut candidate = brep.clone();
        crate::construct::edge_reconcile::snap_vertex_and_adjacent_chords(
            &mut candidate,
            v1,
            s1,
            eh,
        );
        crate::construct::edge_reconcile::snap_vertex_and_adjacent_chords(
            &mut candidate,
            v2,
            s2,
            eh,
        );
        if let Some(e) = candidate.edges.get_mut(eh) {
            e.curve = Some(new_curve);
            // Stale pcurves resurrect the old geometry in later fix
            // passes — null both sides, recompute re-authors from
            // the arc.
            let he0 = e.half_edge;
            let tw = candidate.half_edges[he0].twin;
            candidate.half_edges[he0].pcurve = None;
            candidate.half_edges[tw].pcurve = None;
        }
        if candidate.validate().is_ok() {
            *brep = candidate;
            touched_planes.push(plane_face);
            rewritten += 1;
            if std::env::var_os("PK_CARVE_TRACE").is_some() {
                eprintln!("  [strip-end] {strip:?} edge {eh:?}: exact runout arc authored");
            }
        } else if std::env::var_os("PK_CARVE_TRACE").is_some() {
            eprintln!("  [strip-end] {strip:?} edge {eh:?}: arc ROLLED BACK (validate failed)");
        }
    }
    rewritten
}
