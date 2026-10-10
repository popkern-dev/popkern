//! Per-fragment sample-and-classify for the cellular pipeline —
//! extracted verbatim from `classify_fragments_cellular`'s loop body
//! (PK-496) so seeds, audits, and propagation fallbacks can invoke
//! the FULL tuned machinery (sampling shortcuts, calibration probes,
//! offset pairs, retry ladders) as one callable. Behavior is the
//! per-fragment path, unchanged; the propagation planner decides who
//! calls it.

use std::collections::HashSet;

use popkern_core::geometry::Point3;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use super::fragment_interior::{fragment_interior_sample, planar_outward_normal};
use super::fragment_sample_helpers::{
    non_shared_chord_midpoints, point_near_any_edge, point_near_any_surface, sample_fragment_point,
};
use crate::adaptive_precision::robust_point_in_solid_final;
use crate::error::PointClassification;
use crate::merge::graph::{IntersectionGraph, Operand};
use crate::merge::walk::Fragment;

/// Shared context for one `classify_fragments_cellular` run: the
/// operands, their face AABBs, the coincident-overlap face sets, and
/// the calibrated probe offsets.
pub(super) struct ClassifyCtx<'a> {
    pub graph: &'a IntersectionGraph,
    pub brep_a: &'a BRep,
    pub brep_b: &'a BRep,
    pub aabbs_a: &'a [(FaceHandle, Aabb)],
    pub aabbs_b: &'a [(FaceHandle, Aabb)],
    pub overlap_a: &'a HashSet<FaceHandle>,
    pub overlap_b: &'a HashSet<FaceHandle>,
    pub probe_eps: f64,
    pub probe_eps_coarse: f64,
    pub tolerance: &'a Tolerance,
}

/// Sample a representative point on `fragment` and classify it
/// against the OTHER operand. Returns `(sample, class, cooriented)`
/// — `cooriented` is `Some` only on the offset-pair OnBoundary path.
pub(super) fn sample_and_classify(
    ctx: &ClassifyCtx<'_>,
    fragment: &Fragment,
) -> (Point3, PointClassification, Option<bool>) {
    let ClassifyCtx {
        graph,
        brep_a,
        brep_b,
        aabbs_a,
        aabbs_b,
        overlap_a,
        overlap_b,
        probe_eps,
        probe_eps_coarse,
        tolerance,
    } = *ctx;
    // Sample a point on this fragment. Fragments with inner
    // loops or closed-rim outers need the interior-aware
    // sampler: the edge-midpoint shortcut degenerates to the
    // SEAM VERTEX on a closed rim (maximally unstable for
    // point-in-solid — measured: the bottom cap of a cylinder
    // with a star socket cut through it was dropped and the
    // tool's outside face kept), and a plain centroid of a
    // holed cap lands INSIDE THE HOLE.
    let frag_surface = match fragment.face.operand {
        Operand::A => brep_a.faces.get(fragment.face.face),
        Operand::B => brep_b.faces.get(fragment.face.face),
    }
    .and_then(|f| f.surface.as_ref())
    .map(|s| (s.geometry.clone(), s.orientation));
    // The edge-midpoint shortcut is the tuned default. On a
    // coincident-overlap face it can land ON the other body's
    // boundary EDGE (a flush cutter's flank midpoint sits on the
    // tab's corner line), where the offset-pair probes straddle
    // an edge instead of a face and the coorientation verdict
    // inverts. Resample interior-aware ONLY in that case —
    // blanket interior sampling regressed the flush union, whose
    // arc-bounded overlap fragments sample less reliably through
    // the scanline than through the midpoint shortcut.
    let on_overlap_face_early = match fragment.face.operand {
        Operand::A => overlap_a.contains(&fragment.face.face),
        Operand::B => overlap_b.contains(&fragment.face.face),
    };
    let other_for_edges = match fragment.face.operand {
        Operand::A => brep_b,
        Operand::B => brep_a,
    };
    let shortcut = sample_fragment_point(graph, &fragment.outer_loop);
    let eps_probe = probe_eps;
    // Second trigger: the shortcut point must admit a DECISIVE
    // own-body side probe (one side in, one side out). A shortcut
    // landing on the fragment's own crease edge is indecisive at
    // any offset scale — and every probe made from it, against
    // either body, is noise (measured: a flush cutter's flank
    // fragment probed from the tab's corner edge and the
    // relief-mouth membrane was kept).
    let calibratable_at = |p: &popkern_core::geometry::Point3| -> bool {
        let own_brep_ref = match fragment.face.operand {
            Operand::A => brep_a,
            Operand::B => brep_b,
        };
        let own_aabbs_ref = match fragment.face.operand {
            Operand::A => &aabbs_a,
            Operand::B => &aabbs_b,
        };
        // The probe axis the classifier will use: the planar outward
        // normal, else the chart normal at `p` (a closed-rim cap has
        // one loop vertex and no Newell normal — answering "no
        // opinion" there let a sample ON the cap's own rim through,
        // where the +ε probe grazes the wall and calibration is a
        // parity coin flip: one of two mirror-image lens fragments
        // lost its coorientation and the subtraction kept it).
        // Decisiveness is sign-free, so the chart normal's arbitrary
        // sign is harmless here.
        let Some(n0) = planar_outward_normal(own_brep_ref, fragment.face.face).or_else(|| {
            use popkern_core::geometry::surface::SurfaceEval;
            let sof = own_brep_ref
                .faces
                .get(fragment.face.face)?
                .surface
                .as_ref()?;
            let (u, v) =
                popkern_core::geometry::uv_projection::project_point_to_uv(&sof.geometry, p)?;
            sof.geometry.normal_at(u, v).normalized()
        }) else {
            return true; // no probe axis — pre-calibration path, no opinion
        };
        for e in [eps_probe, probe_eps_coarse] {
            let a =
                robust_point_in_solid_final(&(*p + n0 * e), own_brep_ref, own_aabbs_ref, tolerance);
            let b =
                robust_point_in_solid_final(&(*p - n0 * e), own_brep_ref, own_aabbs_ref, tolerance);
            if matches!(
                (a, b),
                (PointClassification::Inside, PointClassification::Outside)
                    | (PointClassification::Outside, PointClassification::Inside)
            ) {
                return true;
            }
        }
        false
    };
    let sample = if on_overlap_face_early {
        // Both candidate samplers can fail on coincident-overlap
        // fragments: the midpoint shortcut lands on crease edges
        // and the other body's corners; the scanline interior
        // sampler mis-rows on arc-heavy charts. Pick whichever
        // admits a DECISIVE own-body side probe — the tuned
        // shortcut first unless it sits in a known hazard — and
        // keep the shortcut when neither calibrates (pre-calibration
        // behavior, orientation stays unknown downstream).
        let interior = frag_surface
            .as_ref()
            .and_then(|(surf, ori)| fragment_interior_sample(graph, fragment, surf, *ori, true));
        let hazard = point_near_any_edge(other_for_edges, &shortcut, eps_probe * 10.0);
        // HOLED fragments must sample interior-aware regardless
        // of hazard: an all-shared ring fragment's shortcut falls
        // through to the loop-vertex CENTROID, which for an
        // annular region lands in the hole — a point decisively
        // calibratable against the own body yet not ON the
        // fragment at all (measured: a moat cut's tube-material
        // ring on the flush block bottom sampled the tube's void,
        // classified Outside, and the subtract kept the ring as a
        // duplicate skin).
        let ordered = if hazard || !fragment.inner_loops.is_empty() {
            [interior, Some(shortcut)]
        } else {
            [Some(shortcut), interior]
        };
        ordered
            .into_iter()
            .flatten()
            .find(|p| calibratable_at(p))
            .unwrap_or(shortcut)
    } else if let Some(p) = frag_surface
        .as_ref()
        .and_then(|(surf, ori)| fragment_interior_sample(graph, fragment, surf, *ori, false))
    {
        p
    } else {
        // A fragment CUT by the intersection (Shared edge in its
        // loop) can have every surviving boundary edge hugging
        // the other body while only its interior is clear of it
        // (measured on the fanned flush arm: the swoop web's cut
        // curve lies on the plate underside, its flush side edges
        // on the arm walls, and both ends embed in root and
        // plate — which boundary edge the walk enumerated first
        // decided whether the web survived). Vote the non-shared
        // boundary midpoints: unanimity keeps the tuned shortcut
        // (identical to the uncalibrated behavior); disagreement means the
        // boundary straddles the other body, and the one sample
        // that is decisive by construction is a point stepped OFF
        // the cut onto the fragment's own side.
        let other = match fragment.face.operand {
            Operand::A => (brep_b, &aabbs_b),
            Operand::B => (brep_a, &aabbs_a),
        };
        let cut = fragment.outer_loop.iter().any(|&h| {
            matches!(
                graph.half_edges[h].source,
                crate::merge::graph::EdgeSource::Shared { .. }
            )
        });
        let mids = if cut {
            non_shared_chord_midpoints(graph, &fragment.outer_loop, 6)
        } else {
            Vec::new()
        };
        let unanimous = mids.len() < 2 || {
            let first = robust_point_in_solid_final(&mids[0], other.0, other.1, tolerance);
            matches!(
                first,
                PointClassification::Inside | PointClassification::Outside
            ) && mids[1..]
                .iter()
                .all(|p| robust_point_in_solid_final(p, other.0, other.1, tolerance) == first)
        };
        if unanimous {
            // The shortcut itself can land ON the other body's
            // boundary when the fragment is a thin CUT sliver —
            // a bore-cap piece beyond a breached wall has only
            // TWO boundary edges (the cut chord + a rim arc), so
            // the disagreement vote above can never fire
            // (mids.len() < 2) and the first-edge shortcut sits
            // exactly on the cutting plane, a parity coin flip
            // (measured: the sliver classified Inside and the
            // subtract kept a flap outside the body). When the
            // shortcut is hazard-positioned and a non-shared
            // boundary midpoint exists, that midpoint is the
            // decisive sample — it lies on the fragment's OWN
            // original boundary, clear of the cut.
            let hazard = cut
                && !mids.is_empty()
                && point_near_any_edge(other.0, &shortcut, eps_probe * 10.0);
            if hazard {
                mids[0]
            } else if cut && mids.is_empty() {
                // ALL-SHARED outer loop off the tuned overlap
                // path: every boundary point lies ON the other
                // operand (a cut-plane disc inside a wall), so
                // both edge midpoints and the node centroid are
                // knife-edge samples — the centroid of the grind
                // disc's two rim-graze nodes IS the
                // plane×rim×seam triple-degeneracy point, and
                // its parity coin-flip dropped the disc and
                // stranded the wall rim (PK-444, spring end
                // grind). Route through the medial interior
                // sampler graph-merge uses for this class
                // (PK-341/PK-441 hardening); its decline keeps
                // the centroid, stranding the assembly into an
                // honest decline rather than a wrong body.
                crate::merge::fragment_sample_interior::interior_step_sample(
                    graph,
                    &fragment.outer_loop,
                    &fragment.face,
                    brep_a,
                    brep_b,
                )
                .unwrap_or(shortcut)
            } else if !cut && point_near_any_surface(other.0, other.1, &shortcut, eps_probe * 10.0)
            {
                // An UNCUT fragment is one class throughout, but its
                // boundary may lie ON the other body: a box face
                // tangent to a rod along its own edge put the
                // midpoint on the rod and classified the face Inside
                // while its mirror twin sampled its far edge (PK-750).
                frag_surface
                    .as_ref()
                    .and_then(|(s, o)| fragment_interior_sample(graph, fragment, s, *o, true))
                    .unwrap_or(shortcut)
            } else {
                shortcut
            }
        } else {
            frag_surface
                .as_ref()
                .and_then(|surf| {
                    super::fragment_sampling::fragment_cut_offset_sample(graph, fragment, &surf.0)
                })
                .unwrap_or(shortcut)
        }
    };

    // Classify against the OPPOSITE body only.
    let own_operand = fragment.face.operand;
    let (other_brep, other_aabbs) = match own_operand {
        Operand::A => (brep_b, &aabbs_b),
        Operand::B => (brep_a, &aabbs_a),
    };

    // OFFSET-PAIR classification: a fragment can lie exactly on
    // (or extend) the other solid's trim geometry — a socket
    // tool's wall below the part shares its column with the
    // freshly cut hole outline, so every on-face sample grazes
    // the other body's trim and single-point parity is a coin
    // flip (measured: below-plate tool-wall fragments classified
    // Inside and kept, duplicating pocket walls while the true
    // cap was dropped). Classify the pair (sample ± ε·normal)
    // instead: agreement is decisive; disagreement means the
    // fragment sits ON the other boundary.
    // Offset-pair classification is scoped to fragments whose
    // geometry can coincide with the other body's trim curves
    // (closed-rim outers and holed fragments): applying it to
    // every fragment regressed plain cube-notch subtractions
    // (fixtures that pass with single-point classification).
    let on_overlap_face = match fragment.face.operand {
        Operand::A => overlap_a.contains(&fragment.face.face),
        Operand::B => overlap_b.contains(&fragment.face.face),
    };
    let needs_offset_pair = on_overlap_face
        || !fragment.inner_loops.is_empty()
        || fragment.outer_loop.iter().any(|&h| {
            let he = &graph.half_edges[h];
            he.origin == graph.half_edges[he.twin].origin
        });
    // `cooriented`: for OnBoundary fragments with a MATERIAL-TRUE
    // outward normal available, whether the other body's outward
    // normal at the contact agrees with this fragment's. Contact
    // orientation decides keep-one-skin semantics per operation:
    // a COORIENTED subtract contact means the removed region sits
    // exactly behind the fragment (the relief mouth is open — the
    // piece must go), while an ANTI-oriented contact is a tool
    // pressed against the wall from outside (the wall stays).
    // `None` = orientation unknown (no planar outward normal);
    // the uncalibrated rules apply.
    let mut cooriented: Option<bool> = None;
    let class = if !needs_offset_pair {
        robust_point_in_solid_final(&sample, other_brep, other_aabbs, tolerance)
    } else {
        use popkern_core::geometry::surface::SurfaceEval;
        let own_brep = match fragment.face.operand {
            Operand::A => brep_a,
            Operand::B => brep_b,
        };
        let frag_face = own_brep.faces.get(fragment.face.face);
        // Prefer the material-true outward normal (the stored
        // outer-loop Newell vector points INTO the material in
        // this kernel; negate it). The chart normal's sign is
        // arbitrary relative to material, which is fine for the
        // agree/disagree verdict but useless for coorientation.
        let outward = planar_outward_normal(own_brep, fragment.face.face);
        let normal = outward.or_else(|| {
            frag_face
                .and_then(|f| f.surface.as_ref())
                .and_then(|sof| {
                    popkern_core::geometry::uv_projection::project_point_to_uv(
                        &sof.geometry,
                        &sample,
                    )
                    .map(|(u, v)| sof.geometry.normal_at(u, v))
                })
                .and_then(|n| n.normalized())
        });
        match normal {
            Some(n_raw) => {
                let eps = probe_eps;
                // SELF-CALIBRATE the normal against the fragment's
                // OWN body: the loop-Newell sign convention is not
                // uniform across construction paths (measured: a
                // unit-cube face's Newell points into the material,
                // an extruded prism's points out of it). Probe both
                // sides against the own solid; outward is the side
                // that reads Outside. An ambiguous self-probe keeps
                // orientation unknown.
                let (own_brep_ref, own_aabbs) = match fragment.face.operand {
                    Operand::A => (brep_a, &aabbs_a),
                    Operand::B => (brep_b, &aabbs_b),
                };
                let calibrate = |e: f64| {
                    let s_plus = robust_point_in_solid_final(
                        &(sample + n_raw * e),
                        own_brep_ref,
                        own_aabbs,
                        tolerance,
                    );
                    let s_minus = robust_point_in_solid_final(
                        &(sample - n_raw * e),
                        own_brep_ref,
                        own_aabbs,
                        tolerance,
                    );
                    match (s_plus, s_minus) {
                        (PointClassification::Outside, PointClassification::Inside) => Some(n_raw),
                        (PointClassification::Inside, PointClassification::Outside) => {
                            Some(n_raw * -1.0)
                        }
                        _ => None,
                    }
                };
                // Same corner-noise hazard as the classification
                // pair: retry an indecisive calibration at the
                // coarser offset.
                let n_outward = calibrate(eps).or_else(|| calibrate(probe_eps_coarse));
                let n = n_outward.unwrap_or(n_raw);
                let pair_at = |e: f64| {
                    let p_plus = sample + n * e;
                    let p_minus = sample - n * e;
                    let v_plus =
                        robust_point_in_solid_final(&p_plus, other_brep, other_aabbs, tolerance);
                    let v_minus =
                        robust_point_in_solid_final(&p_minus, other_brep, other_aabbs, tolerance);
                    // the pair-scan trace topic: per-probe forensics for the
                    // offset-pair classifier — prints the exact
                    // probe points and their verdicts (the face-selection trace topic
                    // shows only the resolved pair). This is how
                    // the out-of-corridor ×16 retry was caught.
                    if popkern_core::trace::trace_enabled("pair-scan") {
                        eprintln!(
                            "        [pp] e={e:.3e} +({:.3},{:.3},{:.3})={v_plus:?} -({:.3},{:.3},{:.3})={v_minus:?}",
                            p_plus.x, p_plus.y, p_plus.z, p_minus.x, p_minus.y, p_minus.z
                        );
                    }
                    (v_plus, v_minus)
                };
                // Probe at the tuned offset AND at 16× it: 1 µm
                // probes graze corner-adjacent faces when the
                // sample sits on the other body's edge or corner
                // fold, and the ray parity inverts the verdict
                // (measured: a flush cutter's flank fragment read
                // anti-oriented at its corner sample and the
                // relief-mouth membrane was kept). When the two
                // scales disagree, a THIRD pair at 4× breaks the
                // tie by majority: blindly trusting the coarser
                // scale mis-served coarse-tolerance bodies whose
                // sample clearance to the other body sits between
                // the two offsets — at linear 1e-3 the 16× probe
                // reached 0.09, PAST the sphere face 0.053 away,
                // and the straddling verdict dropped every cube
                // annulus (measured: centered-sphere subtraction,
                // six half-claimed rim edges). All-distinct keeps
                // the coarse verdict (prior behavior).
                let (c_plus_f, c_minus_f) = pair_at(eps);
                let (c_plus, c_minus) = {
                    let fine = (c_plus_f, c_minus_f);
                    let coarse = pair_at(probe_eps_coarse);
                    if fine == coarse {
                        fine
                    } else {
                        let mid = pair_at(probe_eps_coarse.sqrt() * eps.sqrt());
                        if mid == fine { fine } else { coarse }
                    }
                };
                if popkern_core::trace::trace_enabled("face-selection") {
                    eprintln!(
                        "      [pair] {:?} outer={} n=({:.2},{:.2},{:.2}) calibrated={} eps={eps:.2e} c_plus={c_plus:?} c_minus={c_minus:?}",
                        fragment.face,
                        fragment.outer_loop.len(),
                        n.x,
                        n.y,
                        n.z,
                        n_outward.is_some()
                    );
                }
                match (c_plus, c_minus) {
                    (PointClassification::Inside, PointClassification::Inside) => {
                        PointClassification::Inside
                    }
                    (PointClassification::Outside, PointClassification::Outside) => {
                        PointClassification::Outside
                    }
                    (c_p, c_m) => {
                        if n_outward.is_some() {
                            // Other body's material on the −n side
                            // (behind this fragment's skin) means
                            // its outward normal ≈ +n: cooriented.
                            cooriented = match (c_p, c_m) {
                                (PointClassification::Outside, PointClassification::Inside) => {
                                    Some(true)
                                }
                                (PointClassification::Inside, PointClassification::Outside) => {
                                    Some(false)
                                }
                                _ => None,
                            };
                        }
                        PointClassification::OnBoundary
                    }
                }
            }
            None => robust_point_in_solid_final(&sample, other_brep, other_aabbs, tolerance),
        }
    };
    (sample, class, cooriented)
}
