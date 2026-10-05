//! Trim-residue post-filters for Stage 1 records (split from
//! `intersect.rs` under the 500-line structural budget): the
//! op-aware kept-flank drop, boundary/redundancy honesty gates, and
//! the residue-swarm cap (PK-465/PK-472).

use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::topology::FaceHandle;

use super::intersect::IntersectionRecord;

/// Drop trim-residue records whose flanking regions would be KEPT
/// under `op` (PK-465). A residue's two flanks always classify
/// identically (crossing the other operand's surface EXTENSION
/// changes nothing); when that shared verdict keeps them, the
/// residue is an informationless slit inside a kept region and must
/// not partition the face (measured: a kept deck walked around such
/// a slit and its assembled loop closed with a 1.2 gap). The verdict
/// is measured at the residue midpoint against the other operand —
/// a point on the surface extension BEYOND the other face's trim,
/// so the ray verdict is not knife-edge.
pub(crate) fn drop_kept_flank_residues(
    records: &mut Vec<IntersectionRecord>,
    op: crate::BoolOperation,
    brep_a: &popkern_core::topology::BRep,
    brep_b: &popkern_core::topology::BRep,
    aabbs_a: &[(FaceHandle, Aabb)],
    aabbs_b: &[(FaceHandle, Aabb)],
    tolerance: &popkern_core::tolerance::Tolerance,
) {
    use crate::error::PointClassification;
    let trace = popkern_core::trace::trace_enabled("boolean");
    // One verdict per residue holds only between the other operand's
    // boundary crossings: cut residues where it changes first (PK-553).
    crate::intersect_residue_split::split_residues_at_verdict_changes(
        records,
        &crate::intersect_residue_split::OtherOperand {
            brep_a,
            brep_b,
            aabbs_a,
            aabbs_b,
        },
        tolerance,
    );
    // Chord samples + bounding boxes of every ORDINARY record, for
    // the cross-record redundancy gate below.
    let ordinary: Vec<(
        popkern_core::geometry::Point3,
        popkern_core::geometry::Point3,
        Vec<popkern_core::geometry::Point3>,
    )> = records
        .iter()
        .filter(|r| r.residue_owner_a.is_none())
        .map(|r| {
            use popkern_core::geometry::curve::CurveEval as _;
            let (t0, t1) = r.curve.curve_3d.parameter_range();
            let pts: Vec<popkern_core::geometry::Point3> = (0..=16)
                .map(|k| r.curve.curve_3d.point_at(t0 + (t1 - t0) * k as f64 / 16.0))
                .collect();
            let mut lo = popkern_core::geometry::Point3::new(f64::MAX, f64::MAX, f64::MAX);
            let mut hi = popkern_core::geometry::Point3::new(f64::MIN, f64::MIN, f64::MIN);
            for p in &pts {
                lo = popkern_core::geometry::Point3::new(
                    lo.x.min(p.x),
                    lo.y.min(p.y),
                    lo.z.min(p.z),
                );
                hi = popkern_core::geometry::Point3::new(
                    hi.x.max(p.x),
                    hi.y.max(p.y),
                    hi.z.max(p.z),
                );
            }
            (lo, hi, pts)
        })
        .collect();
    records.retain(|r| {
        let Some(owner_is_a) = r.residue_owner_a else {
            return true;
        };
        use popkern_core::geometry::curve::CurveEval as _;
        let (t0, t1) = r.curve.curve_3d.parameter_range();
        let mid = r.curve.curve_3d.point_at(0.5 * (t0 + t1));
        // EXPLICIT boundary proximity first: a residue riding the
        // other operand's boundary is the flush class the OnBoundary
        // arm below drops — but multi-ray parity ON a boundary point
        // is a coin flip, and a flush residue voting Inside under a
        // union was RETAINED (measured: the flush gusset's residues
        // ride the arm's bottom edges; the retained slits dead-ended
        // and the loop walk broke with a 29.6 gap). Distance to the
        // other operand's faces is the honest instrument.
        let (other_brep, other_aabbs) = if owner_is_a {
            (brep_b, aabbs_b)
        } else {
            (brep_a, aabbs_a)
        };
        let boundary_gate = tolerance.linear() * 100.0;
        // CROSS-RECORD REDUNDANCY: a residue whose span is already
        // covered by an ORDINARY record's curve partitions nothing —
        // the partition exists as a real two-face edge (measured,
        // PK-472: the impeller's flow-cavity star walls each emit
        // the overhang of their neighbour wall's crossing as a
        // residue — 594 emissions, 346 surviving — and the walk
        // wove the duplicates into the outline loops, breaking a
        // half-edge cycle on a baseline-served boolean). The
        // load-bearing classes (rod over-rim, m6 crest circle, tap
        // in-hole flap) have NO ordinary record along them.
        {
            let pad = boundary_gate;
            let covered = ordinary.iter().any(|(lo, hi, pts)| {
                let inside_box = mid.x >= lo.x - pad
                    && mid.x <= hi.x + pad
                    && mid.y >= lo.y - pad
                    && mid.y <= hi.y + pad
                    && mid.z >= lo.z - pad
                    && mid.z <= hi.z + pad;
                inside_box
                    && pts.windows(2).any(|w| {
                        let ab = w[1] - w[0];
                        let len2 = ab.dot(&ab);
                        let t = if len2 > 0.0 {
                            ((mid - w[0]).dot(&ab) / len2).clamp(0.0, 1.0)
                        } else {
                            0.0
                        };
                        let foot = popkern_core::geometry::Point3::new(
                            w[0].x + ab.x * t,
                            w[0].y + ab.y * t,
                            w[0].z + ab.z * t,
                        );
                        // Polyline distance, widened by the chord
                        // sag of 1/16 sampling (sag ≈ seg²·κ/8).
                        foot.distance_to(&mid) < pad * 4.0 + len2.sqrt() * 0.1
                    })
            });
            if covered {
                if trace {
                    eprintln!(
                        "  [residue-drop] redundant residue at ({:.3},{:.3},{:.3}) removed (covered by an ordinary record)",
                        mid.x, mid.y, mid.z
                    );
                }
                return false;
            }
        }

        // A residue riding its OWNER face's OWN trim boundary
        // duplicates that boundary and partitions nothing — the walk
        // weaves the duplicate pair into the boundary loop (measured,
        // PK-472: the impeller's flow-cavity floor emitted residues
        // ALONG its 406-point star outline; the assembled loops
        // interleaved original and duplicate edges at identical
        // positions and the boolean declined with a broken half-edge
        // cycle). The load-bearing classes cross their owner's
        // INTERIOR — a slit tracing the owner's own boundary never
        // separates regions.
        {
            let (own_brep, own_face) = if owner_is_a {
                (brep_a, r.face_a)
            } else {
                (brep_b, r.face_b)
            };
            let mut min_d = f64::MAX;
            for eh in crate::clip_endpoint_refine::boundary_edges(own_brep, own_face) {
                let Some(c) = own_brep.edges.get(eh).and_then(|e| e.curve.as_ref()) else {
                    continue;
                };
                if let Some(t) = (&c.geometry as &dyn popkern_core::geometry::curve::CurveEval)
                    .closest_point(&mid)
                {
                    let tc = t.clamp(c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                    min_d = min_d.min(c.geometry.point_at(tc).distance_to(&mid));
                }
            }
            if min_d < boundary_gate {
                if trace {
                    eprintln!(
                        "  [residue-drop] own-boundary residue at ({:.3},{:.3},{:.3}) removed (d {min_d:.2e})",
                        mid.x, mid.y, mid.z
                    );
                }
                return false;
            }
        }
        let on_other_boundary = crate::intersect_residue_split::on_operand_boundary(
            &mid,
            other_brep,
            other_aabbs,
            boundary_gate,
            tolerance,
        );
        if on_other_boundary {
            if trace {
                eprintln!(
                    "  [residue-drop] flush residue at ({:.3},{:.3},{:.3}) removed (on other operand's boundary)",
                    mid.x, mid.y, mid.z
                );
            }
            return false;
        }
        // The other operand is a FINAL body (a boolean input, never a
        // mid-split state): classify with the final-body gates — hole,
        // outer-loop and cylinder/cone trim — so a ray from a void
        // (a bore centre) cannot count phantom hits on a face's
        // surface extension (PK-672).
        let verdict = if owner_is_a {
            crate::classify::point_in_solid_final(&mid, brep_b, aabbs_b, tolerance)
        } else {
            crate::classify::point_in_solid_final(&mid, brep_a, aabbs_a, tolerance)
        };
        let keeps = match verdict {
            Ok(PointClassification::Inside) => match (op, owner_is_a) {
                (crate::BoolOperation::Subtraction, true) => false,
                (crate::BoolOperation::Subtraction, false) => true,
                (crate::BoolOperation::Union, _) => false,
                (crate::BoolOperation::Intersection, _) => true,
            },
            Ok(PointClassification::Outside) => match (op, owner_is_a) {
                (crate::BoolOperation::Subtraction, true) => true,
                (crate::BoolOperation::Subtraction, false) => false,
                (crate::BoolOperation::Union, _) => true,
                (crate::BoolOperation::Intersection, _) => false,
            },
            // Boundary / inconclusive: DROP the residue. A residue
            // riding the other operand's boundary is the
            // flush/coincident class, whose topology the boundary
            // machinery already owns — doubling it as a slit
            // self-twins the serve (measured: the half-wiring
            // fixture's deck residue rode the cutter wall plane and
            // its verdict came back OnBoundary). Load-bearing
            // residues sit strictly off the other operand.
            _ => true,
        };
        if keeps && trace {
            eprintln!(
                "  [residue-drop] kept-flank residue at ({:.3},{:.3},{:.3}) removed",
                mid.x, mid.y, mid.z
            );
        }
        if !keeps && trace {
            eprintln!(
                "  [residue-keep] A{:?} x B{:?} owner={} at ({:.3},{:.3},{:.3})",
                r.face_a,
                r.face_b,
                if owner_is_a { "A" } else { "B" },
                mid.x,
                mid.y,
                mid.z
            );
        }
        !keeps
    });

    dedupe_same_face_residues(records, tolerance);

    // RESIDUE-SWARM CAP — temporary containment (PK-472, follow-up
    // filed): a single owner face accumulating a SWARM of surviving
    // residues is the many-face-imprint shape (measured: the
    // impeller's flow-cavity floor collected 26+ owner-side residues,
    // one per star wall, radiating inside an imprint outline already
    // fully covered by ordinary records — the walk wove the swarm
    // into the outline loops and a baseline-served boolean declined
    // with a broken half-edge cycle). The measured load-bearing
    // classes stay far under the cap (tap ≤ 8 per face, m6/rod ≤ 2).
    // The commercial-grade replacement is a per-residue separation
    // witness (does this slit separate a kept region from a
    // discarded one), which needs region knowledge this stage does
    // not have.
    {
        let mut per_face: std::collections::HashMap<(bool, FaceHandle), usize> =
            std::collections::HashMap::new();
        for r in records.iter() {
            if let Some(owner_is_a) = r.residue_owner_a {
                let f = if owner_is_a { r.face_a } else { r.face_b };
                *per_face.entry((owner_is_a, f)).or_insert(0) += 1;
            }
        }
        const RESIDUE_SWARM_CAP: usize = 16;
        let swarmed: Vec<(bool, FaceHandle)> = per_face
            .iter()
            .filter(|&(_, &c)| c > RESIDUE_SWARM_CAP)
            .map(|(&k, _)| k)
            .collect();
        if !swarmed.is_empty() {
            if trace {
                for (a, f) in &swarmed {
                    eprintln!(
                        "  [residue-drop] swarm cap: {} residues on {} face {f:?} exceed {RESIDUE_SWARM_CAP} — all dropped (many-face imprint shape)",
                        per_face[&(*a, *f)],
                        if *a { "A" } else { "B" },
                    );
                }
            }
            records.retain(|r| {
                let Some(owner_is_a) = r.residue_owner_a else {
                    return true;
                };
                let f = if owner_is_a { r.face_a } else { r.face_b };
                !swarmed.contains(&(owner_is_a, f))
            });
        }
    }
}

/// CROSS-PAIR residue dedupe (PK-487): different face pairs can each
/// emit a residue of the SAME span onto the SAME owner face — every
/// pair whose windows overlap the stretch subtracts the shared claim
/// independently (measured on the rotated torus × torus union while
/// its containment vote was inverted: one lens span survived as SIX
/// residues across five pairs, four of them partitioning one A-flank;
/// the duplicate single-face edges collided in the graph and the
/// assembly stranded half-claimed sides). A residue partitions only
/// its owner face, so one copy per owner face per span carries the
/// whole partition; a duplicate whose start, mid and end all lie on a
/// surviving same-owner-face residue's polyline (within the marching
/// slop gate) and whose extent is not larger drops. Ties keep the
/// earlier record — deterministic under the intersect stage's stable
/// record order.
fn dedupe_same_face_residues(
    records: &mut Vec<IntersectionRecord>,
    tolerance: &popkern_core::tolerance::Tolerance,
) {
    use popkern_core::geometry::curve::CurveEval as _;
    let trace = popkern_core::trace::trace_enabled("boolean");
    let residue_idx: Vec<usize> = (0..records.len())
        .filter(|&i| records[i].residue_owner_a.is_some())
        .collect();
    if residue_idx.len() < 2 {
        return;
    }
    let owner_face = |r: &IntersectionRecord| -> (bool, FaceHandle) {
        let a = r.residue_owner_a == Some(true);
        (a, if a { r.face_a } else { r.face_b })
    };
    let sample = |r: &IntersectionRecord, f: f64| -> popkern_core::geometry::Point3 {
        let c = &r.curve.curve_3d;
        let (t0, t1) = c.parameter_range();
        c.point_at(t0 + (t1 - t0) * f)
    };
    let extent = |r: &IntersectionRecord| -> f64 {
        sample(r, 0.0)
            .distance_to(&sample(r, 1.0))
            .max(sample(r, 0.0).distance_to(&sample(r, 0.5)))
    };
    // The same endpoint-slop class the consolidation gates use:
    // independent traces of one stretch land within the marching
    // slop, far above linear tolerance and far below any span's
    // extent.
    let mut drop = vec![false; records.len()];
    for (ai, &ia) in residue_idx.iter().enumerate() {
        if drop[ia] {
            continue;
        }
        let gate = (0.05 * extent(&records[ia])).max(tolerance.linear() * 10.0);
        let poly: Vec<popkern_core::geometry::Point3> = (0..=16)
            .map(|k| sample(&records[ia], k as f64 / 16.0))
            .collect();
        let on_poly = |p: &popkern_core::geometry::Point3| -> bool {
            poly.windows(2).any(|w| {
                let ab = w[1] - w[0];
                let len2 = ab.dot(&ab);
                let t = if len2 > 1e-300 {
                    ((*p - w[0]).dot(&ab) / len2).clamp(0.0, 1.0)
                } else {
                    0.0
                };
                (w[0] + ab * t).distance_to(p) < gate
            })
        };
        for &ib in residue_idx.iter().skip(ai + 1) {
            if drop[ib]
                || owner_face(&records[ia]) != owner_face(&records[ib])
                || extent(&records[ib]) > extent(&records[ia]) + gate
            {
                continue;
            }
            if on_poly(&sample(&records[ib], 0.0))
                && on_poly(&sample(&records[ib], 0.5))
                && on_poly(&sample(&records[ib], 1.0))
            {
                drop[ib] = true;
                if trace {
                    let m = sample(&records[ib], 0.5);
                    eprintln!(
                        "  [residue-drop] cross-pair duplicate at ({:.3},{:.3},{:.3}) removed \
                         (same owner face, covered by pair ({:?},{:?}))",
                        m.x, m.y, m.z, records[ia].face_a, records[ia].face_b
                    );
                }
            }
        }
    }
    let mut it = drop.into_iter();
    records.retain(|_| !it.next().unwrap_or(false));
}

#[cfg(test)]
#[path = "intersect_residue_filter_tests.rs"]
mod tests;
