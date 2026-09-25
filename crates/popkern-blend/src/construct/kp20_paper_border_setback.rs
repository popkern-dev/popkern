//! Corner-setback border extraction (PK-090 setback path, PK-398
//! exact section rows). Nested in `kp20_paper_border_extract` to share
//! its segment/calibration substrate; split out for the file budget.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::gregory::{Kp20PaperGuideBorder, Kp20PaperGuideBorderPiece};
use popkern_core::geometry::{Point3, Vec3};

use super::section_path::SectionPath;
use super::sides::{
    SetbackJunctionPolicy, SideKind, SidePath, bridge_cross_curvatures, setback_sides, unwrap_onto,
};
use super::{
    CornerJet, CrossFieldSamples, Kp20PaperBorderExtractReport, Kp20SetbackArc, OrientedSegment,
    chi_march_calibration, corner_jets, corner_pins, elevate_cubic_row_to_sextic,
    kp20_paper_border_corners, ring_fit, rows_from_cross_fields, section, section_fit,
    split_cubic_row,
};
use crate::fillet_surface::FilletSurfaceResult;

/// Corner-setback border: the vertex-blend boundary as n fillet-end
/// arcs set back to the TRIM-INTERSECTION corners — for each pair of
/// angularly adjacent fillets, the point where their contact curves
/// on the shared support face intersect. The boundary then lies
/// entirely on the fillet surfaces: each sector's two pieces are the
/// halves of ONE arc on ONE G² sweep (C² at the spoke point m by
/// construction), and the sector corners c sit where two sweeps meet
/// tangent to the same face (curvature-compatible corner data — the
/// property the G² sweep buys). This is the classical setback vertex
/// blend, shaped exactly like the characteristic layout.
pub fn kp20_paper_border_from_fillet_setback(
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
) -> Result<(Kp20PaperGuideBorder, Kp20PaperBorderExtractReport), String> {
    kp20_paper_border_from_fillet_setback_with(
        fillet_surfaces,
        vertex_position,
        SetbackJunctionPolicy::CornersOnly,
    )
}

/// [`kp20_paper_border_from_fillet_setback`] under a junction policy.
/// With [`SetbackJunctionPolicy::BridgeReflexFaces`], each reflex face
/// corner of a saddle takes a bridge side on its support (PK-527, see
/// `junctions`). The ring then has one side per strip plus one per
/// bridge, and a bridge side's arc carries its `bridge`.
pub fn kp20_paper_border_from_fillet_setback_with(
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
    policy: SetbackJunctionPolicy,
) -> Result<(Kp20PaperGuideBorder, Kp20PaperBorderExtractReport), String> {
    use popkern_core::geometry::surface::SurfaceEval;
    use popkern_core::geometry::uv_projection::project_point_to_uv;
    let n = fillet_surfaces.len();
    if n < 3 {
        return Err(format!("setback border needs >= 3 fillets, got {n}"));
    }

    // Angular order of the fillets about the vertex (best-fit frame
    // from the trim-curve midpoints).
    let rep_point = |fi: usize| -> Point3 {
        let c = &fillet_surfaces[fi].trim_curve_1;
        let (t0, t1) = c.parameter_range();
        c.point_at(0.5 * (t0 + t1))
    };
    let centroid = {
        let mut acc = [0.0_f64; 3];
        for fi in 0..n {
            let p = rep_point(fi);
            acc[0] += p.x;
            acc[1] += p.y;
            acc[2] += p.z;
        }
        Point3::new(acc[0] / n as f64, acc[1] / n as f64, acc[2] / n as f64)
    };
    let mut axis = Vec3::new(0.0, 0.0, 0.0);
    for fi in 0..n {
        let a = rep_point(fi) - centroid;
        let b = rep_point((fi + 1) % n) - centroid;
        axis = axis + a.cross(&b);
    }
    let axis = axis
        .normalized()
        .ok_or_else(|| "degenerate fillet arrangement".to_string())?;
    let seed = rep_point(0) - centroid;
    let u0 = (seed - axis * seed.dot(&axis))
        .normalized()
        .ok_or_else(|| "degenerate angular frame".to_string())?;
    let v0 = axis.cross(&u0);
    let angle_of = |p: Point3| -> f64 {
        let d = p - centroid;
        d.dot(&v0).atan2(d.dot(&u0))
    };
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| {
        angle_of(rep_point(a))
            .partial_cmp(&angle_of(rep_point(b)))
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // Setback corners between angularly adjacent fillets — see
    // `kp20_paper_border_corners` for the pair-selection and
    // truncated-trim recovery model.
    let corners =
        kp20_paper_border_corners::setback_corners(fillet_surfaces, &order, vertex_position);
    let max_corner_gap = corners.iter().map(|c| c.0).fold(0.0_f64, f64::max);
    // The boundary's sides: one per strip, plus a bridge per reflex
    // face corner under the bridging policy (PK-527).
    let sides = setback_sides(fillet_surfaces, &order, &corners, vertex_position, policy)?;
    let n = sides.len();

    // Corners cut exactly at the trim intersection. (A fractional
    // pullback along both arcs onto the shared tangent plane — the
    // classical setback-distance choice — was probed and retired:
    // production keeps the exact trim-intersection corner.)
    // One arc per fillet between its two corners. Pass 1: the arc's
    // section path on the strip (PK-401, see `section_path`): both
    // ends are cross-section cuts, the axial displacement between the
    // two corner stations (non-zero on mixed radii) is taken up in
    // the arc's interior, so every cap corner turns by the face's
    // apex angle.
    let section_trace = popkern_core::trace::trace_enabled("setback");
    let mut section_axis_fallbacks = 0_usize;
    let mut paths: Vec<SidePath> = Vec::with_capacity(n);
    for (k, side) in sides.iter().enumerate() {
        let SideKind::Strip { fillet_index } = side.kind else {
            if let SideKind::Bridge(bridge) = &side.kind {
                paths.push(SidePath::Bridge(bridge.clone()));
            }
            continue;
        };
        let fillet = &fillet_surfaces[fillet_index];
        let uv_prev = project_point_to_uv(&fillet.surface, &side.start)
            .ok_or_else(|| "corner projection failed".to_string())?;
        let uv_next = project_point_to_uv(&fillet.surface, &side.end)
            .ok_or_else(|| "corner projection failed".to_string())?;
        // A strip spans under a half turn, so its arc is the short way
        // round on a periodic chart: unwrap the end onto the start's
        // branch (the raw projections straddled the cylinder seam on
        // the saddle's convex strips, and the path ran 270° through
        // material, PK-527).
        let uv_next = unwrap_onto(&fillet.surface, uv_prev, uv_next);
        let path = SectionPath::resolve(
            fillet,
            fillet_index,
            (uv_prev, side.start),
            (uv_next, side.end),
        );
        if path.axial.is_none() {
            section_axis_fallbacks += 1;
        }
        if section_trace {
            eprintln!(
                "  [section-path] arc {k}: axial {:?}, displacement {:?}",
                path.axial,
                path.axial_displacement()
            );
        }
        paths.push(SidePath::Strip(path));
    }
    let eval_span = |k: usize, t: f64| -> Point3 {
        let (u, v) = paths[k].uv(t);
        paths[k].surface(fillet_surfaces).point_at(u, v)
    };
    // Pass 2: corners = midpoints of adjacent cut ends (they flank
    // the shared face symmetrically, so the midpoint sits essentially
    // on it).
    let new_corners: Vec<Point3> = (0..n)
        .map(|j| {
            let a_end = eval_span(j, 1.0);
            let b_start = eval_span((j + 1) % n, 0.0);
            Point3::new(
                0.5 * (a_end.x + b_start.x),
                0.5 * (a_end.y + b_start.y),
                0.5 * (a_end.z + b_start.z),
            )
        })
        .collect();
    // Pass 3: sector pieces interpolate arc samples and END exactly
    // at the shared corners.
    let mut sectors = Vec::with_capacity(n);
    let mut first_order_marches = 0_usize;
    let calibration = chi_march_calibration(n)
        .ok_or_else(|| "characteristic march calibration failed".to_string())?;
    let interpolate = |pts: [Point3; 4]| -> [Point3; 4] {
        let inner = |a: f64, b: f64, c: f64, d: f64| {
            Point3::new(
                (a * pts[0].x + b * pts[1].x + c * pts[2].x + d * pts[3].x) / 6.0,
                (a * pts[0].y + b * pts[1].y + c * pts[2].y + d * pts[3].y) / 6.0,
                (a * pts[0].z + b * pts[1].z + c * pts[2].z + d * pts[3].z) / 6.0,
            )
        };
        [
            pts[0],
            inner(-5.0, 18.0, -9.0, 2.0),
            inner(2.0, -9.0, 18.0, -5.0),
            pts[3],
        ]
    };
    let mut setback_arcs: Vec<Kp20SetbackArc> = Vec::with_capacity(n);
    let mut section_fit_residual = 0.0_f64;
    let mut section_fit_fallbacks = 0_usize;
    for k in 0..n {
        let c_start = new_corners[(k + n - 1) % n];
        let c_end = new_corners[k];
        // The coarse arc: ONE cubic per arc through the shared corners
        // (crossing classification and topology splitting consume it).
        let full_arc = interpolate([
            c_start,
            eval_span(k, 1.0 / 3.0),
            eval_span(k, 2.0 / 3.0),
            c_end,
        ]);
        // The rim rows: the strip's section fitted as one C²-at-m
        // two-piece degree-6 spline through the shared corners
        // (PK-398, see `section`). The cubic interpolant's sagitta
        // (measured 7.775e-4 at r = 0.1) was the served cap↔strip
        // tangent-jump class; the fit is accepted only when it beats
        // the interpolant, else the arc serves the elevated cubic.
        let sample = |t: f64| eval_span(k, t);
        let cubic_residual = section_fit::cubic_section_residual(&full_arc, &sample);
        let (first_half, second_half) = match section_fit::fit_section_spline(
            &sample, c_start, c_end,
        ) {
            // `<=`: at a gapped corner both rows interpolate the same
            // off-surface corner point, so their worst residuals
            // coincide there (1.6e-3 on the blunt pentagon G2 sweeps)
            // while the spline is exact everywhere else; a strict
            // comparison served the cubic instead (measured rim
            // tangency 1.8e-1).
            Some(fit) if fit.residual <= cubic_residual => {
                section_fit_residual = section_fit_residual.max(fit.residual);
                if section_trace {
                    eprintln!(
                        "  [section] arc {k}: cubic sagitta {cubic_residual:.3e}, spline residual {:.3e}",
                        fit.residual
                    );
                }
                (fit.halves[0], fit.halves[1])
            }
            fit => {
                section_fit_fallbacks += 1;
                section_fit_residual = section_fit_residual.max(cubic_residual);
                if section_trace {
                    eprintln!(
                        "  [section] arc {k}: fallback to the cubic (sagitta {cubic_residual:.3e}; fit {:?})",
                        fit.map(|f| f.residual)
                    );
                }
                let (a, b) = split_cubic_row(&full_arc);
                (
                    elevate_cubic_row_to_sextic(a),
                    elevate_cubic_row_to_sextic(b),
                )
            }
        };
        // A bridge side names the strip it leaves (its predecessor).
        let fillet_index = paths[k]
            .fillet_index()
            .or_else(|| paths[(k + n - 1) % n].fillet_index())
            .ok_or_else(|| "a bridge side follows another bridge".to_string())?;
        setback_arcs.push(Kp20SetbackArc {
            fillet_index,
            bridge: match &paths[k] {
                SidePath::Bridge(bridge) => Some((**bridge).clone()),
                SidePath::Strip(_) => None,
            },
            cubic: full_arc,
            section_halves: [first_half, second_half],
            corner_prev: c_start,
            corner_next: c_end,
        });
    }
    // Pass 4: cross rows as fitted fields with SHARED junction data
    // (PK-398, PK-400). Each half's corner end takes the corner-
    // compatible jet of the NEIGHBOURING arc's rim (`rim_corner_jet`):
    // arc k's first half starts at c_{k-1} where arc k-1's second half
    // ends; its second half ends at c_k where arc k+1's first half
    // starts. At the spoke point m_k both halves march with ONE vector
    // (the mean of the two χ̃ calibrations there — the first half is
    // the sector's near_m piece, calibration[1]; the second half the
    // near_c piece, calibration[0]; the halves used to take each
    // other's calibration and reached m with χ̃'s skewed CORNER angle).
    // The sampled fields are then fitted ring-wide with the spoke
    // junctions and corner blocks as shared unknowns (`ring_fit`), so
    // the ring is C² across every m and both pieces at every c carry
    // one corner block by construction.
    let shims: Vec<OrientedSegment> = (0..n)
        .map(|k| match &paths[k] {
            SidePath::Strip(path) => {
                let fillet = &fillet_surfaces[path.fillet_index];
                OrientedSegment {
                    curve: fillet.trim_curve_1.clone(),
                    surface: fillet.surface.clone(),
                    is_fillet_end: true,
                    reversed: false,
                }
            }
            SidePath::Bridge(bridge) => OrientedSegment {
                curve: bridge.curve.clone(),
                surface: bridge.plane.clone(),
                is_fillet_end: false,
                reversed: false,
            },
        })
        .collect();
    // Corner classification (see `corner_jets`): per corner c_k the
    // shared corner jets (`.0` serves arc k's second half ending at
    // c_k, `.1` arc k+1's first half), or `None` at a gapped or reflex
    // corner (χ̃-calibrated ends, no shared block), plus the cap
    // boundary's turn at every corner.
    let side_surfaces: Vec<&popkern_core::geometry::SurfaceGeometry> =
        paths.iter().map(|p| p.surface(fillet_surfaces)).collect();
    let classified = corner_jets::classify_corners(
        &setback_arcs,
        &side_surfaces,
        &shims,
        &calibration,
        vertex_position,
        section_trace,
    );
    let corner_jets::CornerClassification {
        jets: corner_jets,
        turns: corner_turns,
        gapped: gapped_corners,
        reflex: reflex_corners,
    } = classified;
    let corner_jet_at = |k: usize| -> (Option<CornerJet>, Option<CornerJet>) {
        (
            corner_jets[(k + n - 1) % n].map(|j| j.1),
            corner_jets[k].map(|j| j.0),
        )
    };
    let mut samples: Vec<[CrossFieldSamples; 2]> = Vec::with_capacity(n);
    for k in 0..n {
        let [first_half, second_half] = setback_arcs[k].section_halves;
        let (corner_prev, corner_next) = corner_jet_at(k);
        let m_vector =
            (shims[k].calibrated_march(&first_half, 6, &calibration[1].end, vertex_position)
                + shims[k].calibrated_march(
                    &second_half,
                    0,
                    &calibration[0].start,
                    vertex_position,
                ))
                * 0.5;
        let start = match corner_prev {
            Some(corner) => -corner.d1,
            None => {
                shims[k].calibrated_march(&first_half, 0, &calibration[1].start, vertex_position)
            }
        };
        let end = match corner_next {
            Some(corner) => -corner.d1,
            None => {
                shims[k].calibrated_march(&second_half, 6, &calibration[0].end, vertex_position)
            }
        };
        let s0 =
            shims[k].march_cross_field(&first_half, [start, m_vector], &mut first_order_marches);
        let s1 =
            shims[k].march_cross_field(&second_half, [m_vector, end], &mut first_order_marches);
        samples.push([s0, s1]);
    }
    // Bridges take the strips' section curvature at their feet (PK-527).
    bridge_cross_curvatures(&paths, &setback_arcs, &shims, &mut samples);
    // Ring-wide fit of every cross row with the shared junction data
    // as unknowns (PK-400, see `ring_fit`): the spoke junctions and
    // corner blocks are single unknowns, so the C²-at-m and shared-
    // block invariants hold by construction, and the rows follow the
    // sampled in-plane fields instead of pins derived from them.
    // Tangency stations per half (see `ring_fit::tangency_stations`).
    let tangency_stations = |k: usize, half: usize, row0: &[Point3; 7], s: &CrossFieldSamples| {
        ring_fit::tangency_stations(
            paths[k].surface(fillet_surfaces),
            &|span_t| paths[k].uv(span_t),
            half,
            row0,
            s,
        )
    };
    let ring: Vec<[ring_fit::RingPiece<'_>; 2]> = (0..n)
        .map(|k| {
            let [first_half, second_half] = &setback_arcs[k].section_halves;
            [
                ring_fit::RingPiece {
                    row0: first_half,
                    samples: &samples[k][0],
                    tangency: tangency_stations(k, 0, first_half, &samples[k][0]),
                },
                ring_fit::RingPiece {
                    row0: second_half,
                    samples: &samples[k][1],
                    tangency: tangency_stations(k, 1, second_half, &samples[k][1]),
                },
            ]
        })
        .collect();
    let mut ring_fit_fallback = false;
    let fitted: Vec<ring_fit::RingRows> = match ring_fit::fit_ring_cross_fields(&ring, &corner_jets)
    {
        Some(rows) => rows,
        None => {
            // Singular normal equations (degenerate sampling): each
            // piece keeps its own pinned fit with the corner jets
            // only; the assembly's overlap average and the fan's rim
            // reconcile remain the containment for the shared data.
            ring_fit_fallback = true;
            if section_trace {
                eprintln!("  [ring-fit] singular normal equations — per-piece pinned fits");
            }
            (0..n)
                .map(|k| {
                    let (corner_prev, corner_next) = corner_jet_at(k);
                    let (p0, p0c) = corner_prev
                        .map(|c| corner_pins(c, false, None))
                        .unwrap_or(([None; 7], [None; 7]));
                    let (p1, p1c) = corner_next
                        .map(|c| corner_pins(c, true, None))
                        .unwrap_or(([None; 7], [None; 7]));
                    let [s0, s1] = &samples[k];
                    [
                        (
                            section::fit_sextic_field_pinned(&s0.jets, p0),
                            section::fit_sextic_field_pinned(&s0.curvatures, p0c),
                        ),
                        (
                            section::fit_sextic_field_pinned(&s1.jets, p1),
                            section::fit_sextic_field_pinned(&s1.curvatures, p1c),
                        ),
                    ]
                })
                .collect()
        }
    };
    drop(ring);
    for k in 0..n {
        let [first_half, second_half] = setback_arcs[k].section_halves;
        let [(f0, f0c), (f1, f1c)] = &fitted[k];
        let piece_a = rows_from_cross_fields(&first_half, f0, f0c);
        let piece_b = rows_from_cross_fields(&second_half, f1, f1c);
        sectors.push((piece_a, piece_b));
    }
    type PiecePair = ([[Point3; 7]; 3], [[Point3; 7]; 3]);
    let pieces: Vec<PiecePair> = sectors;
    let mut grouped = Vec::with_capacity(n);
    for s in 0..n {
        let piece_a = pieces[s].1; // [m_s → c_s]
        let piece_b = pieces[(s + 1) % n].0; // [c_s → m_{s+1}]
        grouped.push([
            Kp20PaperGuideBorderPiece { rows: piece_a },
            Kp20PaperGuideBorderPiece { rows: piece_b },
        ]);
    }

    // The former corner face-tangency clamp (rows 1..2 columns nearest
    // c projected onto the outward side of the shared face plane) is
    // retired with PK-398: it was containment for extraction noise in
    // rows marched from control points. The corner data is now the
    // neighbouring rim's own Hermite block, which on rolling-ball
    // (G1) sections legitimately curves below the face plane next to
    // c; clamping it tilted the jets out of their true tangent planes
    // (measured rim tangency 7e-2..1.9e-1 on the blunt pentagon G1
    // control against 5e-2).
    let _ = &order;
    let corner_block_mismatch =
        kp20_paper_border_corners::corner_block_mismatches(&grouped, section_trace);
    Ok((
        Kp20PaperGuideBorder { sectors: grouped },
        Kp20PaperBorderExtractReport {
            n,
            depth: 0.0,
            max_loop_gap: max_corner_gap,
            first_order_marches,
            reversed_for_orientation: false,
            trim_corners: corners.iter().map(|c| c.1).collect(),
            setback_arcs,
            recovered_face_surfaces: 0,
            section_fit_residual,
            section_fit_fallbacks,
            gapped_corners,
            reflex_corners,
            corner_turns,
            section_axis_fallbacks,
            corner_block_mismatch,
            ring_fit_fallback,
        },
    ))
}
