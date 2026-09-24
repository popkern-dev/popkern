//! Trace-gated fan diagnostics for the setback cap commit: chart
//! landmark positions of sector 0's pieces and the per-sector rim-leg
//! probe measuring (a) distance to the committed topology ring, (b)
//! distance to the nearest fillet strip, (c) distance to the INPUT
//! border rim — discriminating "conditioner slid the rim along the
//! strip" (benign) from "conditioner left the strip" (a G0 wall).
//! Extracted verbatim from `build_cap_face_fan` (structural line
//! budget); active only under `POPKERN_TRACE=setback`.

use super::*;

/// Print the `[chart]`, `[rim]`, and `[ring-edge]` probes for the
/// assembled fan. No-op unless the `setback` trace channel is on.
pub(super) fn trace_fan_chart_and_rim(
    brep: &BRep,
    blend_face: FaceHandle,
    sectors: &[popkern_core::geometry::gregory::Kp20MacroPatch],
    border: &popkern_core::geometry::gregory::Kp20PaperGuideBorder,
    fillet_surfaces: &[FilletSurfaceResult],
) {
    if popkern_core::trace::trace_enabled("setback") {
        let dist_to_loop = |pt: Point3| -> f64 {
            let mut best = f64::INFINITY;
            if let Ok(hes) = face_loop_half_edges(brep, blend_face) {
                for he in hes {
                    let edge = brep.half_edges[he].edge;
                    if let Some(c) = brep.edges[edge].curve.as_ref() {
                        best = best.min(curve_distance(&c.geometry, c.t_start, c.t_end, pt));
                    }
                }
            }
            best
        };
        for (s, sec) in sectors.iter().enumerate().take(1) {
            use popkern_core::geometry::gregory::Kp20PieceId as Pid;
            for (label, pt) in [
                ("P1(0,0) o?", sec.piece(Pid::NearOriginLeft).cps[0][0]),
                ("P1(0,6)", sec.piece(Pid::NearOriginLeft).cps[0][6]),
                ("P1(6,0)", sec.piece(Pid::NearOriginLeft).cps[6][0]),
                ("P1(6,6)", sec.piece(Pid::NearOriginLeft).cps[6][6]),
                ("P4(0,6)", sec.piece(Pid::NearMidpointLeft).cps[0][6]),
                ("P4(6,6)", sec.piece(Pid::NearMidpointLeft).cps[6][6]),
                ("P4(6,0)", sec.piece(Pid::NearMidpointLeft).cps[6][0]),
                ("P2(6,6)", sec.piece(Pid::NearOriginRight).cps[6][6]),
                ("P2(6,0)", sec.piece(Pid::NearOriginRight).cps[6][0]),
                ("P2(0,6)", sec.piece(Pid::NearOriginRight).cps[0][6]),
                ("P3(0,0)", sec.piece(Pid::NearMidpointRight).cps[0][0]),
                ("P3(6,6)", sec.piece(Pid::NearMidpointRight).cps[6][6]),
                ("P3(0,6)", sec.piece(Pid::NearMidpointRight).cps[0][6]),
                ("P3(6,0)", sec.piece(Pid::NearMidpointRight).cps[6][0]),
            ] {
                eprintln!(
                    "  [chart] s{s} {label}: ({:+.3},{:+.3},{:+.3}) loop-dist {:.2e}",
                    pt.x,
                    pt.y,
                    pt.z,
                    dist_to_loop(pt)
                );
            }
        }
        // [rim] probe: per sector, sample the four landed rim legs and
        // measure (a) distance to the committed topology ring, (b)
        // distance to the nearest fillet strip surface, (c) distance
        // to the INPUT border rim (the arc cubics). Discriminates
        // "conditioner slid the rim along the strip" (benign — ring
        // can follow) from "conditioner left the strip" (G0 wall).
        {
            use popkern_core::geometry::gregory::Kp20PieceId as Pid;
            let bez7 = |cps: &[Point3; 7], t: f64| -> Point3 {
                let mut p = *cps;
                for level in (1..7).rev() {
                    for i in 0..level {
                        p[i] = Point3::new(
                            p[i].x + (p[i + 1].x - p[i].x) * t,
                            p[i].y + (p[i + 1].y - p[i].y) * t,
                            p[i].z + (p[i + 1].z - p[i].z) * t,
                        );
                    }
                }
                p[0]
            };
            let strip_dist = |p: Point3| -> f64 {
                let mut best = f64::INFINITY;
                for fillet in fillet_surfaces {
                    if let Some((u, v)) = popkern_core::geometry::uv_projection::project_point_to_uv(
                        &fillet.surface,
                        &p,
                    ) {
                        best = best.min(fillet.surface.point_at(u, v).distance_to(&p));
                    }
                }
                best
            };
            let input_rim_dist = |p: Point3| -> f64 {
                let mut best = f64::INFINITY;
                for sector in &border.sectors {
                    for piece in sector {
                        for i in 0..=128 {
                            let t = i as f64 / 128.0;
                            best = best.min(bez7(&piece.rows[0], t).distance_to(&p));
                        }
                    }
                }
                best
            };
            for (s, sec) in sectors.iter().enumerate() {
                for (label, cps) in [
                    ("A1 P4j6", piece_col(sec.piece(Pid::NearMidpointLeft), 6)),
                    ("A2 P3j6", piece_col(sec.piece(Pid::NearMidpointRight), 6)),
                    ("B1 P3i6", piece_row(sec.piece(Pid::NearMidpointRight), 6)),
                    ("B2 P2i6", piece_row(sec.piece(Pid::NearOriginRight), 6)),
                ] {
                    let mut worst = (0.0_f64, 0.0_f64);
                    for i in 0..=8 {
                        let t = i as f64 / 8.0;
                        let d = dist_to_loop(bez7(&cps, t));
                        if d > worst.0 {
                            worst = (d, t);
                        }
                    }
                    let pw = bez7(&cps, worst.1);
                    eprintln!(
                        "  [rim] s{s} {label}: worst-loop {:.2e} @t={:.2} strip {:.2e} input-rim {:.2e}",
                        worst.0,
                        worst.1,
                        strip_dist(pw),
                        input_rim_dist(pw),
                    );
                    // For badly-off rim points, dump the per-edge ring
                    // distances to localize the missing image span.
                    if worst.0 > 5.0e-3
                        && let Ok(hes) = face_loop_half_edges(brep, blend_face)
                    {
                        for he in hes {
                            let edge = brep.half_edges[he].edge;
                            let v0 = brep.half_edges[he].origin;
                            let p0 = brep.vertices[v0].position;
                            let d = brep.edges[edge].curve.as_ref().map_or(f64::INFINITY, |c| {
                                curve_distance(&c.geometry, c.t_start, c.t_end, pw)
                            });
                            eprintln!(
                                "    [ring-edge] {edge:?} from ({:+.3},{:+.3},{:+.3}) dist {d:.2e}",
                                p0.x, p0.y, p0.z
                            );
                        }
                    }
                }
            }
        }
    }
}
