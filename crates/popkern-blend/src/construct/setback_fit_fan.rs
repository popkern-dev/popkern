//! The variational-fit cap on a carved setback region (PK-574).
//!
//! The carved blend face is bounded by the n setback sides. The core fit
//! (`popkern_core::geometry::variational_cap`) returns n NURBS sectors
//! meeting at a centre o; this module subdivides the blend face into
//! those sectors (spokes from o to each side's midpoint m_s), gives every
//! cap edge its exact iso-curve and every cap half-edge its exact Linear
//! iso pcurve, re-authors the neighbours' pcurves, and certifies the
//! cap's seams (see `setback_fit_certify`). Any failure is a cap-stage
//! decline; the surgery restores the snapshot.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::pcurve::Pcurve;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::variational_cap::{
    CapEnergy, VariationalCap, VariationalCapOptions, VariationalCapTargets, fit_variational_cap,
};
use popkern_core::geometry::{
    CurveGeometry, CurveOnEdge, NurbsSurface, SurfaceGeometry, SurfaceOnFace,
};

use super::*;
use crate::blend_continuity_evidence::ContinuityGates;
use crate::construct::kp20_paper_border_extract::Kp20SetbackArc;
use crate::construct::setback_decline::CapCommitError;
use crate::construct::setback_fit_cap::fit_sides;
use crate::types::{SetbackFitCertified, SetbackFitMode};

#[path = "setback_fit_certify.rs"]
mod setback_fit_certify;

/// A certified fit serve.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct FitServe {
    pub mode: SetbackFitMode,
    pub certified: SetbackFitCertified,
}

/// Fit, build and certify the cap over `blend_face`.
pub(super) fn commit_fit_fan(
    brep: &mut BRep,
    blend_face: FaceHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    arcs: &[Kp20SetbackArc],
    mode: SetbackFitMode,
    tolerance: &Tolerance,
) -> Result<(Vec<FaceHandle>, FitServe), CapCommitError> {
    let gates = ContinuityGates::commercial();
    // Neighbour orientation from the carved topology: the face across
    // each side's middle from the blend face.
    let orientation: Vec<bool> = arcs
        .iter()
        .map(|arc| side_neighbour_orientation(brep, blend_face, arc))
        .collect::<Result<_, String>>()?;
    let outward = |k: usize, _: Point3, _: popkern_core::geometry::Vec3| orientation[k];
    let sides = fit_sides(arcs, fillet_surfaces, &outward)?;
    let (targets, options) = match mode {
        SetbackFitMode::CurvatureWeighted => (
            VariationalCapTargets::g2(tolerance.linear(), 1e-4, gates.g2_curvature),
            VariationalCapOptions {
                max_refinement_levels: 0,
                energy: CapEnergy::MeanCurvatureVariation,
                ..Default::default()
            },
        ),
        SetbackFitMode::Tangent => (
            VariationalCapTargets::g1(tolerance.linear(), 1e-6),
            VariationalCapOptions {
                max_refinement_levels: 2,
                energy: CapEnergy::Bending,
                ..Default::default()
            },
        ),
    };
    let mut cap = fit_variational_cap(&sides, &targets, &options, tolerance.linear())
        .map_err(|e| format!("variational fit: {e}"))?;
    if crate::dispatch_overrides::fit_cap_rim_tilted() {
        tilt_rim_row(&mut cap.sectors[0]);
    }
    if mode == SetbackFitMode::CurvatureWeighted {
        curvature_witness_gate(&cap, gates.g2_curvature)?;
    }
    let cap_faces = build_sector_faces(brep, blend_face, &cap, &sides)?;
    author_cap_pcurves(brep, &cap_faces, &cap)?;
    super::author_setback_pcurves(brep, &[], &cap_faces)?;
    author_neighbour_pcurves(brep, &cap_faces, tolerance)?;
    super::setback_boundary_tolerance::record_boundary_tolerances(brep, &cap_faces);
    let certified = setback_fit_certify::certify_cap(brep, &cap_faces, &gates, &cap.report)?;
    let required = match mode {
        SetbackFitMode::CurvatureWeighted => crate::types::ContinuityLevel::G2,
        SetbackFitMode::Tangent => crate::types::ContinuityLevel::G1,
    };
    if level_rank(certified.level) < level_rank(required) {
        return Err(CapCommitError {
            reason: format!(
                "variational fit ({mode:?}) certified {:?}, needs {required:?}: rim (g0 {:.3e}, g1 {:.3e}, g2 {:.3e}), spoke (g0 {:.3e}, g1 {:.3e}, g2 {:.3e})",
                certified.level,
                certified.rim[0],
                certified.rim[1],
                certified.rim[2],
                certified.spoke[0],
                certified.spoke[1],
                certified.spoke[2],
            ),
            fit_certified: Some(certified),
        });
    }
    if let Err(intersection) = crate::apex_construction::exact_boundary_cap::
        cap_quality_postcondition::validate_candidate_region_self_intersection(brep, &cap_faces)
    {
        return Err(CapCommitError {
            reason: format!("variational fit cap self-intersection postcondition failed: {intersection}"),
            fit_certified: Some(certified),
        });
    }
    Ok((cap_faces, FitServe { mode, certified }))
}

/// Co-parameterise every neighbour half-edge along the rim with the
/// rim edge's new curve (the fit's iso-curve). The rim geometry was
/// replaced, so a neighbour's split or projected pcurve is stale or
/// missing — analytic neighbours (cylinder, torus strips) are not
/// re-authored by the carve's pcurve pass at all. The certifier reads
/// the seam under fraction correspondence, which needs exactly this
/// co-parameterisation. The verification gate is ten linear
/// tolerances: the rim curve lies off the neighbour by the fit's G0
/// residual, which the certifier then bounds.
fn author_neighbour_pcurves(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    tolerance: &Tolerance,
) -> Result<(), String> {
    let gate = 10.0 * tolerance.linear();
    for &face in cap_faces {
        for he in face_loop_half_edges(brep, face)? {
            let twin = brep.half_edges[he].twin;
            let neighbour = brep.loops[brep.half_edges[twin].loop_].face;
            if cap_faces.contains(&neighbour) {
                continue;
            }
            let edge = brep.half_edges[he].edge;
            let curve = brep.edges[edge]
                .curve
                .clone()
                .ok_or("rim edge has no curve")?;
            let authored = popkern_core::topology::edge_pcurve_reauthor::coparameterised_pcurves(
                brep, edge, &curve, gate, tolerance,
            )
            .ok_or_else(|| {
                format!("rim edge {edge:?}: neighbour pcurve does not verify within {gate:.1e}")
            })?;
            for (h, pcurve) in authored {
                if h == twin {
                    brep.half_edges[twin].pcurve = Some(pcurve);
                }
            }
        }
    }
    Ok(())
}

/// The negative-control corruption behind
/// [`crate::dispatch_overrides::with_fit_cap_rim_tilted`]: lift the
/// control row next to the `u = 1` rim along the unit surface normal
/// by a fifth of the row spacing. Only interior columns move — columns
/// 0 and 1 carry the spoke `v = 0` and its cross tangent, the last two
/// the `v = 1` rim and its — so exactly one thing breaks: the `u = 1`
/// rim's cross tangent (position intact).
fn tilt_rim_row(sector: &mut NurbsSurface) {
    let (m_u, m_v) = (sector.num_u, sector.num_v);
    if m_u < 2 || m_v < 5 {
        return;
    }
    for j in 2..m_v - 2 {
        let (rim, inner) = ((m_u - 1) * m_v + j, (m_u - 2) * m_v + j);
        let spacing = sector.control_points[rim].distance_to(&sector.control_points[inner]);
        let v = j as f64 / (m_v - 1) as f64;
        let Some(normal) = sector.normal_at(1.0, v).normalized() else {
            continue;
        };
        sector.control_points[inner] = sector.control_points[inner] + normal * (0.2 * spacing);
    }
}

/// Decline a curvature-weighted fit before any topology when its own
/// rim samples already prove G2 cannot certify. A sample is an attained
/// value, so it bounds the seam's supremum from below; the fit reports
/// the spectral norm of the second-fundamental-form jump, which is at
/// most twice the largest component the certifier bounds (2×2
/// symmetric), hence the factor 2. Measured (E4): certification was
/// 70 % of the rung's time, mostly spent on these doomed caps.
fn curvature_witness_gate(cap: &VariationalCap, g2_gate: f64) -> Result<(), CapCommitError> {
    match cap.report.rim.iter().find(|r| r.g2 > 2.0 * g2_gate) {
        Some(r) => Err(CapCommitError::from(format!(
            "curvature-weighted fit: sampled rim curvature jump {:.3e} on side {} ({:?} half, t = {:.3}) exceeds twice the G2 gate {g2_gate:.1e}; G2 cannot certify",
            r.g2, r.side, r.half, r.parameter
        ))),
        None => Ok(()),
    }
}

fn level_rank(level: crate::types::ContinuityLevel) -> u8 {
    match level {
        crate::types::ContinuityLevel::G0 => 0,
        crate::types::ContinuityLevel::G1 => 1,
        crate::types::ContinuityLevel::G2 => 2,
    }
}

/// `SurfaceOnFace::orientation` of the face across the blend face from
/// the middle of `arc` (true: the neighbour's natural normal is its
/// outward normal).
fn side_neighbour_orientation(
    brep: &BRep,
    blend_face: FaceHandle,
    arc: &Kp20SetbackArc,
) -> Result<bool, String> {
    let curve = match &arc.bridge {
        Some(bridge) => bridge.curve.clone(),
        None => arc.section_curve(),
    };
    let (a, b) = curve.parameter_range();
    let mid = curve.point_at(0.5 * (a + b));
    let mut best: Option<(f64, HalfEdgeHandle)> = None;
    for he in face_loop_half_edges(brep, blend_face)? {
        let edge = brep.half_edges[he].edge;
        let Some(c) = brep.edges[edge].curve.as_ref() else {
            continue;
        };
        let d = curve_distance(&c.geometry, c.t_start, c.t_end, mid);
        if best.is_none_or(|(bd, _)| d < bd) {
            best = Some((d, he));
        }
    }
    let (_, he) = best.ok_or("blend face has no curved boundary edge")?;
    let twin = brep.half_edges[he].twin;
    let face = brep.loops[brep.half_edges[twin].loop_].face;
    brep.faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| s.orientation)
        .ok_or_else(|| "side neighbour face has no surface".to_string())
}

/// Split the blend face into the fit's sectors; returns the sector
/// faces in sector order.
fn build_sector_faces(
    brep: &mut BRep,
    blend_face: FaceHandle,
    cap: &VariationalCap,
    sides: &[popkern_core::geometry::variational_cap::VariationalCapSide],
) -> Result<Vec<FaceHandle>, String> {
    let n = cap.sectors.len();
    // 1. Spoke vertices on the sides, moved exactly onto the fit's m_s.
    let mut m_vertices = Vec::with_capacity(n);
    for &m in &cap.spoke_points {
        let v = split_blend_boundary_at(brep, blend_face, m)?;
        brep.vertices[v].position = m;
        m_vertices.push(v);
    }
    // 2. Spokes: o by mev from m_0, the rest closed onto o by mef.
    let o = cap.centre;
    let mut o_vertex: Option<VertexHandle> = None;
    for (s, &m_v) in m_vertices.iter().enumerate() {
        let spoke = spoke_curve(&cap.sectors[s]);
        match o_vertex {
            None => {
                let he = half_edge_at_vertex_on_face_containing(brep, m_v, &o)?;
                let mev = brep.mev(he, o).map_err(|e| format!("mev spoke 0: {e}"))?;
                brep.edges[mev.edge].curve = Some(spoke);
                o_vertex = Some(mev.vertex);
            }
            Some(o_v) => {
                let (he_m, he_o) = half_edge_pair_in_common_loop(brep, m_v, o_v)?;
                let mef = brep
                    .mef(he_m, he_o)
                    .map_err(|e| format!("mef spoke {s}: {e}"))?;
                brep.edges[mef.edge].curve = Some(spoke);
            }
        }
    }
    let o_vertex = o_vertex.ok_or("no spoke vertex")?;
    brep.vertices[o_vertex].position = o;
    // 3. Sector faces: the faces around o, matched to sectors by their
    //    corners (o, m_s, c_s, m_{s+1}).
    let around: Vec<FaceHandle> = {
        let mut faces: Vec<FaceHandle> = brep
            .half_edges
            .values()
            .filter(|h| h.origin == o_vertex)
            .map(|h| brep.loops[h.loop_].face)
            .collect();
        faces.sort();
        faces.dedup();
        faces
    };
    if around.len() != n {
        return Err(format!(
            "{} faces around the centre, expected {n}",
            around.len()
        ));
    }
    let mut assigned: Vec<Option<FaceHandle>> = vec![None; n];
    for &face in &around {
        let verts: Vec<Point3> = face_loop_half_edges(brep, face)?
            .iter()
            .map(|&he| brep.vertices[brep.half_edges[he].origin].position)
            .collect();
        let score = |s: usize| {
            let sec = &cap.sectors[s];
            [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)]
                .iter()
                .map(|&(u, v)| {
                    let p = sec.point_at(u, v);
                    verts
                        .iter()
                        .map(|q| q.distance_to(&p))
                        .fold(f64::INFINITY, f64::min)
                })
                .fold(0.0, f64::max)
        };
        let (best, d) = (0..n)
            .map(|s| (s, score(s)))
            .min_by(|a, b| a.1.total_cmp(&b.1))
            .ok_or("no sector")?;
        if d > 1e-6 || assigned[best].is_some() {
            return Err(format!(
                "face {face:?} matches no free sector (best {best} at {d:.3e})"
            ));
        }
        assigned[best] = Some(face);
    }
    let faces: Vec<FaceHandle> = assigned
        .into_iter()
        .collect::<Option<_>>()
        .ok_or("unassigned sector")?;
    // 4. Surfaces, oriented by the neighbour's outward normal at m_s.
    for (s, &face) in faces.iter().enumerate() {
        let sec = &cap.sectors[s];
        let side = &sides[s];
        let p = sec.point_at(1.0, 0.5);
        let n_cap = sec.partial_u(1.0, 0.5).cross(&sec.partial_v(1.0, 0.5));
        let (u, v) =
            popkern_core::geometry::uv_projection::project_point_to_uv(&side.neighbour, &p)
                .ok_or("sector rim does not project onto its neighbour")?;
        let mut n_nb = side.neighbour.normal_at(u, v);
        if side.neighbour_reversed {
            n_nb = n_nb * -1.0;
        }
        brep.faces[face].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Nurbs(sec.clone()),
            orientation: n_cap.dot(&n_nb) > 0.0,
        });
    }
    // 5. Rim edges: each non-spoke cap edge takes its sector's rim
    //    iso-curve over its own parameter window; rim vertices snap onto
    //    it.
    for (s, &face) in faces.iter().enumerate() {
        let sec = &cap.sectors[s];
        let rims = [sec.iso_curve_at_u(1.0), sec.iso_curve_at_v(1.0)];
        for he in face_loop_half_edges(brep, face)? {
            let edge = brep.half_edges[he].edge;
            let (a, b) = (
                brep.half_edges[he].origin,
                brep.half_edges[brep.half_edges[he].twin].origin,
            );
            if a == o_vertex || b == o_vertex {
                continue;
            }
            let (pa, pb) = (brep.vertices[a].position, brep.vertices[b].position);
            // A rim piece lies on one iso-curve (its corner end sits on
            // both): take the rim that carries both ends.
            let (rim, ta, tb, d) = rim_carrying(&rims, pa, pb);
            if d > 1e-3 {
                return Err(format!(
                    "rim edge {d:.3e} off sector {s}'s rims (ends {pa:?}, {pb:?})"
                ));
            }
            brep.vertices[a].position = rims[rim].point_at(ta);
            brep.vertices[b].position = rims[rim].point_at(tb);
            brep.edges[edge].curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Nurbs(rims[rim].clone()),
                t_start: ta.min(tb),
                t_end: ta.max(tb),
            });
        }
    }
    Ok(faces)
}

/// Spoke `s` as its sector's `v = 0` iso-curve (o at u = 0, m_s at 1).
fn spoke_curve(sector: &NurbsSurface) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Nurbs(sector.iso_curve_at_v(0.0)),
        t_start: 0.0,
        t_end: 1.0,
    }
}

/// The rim iso-curve carrying both ends of a rim edge: the one whose
/// worse end distance is smallest, with both end parameters and that
/// distance.
fn rim_carrying(
    rims: &[popkern_core::geometry::NurbsCurve; 2],
    pa: Point3,
    pb: Point3,
) -> (usize, f64, f64, f64) {
    let fit = |k: usize| {
        let (ta, tb) = (param_on(&rims[k], pa), param_on(&rims[k], pb));
        let d = rims[k]
            .point_at(ta)
            .distance_to(&pa)
            .max(rims[k].point_at(tb).distance_to(&pb));
        (k, ta, tb, d)
    };
    let (first, second) = (fit(0), fit(1));
    if first.3 <= second.3 { first } else { second }
}

/// Parameter of the point of `curve` (on [0, 1]) nearest `p`: a dense
/// scan refined by ternary search, snapped to the ends within 1e-12.
fn param_on(curve: &popkern_core::geometry::NurbsCurve, p: Point3) -> f64 {
    let mut best = (f64::INFINITY, 0.0);
    for i in 0..=256 {
        let t = i as f64 / 256.0;
        let d = curve.point_at(t).distance_to(&p);
        if d < best.0 {
            best = (d, t);
        }
    }
    let (mut lo, mut hi) = (
        (best.1 - 1.0 / 256.0).max(0.0),
        (best.1 + 1.0 / 256.0).min(1.0),
    );
    for _ in 0..80 {
        let t1 = lo + (hi - lo) / 3.0;
        let t2 = hi - (hi - lo) / 3.0;
        if curve.point_at(t1).distance_to(&p) < curve.point_at(t2).distance_to(&p) {
            hi = t2;
        } else {
            lo = t1;
        }
    }
    let t = 0.5 * (lo + hi);
    if t < 1e-12 {
        0.0
    } else if t > 1.0 - 1e-12 {
        1.0
    } else {
        t
    }
}

/// Exact Linear iso pcurves on every cap half-edge: each vertex of a
/// sector face maps to its sector parameters (o → (0,0), a spoke point
/// → an edge of the unit square, rim points → u = 1 or v = 1).
fn author_cap_pcurves(
    brep: &mut BRep,
    cap_faces: &[FaceHandle],
    cap: &VariationalCap,
) -> Result<(), String> {
    for (s, &face) in cap_faces.iter().enumerate() {
        let sec = &cap.sectors[s];
        let rims = [sec.iso_curve_at_u(1.0), sec.iso_curve_at_v(1.0)];
        let spokes = [sec.iso_curve_at_v(0.0), sec.iso_curve_at_u(0.0)];
        let uv_of = |p: Point3| -> (f64, f64) {
            let cands = [
                (param_on(&rims[0], p), 0usize),
                (param_on(&rims[1], p), 1),
                (param_on(&spokes[0], p), 2),
                (param_on(&spokes[1], p), 3),
            ];
            let curves = [&rims[0], &rims[1], &spokes[0], &spokes[1]];
            let (t, k) = cands
                .iter()
                .copied()
                .min_by(|a, b| {
                    curves[a.1]
                        .point_at(a.0)
                        .distance_to(&p)
                        .total_cmp(&curves[b.1].point_at(b.0).distance_to(&p))
                })
                .unwrap_or((0.0, 0));
            match k {
                0 => (1.0, t),
                1 => (t, 1.0),
                2 => (t, 0.0),
                _ => (0.0, t),
            }
        };
        for he in face_loop_half_edges(brep, face)? {
            let a = brep.vertices[brep.half_edges[he].origin].position;
            let b = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
            let (uv_a, uv_b) = (uv_of(a), uv_of(b));
            brep.half_edges[he].pcurve = Some(Pcurve::linear(uv_a, uv_b, 0.0, 1.0));
        }
    }
    Ok(())
}
