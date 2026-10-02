//! The corner-patched variational cap on a carved setback region
//! (PK-619): the core's `fit_corner_patched_cap` gives a main cap and one
//! four-sided corner patch per thin side, joined along seams K. Here the
//! blend face is cut along each seam, each region is subdivided into its
//! piece's sectors, the seams take the main cap's rim iso-curves (each
//! side keeps its own piece's exact iso pcurve), and every
//! seam touching the cap, the K seams included, is certified. Any
//! failure is a cap-stage decline; the rung restores the snapshot.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::variational_cap::{
    CornerPatchOptions, VariationalCapOptions, VariationalCapSide, VariationalCapTargets,
    fit_corner_patched_cap,
};
use popkern_core::geometry::{CurveGeometry, CurveOnEdge};

use super::*;

/// Fit, build and certify the corner-patched cap over `blend_face`.
pub(super) fn commit_patched_fan(
    brep: &mut BRep,
    blend_face: FaceHandle,
    fillet_surfaces: &[FilletSurfaceResult],
    sides: &[VariationalCapSide],
    (targets, options): (&VariationalCapTargets, &VariationalCapOptions),
    tolerance: &Tolerance,
) -> Result<(Vec<FaceHandle>, FitServe), CapCommitError> {
    let gates = ContinuityGates::commercial();
    let mut cap = fit_corner_patched_cap(
        sides,
        targets,
        options,
        &CornerPatchOptions::default(),
        tolerance.linear(),
    )
    .map_err(|e| format!("variational fit: {e}"))?;
    if crate::dispatch_overrides::fit_cap_rim_tilted() {
        tilt_rim_row(&mut cap.main.cap.sectors[0]);
        tilted_face(&mut cap.main.cap.faces[0], &cap.main.cap.sectors[0]);
    }
    // 1. Cut the region along each seam; the rest stays the main region.
    let mut main_region = blend_face;
    let mut patch_regions = Vec::with_capacity(cap.corners.len());
    for corner in &cap.corners {
        let (patch_region, rest) = cut_along_seam(brep, main_region, corner, sides)?;
        patch_regions.push(patch_region);
        main_region = rest;
    }
    // 2. Subdivide the patches, then the main cap: the seams end on the
    //    main cap's rim iso-curves.
    let mut cap_faces = Vec::new();
    let mut pieces = Vec::new();
    let mut patch_faces = Vec::with_capacity(cap.corners.len());
    for (k, (corner, &region)) in cap.corners.iter().zip(&patch_regions).enumerate() {
        let faces = build_sector_faces(brep, region, &corner.piece.cap, &corner.piece.sides)?;
        pieces.extend(face_pieces(k + 1, faces.len()));
        patch_faces.push(faces.clone());
        cap_faces.extend(faces);
    }
    let main_faces = build_sector_faces(brep, main_region, &cap.main.cap, &cap.main.sides)?;
    pieces.extend(face_pieces(0, main_faces.len()));
    cap_faces.extend(main_faces.iter().copied());
    // 3. Pcurves: every cap half-edge, the seams' included, takes its own
    //    piece's exact iso pcurve. A seam half-edge co-parameterised
    //    onto a patch instead reaches the patch's dual-jet corner off
    //    its iso, where the rational blend is 0/0 and the certifier's
    //    boundary surrogate cannot read it (measured: every seam
    //    inconclusive, `DenominatorContainsZero` at the seam's end).
    //    The outer neighbours are co-parameterised against the rims.
    for (corner, faces) in cap.corners.iter().zip(&patch_faces) {
        author_cap_pcurves(brep, faces, &corner.piece.cap)?;
    }
    author_cap_pcurves(brep, &main_faces, &cap.main.cap)?;
    super::super::author_setback_pcurves(brep, &[], &cap_faces)?;
    author_neighbour_pcurves_outside(brep, &cap_faces, &cap_faces, tolerance)?;
    super::super::setback_boundary_tolerance::record_boundary_tolerances(brep, &cap_faces);
    // 4. Certify every seam touching the cap; the record carries the
    //    main cap's report with the thinnest piece's regularity.
    let mut report = cap.main.cap.report.clone();
    report.min_regularity = cap
        .corners
        .iter()
        .map(|c| c.piece.cap.report.min_regularity)
        .fold(report.min_regularity, f64::min);
    let certified = setback_fit_certify::certify_cap(brep, &cap_faces, &gates, &report)?;
    let mode = SetbackFitMode::TangentCornerPatched;
    let required = mode.level();
    if level_rank(certified.level) < level_rank(required) {
        return Err(CapCommitError {
            reason: format!(
                "corner-patched fit certified {:?}, needs {required:?}: rim (g0 {:.3e}, g1 {:.3e}), seams and spokes (g0 {:.3e}, g1 {:.3e})",
                certified.level,
                certified.rim[0],
                certified.rim[1],
                certified.spoke[0],
                certified.spoke[1],
            ),
            fit_certified: Some(certified),
        });
    }
    if let Err(intersection) = fit_cap_region_self_intersection(brep, &cap_faces, fillet_surfaces) {
        return Err(CapCommitError {
            reason: format!("fit cap self-intersection postcondition failed: {intersection}"),
            fit_certified: Some(certified),
        });
    }
    Ok((
        cap_faces,
        FitServe {
            mode,
            certified,
            pieces,
        },
    ))
}

/// Piece `piece`'s faces, in sector order.
fn face_pieces(piece: usize, sectors: usize) -> Vec<FacePiece> {
    (0..sectors)
        .map(|sector| FacePiece {
            piece,
            sector,
            sectors,
        })
        .collect()
}

/// Split `region`'s boundary at the seam's ends and cut it along the
/// seam with `mef`. Returns (the patch region, the rest): the patch
/// region is the one whose boundary carries the thin side.
fn cut_along_seam(
    brep: &mut BRep,
    region: FaceHandle,
    corner: &popkern_core::geometry::variational_cap::CornerPatch,
    sides: &[VariationalCapSide],
) -> Result<(FaceHandle, FaceHandle), String> {
    let seam = &corner.seam;
    let (s0, s1) = seam.parameter_range();
    let v_b = split_blend_boundary_at(brep, region, seam.point_at(s0))?;
    let v_a = split_blend_boundary_at(brep, region, seam.point_at(s1))?;
    let h_b = loop_half_edge_from(brep, region, v_b)?;
    let h_a = loop_half_edge_from(brep, region, v_a)?;
    let mef = brep
        .mef(h_b, h_a)
        .map_err(|e| format!("mef seam over side {}: {e}", corner.side))?;
    brep.vertices[v_b].position = seam.point_at(s0);
    brep.vertices[v_a].position = seam.point_at(s1);
    brep.edges[mef.edge].curve = Some(CurveOnEdge {
        geometry: CurveGeometry::Nurbs(seam.clone()),
        t_start: s0,
        t_end: s1,
    });
    let he = brep.edges[mef.edge].half_edge;
    let twin = brep.half_edges[he].twin;
    let faces = [
        brep.loops[brep.half_edges[he].loop_].face,
        brep.loops[brep.half_edges[twin].loop_].face,
    ];
    let thin = &sides[corner.side].curve;
    let (t0, t1) = thin.parameter_range();
    let thin_mid = thin.point_at(0.5 * (t0 + t1));
    let carries_thin = |face: FaceHandle| -> Result<bool, String> {
        Ok(face_loop_half_edges(brep, face)?.into_iter().any(|h| {
            brep.edges[brep.half_edges[h].edge]
                .curve
                .as_ref()
                .is_some_and(|c| curve_distance(&c.geometry, c.t_start, c.t_end, thin_mid) < 1e-4)
        }))
    };
    match (carries_thin(faces[0])?, carries_thin(faces[1])?) {
        (true, false) => Ok((faces[0], faces[1])),
        (false, true) => Ok((faces[1], faces[0])),
        other => Err(format!(
            "seam over side {}: thin side on both or neither region ({other:?})",
            corner.side
        )),
    }
}

/// The half-edge leaving `vertex` on `face`'s loop.
fn loop_half_edge_from(
    brep: &BRep,
    face: FaceHandle,
    vertex: VertexHandle,
) -> Result<HalfEdgeHandle, String> {
    face_loop_half_edges(brep, face)?
        .into_iter()
        .find(|&h| brep.half_edges[h].origin == vertex)
        .ok_or_else(|| "seam end not on the region's loop".to_string())
}

/// `author_neighbour_pcurves` for `faces`, skipping twins on any face of
/// `cap`: the outer neighbours only.
fn author_neighbour_pcurves_outside(
    brep: &mut BRep,
    faces: &[FaceHandle],
    cap: &[FaceHandle],
    tolerance: &Tolerance,
) -> Result<(), String> {
    let gate = 10.0 * tolerance.linear();
    for &face in faces {
        for he in face_loop_half_edges(brep, face)? {
            let twin = brep.half_edges[he].twin;
            if cap.contains(&brep.loops[brep.half_edges[twin].loop_].face) {
                continue;
            }
            let edge = brep.half_edges[he].edge;
            let curve: CurveOnEdge = brep.edges[edge]
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
