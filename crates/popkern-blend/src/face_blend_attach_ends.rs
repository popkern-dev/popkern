//! Strip-end classification for the TrimBoth attach: the
//! closed-pocket end split and the PK-163 through-end recognition.
//! Split from `face_blend_attach_surgery.rs` (line ratchet).

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use popkern_core::geometry::Point3;

use crate::contact::ContactSolution;
use crate::error::BlendError;
use crate::face_blend_attach::lerp;
use crate::face_blend_attach_geometry::{
    decline, face_surface, face_uv_polygon, pick_doomed, point_in_polygon, project_uv,
    uniform_polyline_nurbs, uv_polyline_nurbs,
};

use super::face_blend_attach_surgery::reconcile_edge_curve;

pub(crate) struct EndSplit {
    pub(crate) kept: FaceHandle,
    pub(crate) doomed: FaceHandle,
    pub(crate) edge: EdgeHandle,
}

/// One strip end: closed-pocket split, or a THROUGH end whose arc
/// escaped the carrying face's trim (PK-163 cap extension).
pub(crate) enum EndOutcome {
    Closed(EndSplit),
    Through(crate::face_blend_attach_through::ThroughEnd),
}

pub(crate) fn split_end_face(
    work: &mut BRep,
    arc: &[Point3],
    section: &ContactSolution,
    exclude: &[FaceHandle; 4],
    tolerance: &Tolerance,
) -> Result<EndOutcome, BlendError> {
    let gate = tolerance.linear().max(1.0e-9) * 100.0;
    let mut best: Option<(FaceHandle, f64)> = None;
    let mut through: Option<(FaceHandle, f64)> = None;
    for (fh, f) in work.faces.iter() {
        if exclude.contains(&fh) {
            continue;
        }
        let Some(sof) = f.surface.as_ref() else {
            continue;
        };
        let mut max_d = 0.0_f64;
        let mut projectable = true;
        let mut arc_uv = Vec::with_capacity(arc.len());
        for p in arc {
            let Some(uv) = sof.geometry.closest_point(p) else {
                projectable = false;
                break;
            };
            max_d = max_d.max(sof.geometry.point_at(uv.0, uv.1).distance_to(p));
            arc_uv.push(uv);
        }
        if !(projectable && max_d < gate) {
            continue;
        }
        // The arc must lie inside the face's TRIM, not merely on its
        // surface — a through-gap's end arc sits on the cap plane but
        // over the void (interior samples outside the trim polygon;
        // endpoints legitimately sit ON the boundary, so only the
        // interior samples are tested).
        let Ok(polygon) = face_uv_polygon(work, fh) else {
            continue;
        };
        let contained = arc_uv[1..arc_uv.len().saturating_sub(1)]
            .iter()
            .all(|&uv| point_in_polygon(uv, &polygon));
        if contained && best.map(|(_, bd)| max_d < bd).unwrap_or(true) {
            best = Some((fh, max_d));
        }
        if !contained && through.map(|(_, bd)| max_d < bd).unwrap_or(true) {
            through = Some((fh, max_d));
        }
    }
    // THROUGH end (PK-163): the arc lies on a face's surface but
    // over the void — the material-adding cap extension takes over.
    if best.is_none()
        && let Some((end_face, _)) = through
    {
        return Ok(EndOutcome::Through(
            crate::face_blend_attach_through::ThroughEnd {
                end_face,
                arc: arc.to_vec(),
            },
        ));
    }
    let (cap, _) = best.ok_or_else(|| {
        decline(
            "end-arc support",
            format!(
                "no face carries the strip end arc inside its trim within {gate:.3e} \
                 (and none carries it on-surface for the through cap extension)"
            ),
        )
    })?;

    let surface = face_surface(work, cap)?.clone();
    let uv_pts = project_uv(&surface, arc)?;
    let uv_curve = uv_polyline_nurbs(&uv_pts);
    let curve3 = uniform_polyline_nurbs(arc);
    let split = work
        .split_face(cap, &uv_curve, &curve3, tolerance)
        .map_err(|e| decline("end split", format!("{e}")))?;
    let edge = *split
        .new_edges
        .first()
        .ok_or_else(|| decline("end split", "split_face produced no trim edge".into()))?;
    reconcile_edge_curve(work, edge, &curve3, tolerance);

    // Witness: anchor at the arc CHORD midpoint with the sag
    // direction — the doomed wedge's vertices all sit on the sag
    // side of the chord, while the arc-bottom anchor would put both
    // centroids on the same side (measured on the channel end wall).
    let mid = arc[arc.len() / 2];
    let chord_mid = lerp(arc[0], arc[arc.len() - 1], 0.5);
    let w3 = (mid - chord_mid).normalized().ok_or_else(|| {
        decline(
            "end trim side",
            "end arc has no sag — witness degenerate".into(),
        )
    })?;
    let uv_mid = surface.closest_point(&mid).ok_or_else(|| {
        decline(
            "end trim side",
            "arc midpoint does not project onto the end face".into(),
        )
    })?;
    let n = surface.normal_at(uv_mid.0, uv_mid.1);
    let w = (w3 - n * w3.dot(&n)).normalized().ok_or_else(|| {
        decline(
            "end trim side",
            "arc sag is normal to the end face — witness degenerate".into(),
        )
    })?;
    let (doomed, kept) = pick_doomed(
        work,
        split.face_a,
        split.face_b,
        chord_mid,
        w,
        "end trim side",
    )?;
    let _ = section;
    Ok(EndOutcome::Closed(EndSplit { kept, doomed, edge }))
}
