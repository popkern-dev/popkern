//! Certified continuity of a BRep edge seam (PK-499): resolves the
//! edge's two faces, their surfaces and half-edge pcurves, orients
//! each pcurve to the edge curve, and hands the pair to
//! [`certify_seam`]. Topology references geometry here, never the
//! reverse (the certifier itself is geometry-only).

use crate::geometry::certify::{
    CertificateStatus, CertifyDecline, InconclusiveReason, SeamCertificate, SeamCertifyOptions,
    SeamSideSource, certify_seam,
};
use crate::geometry::curve::{CurveEval, CurveOnEdge};
use crate::geometry::surface::SurfaceEval;
use crate::topology::{BRep, EdgeHandle, HalfEdgeHandle};

/// Why an edge could not be presented to the certifier as a seam.
fn not_a_seam() -> SeamCertificate {
    SeamCertificate::unbounded()
}

/// One side of the seam: the half-edge's face surface and pcurve,
/// with the pcurve's orientation against the edge curve resolved by
/// matching its end images to the curve's start point.
fn side_source(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    curve: &CurveOnEdge,
) -> Result<SeamSideSource, Option<CertifyDecline>> {
    let data = brep.half_edges.get(half_edge).ok_or(None)?;
    let pcurve = data.pcurve.as_ref().ok_or(None)?;
    let face = brep.loops.get(data.loop_).ok_or(None)?.face;
    let surface = &brep
        .faces
        .get(face)
        .ok_or(None)?
        .surface
        .as_ref()
        .ok_or(None)?
        .geometry;
    let start = (&curve.geometry as &dyn CurveEval).point_at(curve.t_start);
    let (uv0, uv1) = (
        pcurve.point_at(pcurve.t_start),
        pcurve.point_at(pcurve.t_end),
    );
    let reversed = surface.point_at(uv1.x, uv1.y).distance_to(&start)
        < surface.point_at(uv0.x, uv0.y).distance_to(&start);
    SeamSideSource::new(surface, pcurve, reversed).map_err(Some)
}

/// Certify the seam `edge` of `brep`. Boundary edges, edges whose two
/// half-edges lie on one face, and edges missing a curve, pcurve or
/// surface are `Inconclusive(NotASeam)`; sides whose geometry cannot
/// be enclosed are `Inconclusive(Side(..))`.
pub fn certify_edge_seam(
    brep: &BRep,
    edge: EdgeHandle,
    options: &SeamCertifyOptions,
) -> SeamCertificate {
    let Some(data) = brep.edges.get(edge) else {
        return not_a_seam();
    };
    let Some(curve) = data.curve.as_ref() else {
        return not_a_seam();
    };
    let half_edge = data.half_edge;
    let Some(twin) = brep.half_edges.get(half_edge).map(|h| h.twin) else {
        return not_a_seam();
    };
    let face_of = |h: HalfEdgeHandle| {
        brep.half_edges
            .get(h)
            .and_then(|d| brep.loops.get(d.loop_))
            .map(|l| l.face)
    };
    match (face_of(half_edge), face_of(twin)) {
        (Some(a), Some(b)) if a != b => {}
        _ => return not_a_seam(),
    }
    let sides = (
        side_source(brep, half_edge, curve),
        side_source(brep, twin, curve),
    );
    let declined = |index: usize, decline: CertifyDecline| {
        let mut out = SeamCertificate::unbounded();
        out.status = CertificateStatus::Inconclusive {
            reason: InconclusiveReason::Side(index, decline),
            at: (0.0, 1.0),
        };
        out
    };
    match sides {
        (Ok(mut a), Ok(mut b)) => certify_seam(&mut a, &mut b, options),
        (Err(Some(decline)), _) => declined(0, decline),
        (_, Err(Some(decline))) => declined(1, decline),
        _ => not_a_seam(),
    }
}

/// Certification options for `brep`: characteristic length from its
/// vertex bounding-box diagonal (floor 1), correspondence and gap from
/// its own linear tolerance (see [`SeamCertifyOptions::for_model`]).
pub fn seam_certify_options(brep: &BRep) -> SeamCertifyOptions {
    let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
    for (_, vertex) in brep.vertices.iter() {
        let p = [vertex.position.x, vertex.position.y, vertex.position.z];
        for k in 0..3 {
            lo[k] = lo[k].min(p[k]);
            hi[k] = hi[k].max(p[k]);
        }
    }
    let diagonal = (0..3).map(|k| (hi[k] - lo[k]).powi(2)).sum::<f64>().sqrt();
    let scale = if diagonal.is_finite() {
        diagonal.max(1.0)
    } else {
        1.0
    };
    SeamCertifyOptions::for_model(scale, brep.tolerance.linear())
}
