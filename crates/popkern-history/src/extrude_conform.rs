//! Extrude ends that meet a face not square to the extrude (brep-bench
//! §4.4): a pad run down to a drum, a post up to a dome's inner face.
//!
//! The end of an up-to-surface extrude used to be one distance along
//! the ray from the section's centroid, so it came out flat: under a
//! curved (or tilted) face it stopped short at the section's rim or
//! ran into the face's body at its middle. Here every point of the
//! section casts its own ray. When the hits agree to the tolerance the
//! face is square to the extrude and the flat end at the centroid's
//! distance IS the face (planes, and a face the section meets flat).
//!
//! Otherwise a boss into the face's own body extrudes past the
//! farthest hit: what lies beyond the face is inside that body, so the
//! union ends exactly on the face. That claim is checked, not assumed:
//! points beyond the face along every ray must lie in the body (a
//! thin wall behind the face would let the boss come through). Every
//! other mode (a new body, a cut) and any non-zero offset would need
//! the tool trimmed by the face's surface; they decline with a typed
//! failure rather than serve the flat end.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::error::ReplayError;
use crate::failure::FeatureFailure;
use crate::feature::{ExtrudeEnd, OperationMode};
use crate::parameter::ParameterTable;
use crate::resolve::{resolve_extrude_distance, resolve_param};

pub(crate) fn declined(reason: &str, near: f64, far: f64) -> ReplayError {
    ReplayError::FeatureFailed {
        label: "extrude".into(),
        error: Box::new(FeatureFailure::UpToFaceDeclined {
            reason: reason.into(),
            near,
            far,
        }),
    }
}

/// At most this many section points cast a ray (each is a surface
/// intersection); a larger mesh is thinned evenly.
const MAX_SECTION_POINTS: usize = 128;

/// Points of the section the end is checked from: the profile face's
/// loop vertices and its mesh vertices (boundary and interior),
/// thinned evenly to [`MAX_SECTION_POINTS`].
pub(crate) fn section_points(brep: &BRep, face: FaceHandle) -> Vec<Point3> {
    let mut points = brep.face_outer_loop_positions(face);
    let options = popkern_mesh::TessellationOptions::default();
    if let Ok(mesh) = popkern_mesh::tessellate_face(brep, face, &options) {
        let step = mesh.vertices.len().div_ceil(MAX_SECTION_POINTS).max(1);
        points.extend(mesh.vertices.iter().step_by(step).map(|v| v.position));
    }
    points
}

/// The distance the primary end of an extrude of `face` along `dir`
/// runs to. Ends other than up-to / offset-from a surface resolve as
/// before; those do too when the face is square to the extrude.
#[allow(clippy::too_many_arguments)]
pub(crate) fn primary_end_distance(
    end: &ExtrudeEnd,
    mode: OperationMode,
    brep: &BRep,
    face: FaceHandle,
    centroid: Point3,
    dir: Vec3,
    host: Option<&BRep>,
    param_table: &ParameterTable,
    tolerance: &Tolerance,
) -> Result<f64, ReplayError> {
    let at_centroid =
        resolve_extrude_distance(end, centroid, dir, host, brep, param_table, tolerance)?;
    let (target, offset) = match end {
        ExtrudeEnd::UpToSurface { face } => (face, 0.0),
        ExtrudeEnd::OffsetFromSurface { face, offset } => {
            (face, resolve_param(param_table, *offset)?)
        }
        _ => return Ok(at_centroid),
    };
    let points = section_points(brep, face);
    let target_brep = host.unwrap_or(brep);
    let target_face = crate::resolve::resolve_face(target, target_brep)?;
    let surface = target_brep.faces[target_face]
        .surface
        .as_ref()
        .map(|s| &s.geometry)
        .ok_or_else(|| declined("face_out_of_reach", at_centroid, at_centroid))?;
    // The end is the face's SURFACE (a section wider than a planar
    // face stops on its plane): the rays search a parameter range
    // covering the face and where the section lands on it.
    let mut spread = target_brep.face_outer_loop_positions(target_face);
    spread.extend(points.iter().map(|p| *p + dir * at_centroid));
    let range = popkern_core::geometry::uv_projection::bounded_parameter_range(&spread, surface);
    let mut hits = Vec::with_capacity(points.len());
    for p in &points {
        match ray_hit(surface, range, *p, dir, tolerance)? {
            Some(d) => hits.push(d),
            None => return Err(declined("face_out_of_reach", at_centroid, at_centroid)),
        }
    }
    let near = hits.iter().copied().fold(f64::INFINITY, f64::min);
    let far = hits.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if far - near <= tolerance.linear() {
        return Ok(at_centroid);
    }
    let Some(body) = host else {
        return Err(declined("curved_end_needs_boss", near, far));
    };
    if !matches!(mode, OperationMode::Boss) || offset != 0.0 {
        return Err(declined("curved_end_needs_boss", near, far));
    }
    // Past the farthest hit by a share of the spread: the mesh's
    // points only sample the face's depth, and an end short of it
    // would leave a void between the cap and the face.
    let end_at = far + 0.05 * (far - near) + 10.0 * tolerance.linear();
    let boxes = popkern_bool::face_bounding_boxes(body);
    let inside = |q: Point3| {
        !matches!(
            popkern_bool::point_in_solid(&q, body, &boxes, tolerance),
            Ok(popkern_bool::PointClassification::Outside) | Err(_)
        )
    };
    for (p, d) in points.iter().zip(&hits) {
        let beyond = [(d + end_at) / 2.0, end_at];
        if beyond
            .iter()
            .any(|t| *t > d + tolerance.linear() && !inside(*p + dir * *t))
        {
            return Err(declined("thin_behind_face", near, far));
        }
    }
    Ok(end_at)
}

/// The nearest forward distance from `origin` along `dir` to `surface`
/// within `range`.
fn ray_hit(
    surface: &popkern_core::geometry::surface::SurfaceGeometry,
    range: ((f64, f64), (f64, f64)),
    origin: Point3,
    dir: Vec3,
    tolerance: &Tolerance,
) -> Result<Option<f64>, ReplayError> {
    let ray = popkern_core::geometry::CurveGeometry::Line {
        origin,
        direction: dir,
    };
    let result = popkern_core::geometry::intersection::intersect_surface_curve(
        surface,
        range.0,
        range.1,
        &ray,
        (0.0, 1e6),
        tolerance,
    )?;
    Ok(result
        .intersections
        .iter()
        .filter_map(|i| match i {
            popkern_core::geometry::intersection::SurfaceCurveIntersection::Point { t, .. } => {
                Some(*t)
            }
            _ => None,
        })
        .filter(|t| *t > tolerance.linear())
        .fold(None, |best: Option<f64>, t| {
            Some(best.map_or(t, |b| b.min(t)))
        }))
}
