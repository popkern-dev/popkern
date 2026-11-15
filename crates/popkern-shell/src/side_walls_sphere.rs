//! A wall's vertical edges on a sphere follow the sphere (PK-1216).
//!
//! A wall lying in a removed face's surface joins the outer rim to the
//! inner one; its vertical edges were straight lines, exact where the
//! wall is a plane, a cylinder or a cone (the rim vertex moves along a
//! ruling) but off a sphere: a capsule opened at its dome has its wall on
//! the dome's sphere from the rim at z = 10 to the cavity's at z = 13,
//! and the line between them left the sphere by 0.13 (an edge-precision
//! decline). On a sphere the vertical is the great-circle arc between
//! its ends, in the plane through them and the centre.

use std::collections::HashMap;

use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::circle_recognition::{CircleFit, arc_between};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::VerticalEdge;

/// Re-author every vertical edge whose wall lies on a sphere as the
/// great-circle arc between its ends.
pub(super) fn arc_verticals_on_spheres(
    builder: &mut BRepBuilder,
    verticals: &HashMap<(VertexHandle, VertexHandle), VerticalEdge>,
    vertical_faces: &HashMap<(VertexHandle, VertexHandle), FaceHandle>,
) {
    let lin = builder.brep().tolerance.linear();
    for (key, vertical) in verticals {
        let Some(face) = vertical_faces.get(key) else {
            continue;
        };
        let Some(SurfaceGeometry::Sphere { center, radius, .. }) = builder
            .brep()
            .faces
            .get(*face)
            .and_then(|f| f.surface.as_ref())
            .map(|s| s.geometry.clone())
        else {
            continue;
        };
        let brep = builder.brep();
        let (p_out, p_in) = (brep.vertices[key.0].position, brep.vertices[key.1].position);
        let Some(edge) = brep.half_edges.get(vertical.down).map(|h| h.edge) else {
            continue;
        };
        // Both ends on the sphere, apart, and not antipodal (the plane
        // through them and the centre must be determined).
        let (a, b) = (p_out - center, p_in - center);
        if (a.length() - radius).abs() > lin || (b.length() - radius).abs() > lin {
            continue;
        }
        let Some(normal) = a.cross(&b).normalized() else {
            continue;
        };
        let fit = CircleFit {
            center,
            normal,
            radius,
            max_deviation: 0.0,
            closed: false,
        };
        builder.set_edge_curve(edge, arc_between(&fit, normal, p_out, p_in));
    }
}
