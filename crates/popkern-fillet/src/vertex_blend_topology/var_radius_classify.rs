use super::*;

/// One cap face's outer-loop classification, used to identify
/// which edges need to be re-routed into cap-vs-cap form.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct CapLoopClassification {
    /// The cap face this classification describes.
    pub(crate) cap_face: FaceHandle,
    /// Per outer-loop half-edge in walk order: the half-edge
    /// handle, its origin vertex, the edge handle, and the
    /// twin-face's classification.
    pub(crate) entries: Vec<CapLoopEntry>,
}

/// Single outer-loop entry on a cap face.
#[allow(dead_code)]
#[derive(Debug, Clone)]
pub(crate) struct CapLoopEntry {
    /// Half-edge in the cap face's outer loop.
    pub(crate) half_edge: HalfEdgeHandle,
    /// Origin vertex of `half_edge`.
    pub(crate) origin: VertexHandle,
    /// Underlying edge.
    pub(crate) edge: EdgeHandle,
    /// Face on the *opposite* side of `edge`.
    pub(crate) twin_face: FaceHandle,
    /// Classification of `twin_face` for surgery purposes.
    pub(crate) twin_kind: TwinKind,
}

/// What lies on the *other* side of a cap face's outer-loop edge.
#[allow(dead_code)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum TwinKind {
    /// The cap's own per-fillet cylinder (the "good" cap-vs-
    /// cylinder seam — this edge stays as is, it IS the cap-tangent
    /// arc).
    OwnCylinder,
    /// A cylinder belonging to a *neighbour* fillet (the "bad"
    /// cap-vs-cylinder seam — surgery must re-route this).
    NeighbourCylinder,
    /// Another cap face (expected: cap-vs-cap inner spoke).
    Cap,
    /// Anything else (cube planar face, NURBS, …) — surgery
    /// declines on these.
    Other,
}

/// Detect each cap sub-face's outer-loop classification.
///
/// Reads `cap_faces` and `per_fillet_caps` (from the substituted
/// `BlendBoundary`) and returns one `CapLoopClassification` per
/// cap face, in the same order as `cap_faces`.  Pure read of the
/// BRep — no mutation.
#[allow(dead_code)]
pub(crate) fn classify_cap_loops(
    brep: &BRep,
    cap_faces: &[FaceHandle],
    per_fillet_caps: &[(popkern_core::geometry::Point3, f64)],
) -> Vec<CapLoopClassification> {
    use popkern_core::geometry::surface::SurfaceGeometry;

    let mut out: Vec<CapLoopClassification> = Vec::with_capacity(cap_faces.len());
    let cap_set: std::collections::BTreeSet<FaceHandle> = cap_faces.iter().copied().collect();

    for (k, &cap_face) in cap_faces.iter().enumerate() {
        let Some(face_data) = brep.faces.get(cap_face) else {
            continue;
        };
        let outer_loop = face_data.outer_loop;
        let Some(loop_data) = brep.loops.get(outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;

        // The cap's "own" sphere — assigned by the sphere-cap rung
        // from per_fillet_caps[k].  If the cap has a different
        // surface (e.g., the rung wasn't enabled), classification
        // still works topologically but `OwnCylinder` falls back
        // to a position-based distance check.
        let (own_centre, own_radius) = per_fillet_caps
            .get(k)
            .copied()
            .unwrap_or((popkern_core::geometry::Point3::new(0.0, 0.0, 0.0), 0.0));

        let mut entries: Vec<CapLoopEntry> = Vec::new();
        let mut cursor = start;
        let mut safety = 16usize;
        loop {
            if safety == 0 {
                break;
            }
            safety -= 1;
            let Some(he) = brep.half_edges.get(cursor) else {
                break;
            };
            let origin = he.origin;
            let edge = he.edge;
            let twin = he.twin;
            let next = he.next;
            let twin_face = brep
                .half_edges
                .get(twin)
                .and_then(|t| brep.loops.get(t.loop_))
                .map(|l| l.face)
                .unwrap_or(cap_face);

            let twin_kind = if cap_set.contains(&twin_face) {
                TwinKind::Cap
            } else {
                let twin_surf = brep
                    .faces
                    .get(twin_face)
                    .and_then(|f| f.surface.as_ref())
                    .map(|s| &s.geometry);
                match twin_surf {
                    Some(SurfaceGeometry::Cylinder {
                        origin: cyl_o,
                        axis,
                        radius,
                    }) => {
                        // OwnCylinder iff the cap's sphere is tangent
                        // to this cylinder along the cap-tangent
                        // circle: radius matches AND sphere centre
                        // projects onto the cylinder axis.
                        let to_centre = popkern_core::geometry::Vec3::new(
                            own_centre.x - cyl_o.x,
                            own_centre.y - cyl_o.y,
                            own_centre.z - cyl_o.z,
                        );
                        let axis_unit = axis.normalized().unwrap_or(*axis);
                        let along = to_centre.dot(&axis_unit);
                        let projection = popkern_core::geometry::Point3::new(
                            cyl_o.x + axis_unit.x * along,
                            cyl_o.y + axis_unit.y * along,
                            cyl_o.z + axis_unit.z * along,
                        );
                        let radial = own_centre.distance_to(&projection);
                        // Tangency conditions: radial ≈ 0 (centre
                        // on axis) AND own_radius ≈ cylinder.radius.
                        let radial_ok = radial < own_radius * 0.05;
                        let radius_ok = (own_radius - *radius).abs() < own_radius * 0.05;
                        if radial_ok && radius_ok {
                            TwinKind::OwnCylinder
                        } else {
                            TwinKind::NeighbourCylinder
                        }
                    }
                    _ => TwinKind::Other,
                }
            };

            entries.push(CapLoopEntry {
                half_edge: cursor,
                origin,
                edge,
                twin_face,
                twin_kind,
            });

            cursor = next;
            if cursor == start {
                break;
            }
        }

        out.push(CapLoopClassification { cap_face, entries });
    }
    out
}
