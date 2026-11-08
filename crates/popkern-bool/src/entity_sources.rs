//! Which operand entities a boolean result's face or edge is a piece
//! of (PK-633, PK-1121): the relation lineage records and naming
//! scheme 5 names by.
//!
//! Decided geometrically and independently of how the pipeline made
//! the entity: an interior point of the output face lies on an input
//! face's surface (within ten tolerances, normals parallel) and inside
//! its trim; an output edge's midpoint lies on an input edge's curve
//! within its range, tangents parallel. A fragment satisfies this
//! however much the boolean cut it, which a centroid or midpoint match
//! misses. Coincident faces of both operands both match.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::intersection::surface_subdivision::Aabb;
use popkern_core::naming::{EntityName, NameId};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, Edge, EdgeHandle, FaceHandle};

use crate::entity_lineage::{contains, foot, in_range, interior_point, lies_on, padded};

/// An operand entity a result entity is a piece of: the operand and
/// the entity's name id in that operand's pool.
pub(crate) type Source<'a> = (&'a BRep, NameId);

/// The names and ancestry, as values, of `sources` (the lineage a
/// result entity records).
pub(crate) fn lineage_of(sources: &[Source<'_>]) -> Vec<EntityName> {
    let mut lineage = Vec::new();
    for (input, id) in sources {
        if let Some(value) = input.name_pool.get(*id) {
            lineage.push(value.clone());
        }
        lineage.extend(input.name_pool.ancestry(*id));
    }
    lineage
}

/// The operand faces (of `a`, then each tool) with padded bounds, the
/// cheap reject before a projection.
pub(crate) struct FaceSources<'a> {
    inputs: Vec<(&'a BRep, Vec<(FaceHandle, Aabb)>)>,
    on: f64,
}

impl<'a> FaceSources<'a> {
    /// The faces of `a` and `tools`; `a`'s bounds from the carry index
    /// when it has them.
    pub(crate) fn new(
        a: &'a BRep,
        tools: &[&'a BRep],
        tol: &Tolerance,
        carry: Option<&crate::entity_lineage_carry::CarryIndex>,
    ) -> Self {
        // The names' own match tolerance: the boolean's, widened for the
        // projection round-off on curved inputs.
        let on = 10.0 * tol.linear();
        let boxes_of = |brep: &BRep| match carry {
            Some(c) if std::ptr::eq(brep, a) => c.a_face_boxes().to_vec(),
            _ => crate::intersect::compute_face_aabbs(brep),
        };
        let inputs = std::iter::once(a)
            .chain(tools.iter().copied())
            .map(|brep| {
                let boxes = boxes_of(brep)
                    .into_iter()
                    .map(|(f, bb)| (f, padded(&bb, on)))
                    .collect();
                (brep, boxes)
            })
            .collect();
        Self { inputs, on }
    }

    /// The named operand faces `face` of `result` is a piece of.
    pub(crate) fn of(&self, result: &BRep, face: FaceHandle) -> Vec<Source<'a>> {
        self.sensed(result, face)
            .into_iter()
            .map(|(s, _)| s)
            .collect()
    }

    /// [`FaceSources::of`], each with whether `face` faces the way its
    /// source does (a coincident face of the other side faces opposite).
    pub(crate) fn sensed(&self, result: &BRep, face: FaceHandle) -> Vec<(Source<'a>, bool)> {
        let Some((point, normal)) = interior_point(result, face) else {
            return Vec::new();
        };
        let mut out = Vec::new();
        for (input, boxes) in &self.inputs {
            for (candidate, aabb) in boxes {
                if !contains(aabb, &point) || !lies_on(input, *candidate, &point, &normal, self.on)
                {
                    continue;
                }
                if let Some(id) = input.faces.get(*candidate).and_then(|f| f.name) {
                    let aligned =
                        facing(input, *candidate, &point).is_none_or(|m| m.dot(&normal) > 0.0);
                    out.push(((*input, id), aligned));
                }
            }
        }
        out
    }
}

/// The operand edges (of `a`, then each tool) with padded curve bounds.
pub(crate) struct EdgeSources<'a> {
    inputs: Vec<(&'a BRep, Vec<(&'a Edge, Aabb)>)>,
    on: f64,
}

impl<'a> EdgeSources<'a> {
    /// The edges of `a` and `tools`.
    pub(crate) fn new(a: &'a BRep, tools: &[&'a BRep], tol: &Tolerance) -> Self {
        let on = 10.0 * tol.linear();
        let inputs = std::iter::once(a)
            .chain(tools.iter().copied())
            .map(|brep| {
                let boxes = brep
                    .edges
                    .values()
                    .filter_map(|e| {
                        let c = e.curve.as_ref()?;
                        let range = (c.t_start.min(c.t_end), c.t_start.max(c.t_end));
                        Some((e, padded(&c.geometry.curve_bounding_box(range), on)))
                    })
                    .collect();
                (brep, boxes)
            })
            .collect();
        Self { inputs, on }
    }

    /// The named operand edges `edge` of `result` is a piece of.
    pub(crate) fn of(&self, result: &BRep, edge: EdgeHandle) -> Vec<Source<'a>> {
        let Some(curve) = result.edges.get(edge).and_then(|e| e.curve.as_ref()) else {
            return Vec::new();
        };
        let t = 0.5 * (curve.t_start + curve.t_end);
        let point = curve.geometry.point_at(t);
        let raw = curve.geometry.tangent_at(t);
        let tangent = raw.normalized().unwrap_or(raw);
        let mut out = Vec::new();
        for (input, boxes) in &self.inputs {
            for (data, aabb) in boxes {
                if !contains(aabb, &point) {
                    continue;
                }
                let (Some(c), Some(id)) = (data.curve.as_ref(), data.name) else {
                    continue;
                };
                let Some(s) = foot(&c.geometry, &point) else {
                    continue;
                };
                if (c.geometry.point_at(s) - point).length() > self.on
                    || !in_range(&c.geometry, s, c.t_start, c.t_end)
                {
                    continue;
                }
                let along = c.geometry.tangent_at(s);
                let along = along.normalized().unwrap_or(along);
                if along.dot(&tangent).abs() < 1.0 - 1e-6 {
                    continue;
                }
                out.push((*input, id));
            }
        }
        out
    }
}

/// The outward normal of `face` of `input` at `p` (its orientation
/// applied).
fn facing(input: &BRep, face: FaceHandle, p: &popkern_core::geometry::Point3) -> Option<Vec3> {
    use popkern_core::geometry::surface::SurfaceEval;
    let surface = input.faces.get(face)?.surface.as_ref()?;
    let (u, v) = popkern_core::geometry::uv_projection::project_point_to_uv(&surface.geometry, p)?;
    let raw = surface.geometry.normal_at(u, v);
    let n = raw.normalized().unwrap_or(raw);
    Some(if surface.orientation { n } else { -n })
}
