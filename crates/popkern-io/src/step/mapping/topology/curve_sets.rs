//! `GEOMETRIC_CURVE_SET` and `GEOMETRIC_SET` items: wire edges and
//! acorns (PK-1114).
//!
//! A curve set lists bounded curves and points, as OCC writes free
//! edges and vertices (and `export_general` a general body's wires and
//! acorns). Each set imports as one body in the void: a wire edge per
//! curve between its end points (a closed curve as two edges, a wire
//! edge's ends being distinct vertices), a polyline's segments as line
//! edges, an acorn per point; end points within the linear tolerance are
//! one vertex. Surfaces in a `GEOMETRIC_SET` and unbounded lines are
//! skipped with a diagnostic.

use popkern_core::geometry::Point3;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::VertexHandle;
use popkern_core::topology::brep::BRep;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::mapping::geometry::GeometryCache;
use crate::step::parser::graph::EntityGraph;
use crate::step::parser::types::{EntityId, Parameter};

/// Build one wire-and-acorn body per curve set holding anything.
pub(super) fn map_curve_sets(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    log: &mut DiagnosticLog,
    tolerance: Tolerance,
) -> Vec<(EntityId, BRep)> {
    let mut ids: Vec<EntityId> = graph
        .iter()
        .filter(|(_, e)| {
            matches!(
                e.type_name.as_str(),
                "GEOMETRIC_CURVE_SET" | "GEOMETRIC_SET"
            )
        })
        .map(|(&id, _)| id)
        .collect();
    ids.sort_unstable();
    let mut results = Vec::new();
    for id in ids {
        let Some(Parameter::List(elements)) = graph.get(id).and_then(|e| e.params.get(1)) else {
            log.error(Phase::MapTopology, Some(id), "curve set without elements");
            continue;
        };
        let mut set = CurveSet {
            body: BRep::new(tolerance),
            vertices: Vec::new(),
            tolerance: tolerance.linear(),
        };
        for element in elements {
            let Parameter::EntityRef(element) = element else {
                continue;
            };
            set.add(graph, geo, *element, log);
        }
        if !set.body.vertices.is_empty() {
            results.push((id, set.body));
        }
    }
    results
}

/// A curve set's body under construction.
struct CurveSet {
    body: BRep,
    vertices: Vec<(Point3, VertexHandle)>,
    tolerance: f64,
}

impl CurveSet {
    fn add(
        &mut self,
        graph: &EntityGraph,
        geo: &mut GeometryCache<'_>,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) {
        let Some(entity) = graph.get(id) else {
            return;
        };
        match entity.type_name.as_str() {
            "CARTESIAN_POINT" => {
                if let Some(p) = geo.resolve_point(id, log) {
                    self.vertex(p);
                }
            }
            "POLYLINE" => {
                let Some(Parameter::List(points)) = entity.params.get(1) else {
                    return;
                };
                let points: Vec<Point3> = points
                    .iter()
                    .filter_map(|p| match p {
                        Parameter::EntityRef(p) => geo.resolve_point(*p, log),
                        _ => None,
                    })
                    .collect();
                for pair in points.windows(2) {
                    let line = CurveOnEdge {
                        geometry: CurveGeometry::Line {
                            origin: pair[0],
                            direction: pair[1] - pair[0],
                        },
                        t_start: 0.0,
                        t_end: 1.0,
                    };
                    self.edge(line, id, log);
                }
            }
            "TRIMMED_CURVE" => {
                if let Some(curve) = trimmed(graph, geo, id, &entity.params, log) {
                    self.edge(curve, id, log);
                }
            }
            // A line has no bounds of its own.
            "LINE" => log.warn(
                Phase::MapTopology,
                Some(id),
                "unbounded line in a curve set skipped",
            ),
            // Bounded curves (splines) and closed ones (circles, ellipses)
            // run over their own parameter range; anything else that is
            // not a curve (a surface in a geometric set) is skipped.
            _ => match geo.resolve_curve(id, log) {
                Some(curve) if (curve.t_end - curve.t_start).is_finite() => {
                    self.edge(curve, id, log)
                }
                _ => log.warn(
                    Phase::MapTopology,
                    Some(id),
                    format!("{} in a curve set skipped", entity.type_name),
                ),
            },
        }
    }

    /// The vertex at `p`, made an acorn when new.
    fn vertex(&mut self, p: Point3) -> Option<VertexHandle> {
        if let Some(&(_, v)) = self
            .vertices
            .iter()
            .find(|(q, _)| q.distance_to(&p) <= self.tolerance)
        {
            return Some(v);
        }
        let outside = self.body.outside;
        let v = self.body.make_acorn(outside, p).ok()?.vertex;
        self.vertices.push((p, v));
        Some(v)
    }

    /// A wire edge along `curve` (two when it closes on itself).
    fn edge(&mut self, curve: CurveOnEdge, id: EntityId, log: &mut DiagnosticLog) {
        let start = curve.geometry.point_at(curve.t_start);
        let end = curve.geometry.point_at(curve.t_end);
        if start.distance_to(&end) <= self.tolerance {
            let mid = 0.5 * (curve.t_start + curve.t_end);
            let (mut a, mut b) = (curve.clone(), curve);
            a.t_end = mid;
            b.t_start = mid;
            if a.geometry.point_at(mid).distance_to(&start) <= self.tolerance {
                log.warn(Phase::MapTopology, Some(id), "degenerate curve skipped");
                return;
            }
            self.edge(a, id, log);
            self.edge(b, id, log);
            return;
        }
        let (Some(v0), Some(v1)) = (self.vertex(start), self.vertex(end)) else {
            return;
        };
        let outside = self.body.outside;
        match self.body.make_wire_edge(outside, v0, v1) {
            Ok(made) => self.body.edges[made.edge].curve = Some(curve),
            Err(e) => log.warn(
                Phase::MapTopology,
                Some(id),
                format!("wire edge skipped: {e}"),
            ),
        }
    }
}

/// A `TRIMMED_CURVE`'s basis trimmed to its run: each trim a point
/// (preferred, projected onto the basis) or a parameter value, the run
/// from the first to the second along the curve's sense.
fn trimmed(
    graph: &EntityGraph,
    geo: &mut GeometryCache<'_>,
    id: EntityId,
    params: &[Parameter],
    log: &mut DiagnosticLog,
) -> Option<CurveOnEdge> {
    let Some(Parameter::EntityRef(basis)) = params.get(1) else {
        return None;
    };
    let mut curve = geo.resolve_curve(*basis, log)?;
    let mut trim = |p: Option<&Parameter>| -> Option<f64> {
        let Some(Parameter::List(select)) = p else {
            return None;
        };
        let point = select.iter().find_map(|s| match s {
            Parameter::EntityRef(r)
                if graph
                    .get(*r)
                    .is_some_and(|e| e.type_name == "CARTESIAN_POINT") =>
            {
                Some(*r)
            }
            _ => None,
        });
        if let Some(point) = point {
            let p = geo.resolve_point(point, log)?;
            return curve.geometry.closest_point(&p);
        }
        select.iter().find_map(|s| match s {
            Parameter::Typed(name, value) if name == "PARAMETER_VALUE" => match **value {
                Parameter::Real(t) => Some(t),
                Parameter::Int(t) => Some(t as f64),
                _ => None,
            },
            _ => None,
        })
    };
    let (Some(mut t0), Some(mut t1)) = (trim(params.get(2)), trim(params.get(3))) else {
        log.warn(Phase::MapTopology, Some(id), "TRIMMED_CURVE trims not read");
        return None;
    };
    let sense = !matches!(params.get(4), Some(Parameter::Enum(e)) if e == "F");
    super::edge_trim::normalize_periodic_parameter_order(&curve.geometry, sense, &mut t0, &mut t1);
    curve.t_start = t0;
    curve.t_end = t1;
    Some(curve)
}
