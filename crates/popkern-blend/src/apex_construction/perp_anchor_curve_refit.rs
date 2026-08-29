//! Curve refit helpers for perp-anchor relocation.

use popkern_core::geometry::curve::CurveOnEdge;
use popkern_core::geometry::{CurveGeometry, Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

pub(super) fn relocate_vertex_with_curve_coherence(
    brep: &mut BRep,
    vertex: VertexHandle,
    target: Point3,
) -> Result<(), String> {
    let mut incidents: Vec<(EdgeHandle, EndChoice)> = Vec::new();
    for (eh, edge) in &brep.edges {
        let primary = edge.half_edge;
        let Some(he) = brep.half_edges.get(primary) else {
            continue;
        };
        let twin_he = brep
            .half_edges
            .get(he.twin)
            .ok_or_else(|| format!("edge {eh:?}: twin half-edge missing"))?;
        if he.origin == vertex {
            incidents.push((eh, EndChoice::PrimaryOrigin));
        } else if twin_he.origin == vertex {
            incidents.push((eh, EndChoice::TwinOrigin));
        }
    }

    brep.vertices
        .get_mut(vertex)
        .ok_or_else(|| format!("vertex {vertex:?} disappeared mid-relocate"))?
        .position = target;

    for (eh, end_choice) in incidents {
        refit_edge_curve_for_relocated_endpoint(brep, eh, end_choice)?;
    }
    Ok(())
}

#[derive(Debug, Clone, Copy)]
enum EndChoice {
    PrimaryOrigin,
    TwinOrigin,
}

fn refit_edge_curve_for_relocated_endpoint(
    brep: &mut BRep,
    edge: EdgeHandle,
    end_choice: EndChoice,
) -> Result<(), String> {
    let Some(edge_data) = brep.edges.get(edge) else {
        return Ok(());
    };
    let primary = edge_data.half_edge;
    let Some(he) = brep.half_edges.get(primary) else {
        return Ok(());
    };
    let v_origin = he.origin;
    let twin = he.twin;
    let Some(twin_he) = brep.half_edges.get(twin) else {
        return Ok(());
    };
    let v_dest = twin_he.origin;
    let p_origin = brep.vertices[v_origin].position;
    let p_dest = brep.vertices[v_dest].position;
    let _ = end_choice;

    let existing = brep.edges.get(edge).and_then(|e| e.curve.clone());
    let new_curve = match existing {
        None => return Ok(()),
        Some(coe) => refit_curve_through_endpoints(&coe, p_origin, p_dest),
    };

    if let Some(edge_data) = brep.edges.get_mut(edge) {
        edge_data.curve = Some(new_curve);
    }
    Ok(())
}

fn refit_curve_through_endpoints(coe: &CurveOnEdge, p0: Point3, p1: Point3) -> CurveOnEdge {
    match &coe.geometry {
        CurveGeometry::Line { .. } => line_curve_on_edge(p0, p1),
        CurveGeometry::Arc {
            center,
            normal,
            radius,
            ..
        } => {
            let (start_angle, end_angle) = arc_endpoint_angles(*center, *normal, p0, p1);
            CurveOnEdge {
                geometry: CurveGeometry::Arc {
                    center: *center,
                    normal: *normal,
                    radius: *radius,
                    start_angle,
                    end_angle,
                },
                t_start: start_angle,
                t_end: end_angle,
            }
        }
        CurveGeometry::Circle {
            center,
            normal,
            radius,
        } => {
            let (start_angle, end_angle) = arc_endpoint_angles(*center, *normal, p0, p1);
            CurveOnEdge {
                geometry: CurveGeometry::Arc {
                    center: *center,
                    normal: *normal,
                    radius: *radius,
                    start_angle,
                    end_angle,
                },
                t_start: start_angle,
                t_end: end_angle,
            }
        }
        CurveGeometry::Ellipse {
            center,
            u_basis,
            v_basis,
            ..
        } => {
            let (start_angle, end_angle) =
                ellipse_endpoint_angles(*center, *u_basis, *v_basis, p0, p1);
            CurveOnEdge {
                geometry: CurveGeometry::Ellipse {
                    center: *center,
                    u_basis: *u_basis,
                    v_basis: *v_basis,
                    start_angle,
                    end_angle,
                },
                t_start: start_angle,
                t_end: end_angle,
            }
        }
        CurveGeometry::Parabola { .. }
        | CurveGeometry::Hyperbola { .. }
        | CurveGeometry::Nurbs(_) => line_curve_on_edge(p0, p1),
    }
}

fn line_curve_on_edge(p0: Point3, p1: Point3) -> CurveOnEdge {
    CurveOnEdge {
        geometry: CurveGeometry::Line {
            origin: p0,
            direction: Vec3::new(p1.x - p0.x, p1.y - p0.y, p1.z - p0.z),
        },
        t_start: 0.0,
        t_end: 1.0,
    }
}

fn arc_endpoint_angles(center: Point3, normal: Vec3, p0: Point3, p1: Point3) -> (f64, f64) {
    let (u, v) = in_plane_basis(normal);
    let pa = Vec3::new(p0.x - center.x, p0.y - center.y, p0.z - center.z);
    let pb = Vec3::new(p1.x - center.x, p1.y - center.y, p1.z - center.z);
    let a0 = pa.dot(&v).atan2(pa.dot(&u));
    let a1 = pb.dot(&v).atan2(pb.dot(&u));
    (a0, a1)
}

fn ellipse_endpoint_angles(
    center: Point3,
    u_basis: Vec3,
    v_basis: Vec3,
    p0: Point3,
    p1: Point3,
) -> (f64, f64) {
    let pa = Vec3::new(p0.x - center.x, p0.y - center.y, p0.z - center.z);
    let pb = Vec3::new(p1.x - center.x, p1.y - center.y, p1.z - center.z);
    let u_len2 = u_basis.dot(&u_basis);
    let v_len2 = v_basis.dot(&v_basis);
    let a0 = (pa.dot(&v_basis) / v_len2.sqrt()).atan2(pa.dot(&u_basis) / u_len2.sqrt());
    let a1 = (pb.dot(&v_basis) / v_len2.sqrt()).atan2(pb.dot(&u_basis) / u_len2.sqrt());
    (a0, a1)
}

fn in_plane_basis(normal: Vec3) -> (Vec3, Vec3) {
    let n = normal.normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
    let world = if n.x.abs() < n.y.abs() && n.x.abs() < n.z.abs() {
        Vec3::new(1.0, 0.0, 0.0)
    } else if n.y.abs() < n.z.abs() {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let u = n
        .cross(&world)
        .normalized()
        .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
    let v = n.cross(&u).normalized().unwrap_or(Vec3::new(0.0, 1.0, 0.0));
    (u, v)
}
