//! Exact-boundary support-face lifting for paired cap/support seams.

use super::apex_face_lift::{
    lift_apex_face_with_bridge, refit_face_outer_loop_pcurves,
    write_analytic_boundary_pcurve_on_support_side,
};
use super::apex_face_lift_boundary_curves::{
    compatible_unit_boundary_curves_for_support, write_edge_curve_from_boundary_row,
};
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};
use std::collections::BTreeMap;

pub(crate) fn lift_support_faces_for_bridge_edges(
    brep: &mut BRep,
    jobs: &[(FaceHandle, EdgeHandle)],
) -> usize {
    let trace = std::env::var_os("POPKERN_APEX_FACE_LIFT_TRACE").is_some();
    let mut jobs_by_face: BTreeMap<FaceHandle, Vec<EdgeHandle>> = BTreeMap::new();
    for &(support_face, bridge_edge) in jobs {
        let edges = jobs_by_face.entry(support_face).or_default();
        if !edges.contains(&bridge_edge) {
            edges.push(bridge_edge);
        }
    }

    let mut count = 0;
    for (support_face, bridge_edges) in jobs_by_face {
        if bridge_edges.len() >= 2 {
            let edge_a = bridge_edges[0];
            let edge_b = bridge_edges[1];
            match lift_support_face_with_boundary_pair(brep, support_face, edge_a, edge_b) {
                Ok(()) => {
                    write_analytic_boundary_pcurve_on_support_side(brep, support_face, edge_a, 0.0);
                    write_analytic_boundary_pcurve_on_support_side(brep, support_face, edge_b, 1.0);
                    if trace {
                        eprintln!(
                            "[apex-face-lift]   OK exact-boundary support={support_face:?} boundaries={edge_a:?},{edge_b:?}"
                        );
                    }
                    count += 1;
                    continue;
                }
                Err(msg) => {
                    if trace {
                        eprintln!(
                            "[apex-face-lift]   ERR exact-boundary pair support={support_face:?} boundaries={edge_a:?},{edge_b:?} : {msg}"
                        );
                    }
                }
            }
        }

        for bridge_edge in bridge_edges {
            match lift_apex_face_with_bridge(brep, support_face, bridge_edge) {
                Ok(()) => {
                    refit_face_outer_loop_pcurves(brep, support_face);
                    write_analytic_boundary_pcurve_on_support_side(
                        brep,
                        support_face,
                        bridge_edge,
                        0.0,
                    );
                    if trace {
                        eprintln!(
                            "[apex-face-lift]   OK exact-boundary support={support_face:?} bridge={bridge_edge:?}"
                        );
                    }
                    count += 1;
                }
                Err(msg) => {
                    if trace {
                        eprintln!(
                            "[apex-face-lift]   ERR exact-boundary support={support_face:?} bridge={bridge_edge:?} : {msg}"
                        );
                    }
                }
            }
        }
    }
    count
}

fn lift_support_face_with_boundary_pair(
    brep: &mut BRep,
    support_face: FaceHandle,
    edge_a: EdgeHandle,
    edge_b: EdgeHandle,
) -> Result<(), String> {
    let face_data = brep
        .faces
        .get(support_face)
        .ok_or_else(|| "support face not found".to_string())?;
    let sof = face_data
        .surface
        .as_ref()
        .ok_or_else(|| "support face has no surface".to_string())?;

    let (curve_a, curve_b) =
        compatible_unit_boundary_curves_for_support(brep, support_face, edge_a, edge_b)?;
    let num_v = curve_a.control_points.len();
    if num_v != curve_b.control_points.len()
        || curve_a.weights.len() != num_v
        || curve_b.weights.len() != num_v
        || curve_a.knots.len() != curve_b.knots.len()
    {
        return Err("boundary curves could not be made compatible".into());
    }

    let mut control_points = Vec::with_capacity(num_v * 2);
    control_points.extend(curve_a.control_points.iter().copied());
    control_points.extend(curve_b.control_points.iter().copied());

    let mut weights = Vec::with_capacity(num_v * 2);
    weights.extend(curve_a.weights.iter().copied());
    weights.extend(curve_b.weights.iter().copied());

    let surface = NurbsSurface {
        degree_u: 1,
        degree_v: curve_a.degree,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: curve_a.knots.clone(),
        control_points,
        num_u: 2,
        num_v,
        weights,
    };

    let mut new_sof = sof.clone();
    new_sof.geometry = SurfaceGeometry::Nurbs(surface);
    brep.faces[support_face].surface = Some(new_sof);
    write_edge_curve_from_boundary_row(brep, edge_a, curve_a);
    write_edge_curve_from_boundary_row(brep, edge_b, curve_b);

    Ok(())
}
