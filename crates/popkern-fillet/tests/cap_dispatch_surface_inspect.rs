//! Path-1 Step B — inspect the post-rollback BRep for the cap_face
//! and support_face handles surfaced by the decline diagnostics on
//! hexagon6_variable, and dump their surfaces in detail. If the
//! SlotMap keys survive the rollback, we can read the support
//! plane's origin / u_axis / v_axis (and the cap surface kind /
//! CPs) without instrumenting production.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{
    ApexConstructionDiagnostic, BlendResult, FilletOptions, FilletRadius, FilletSpec, fillet_edges,
};

#[path = "kp_gate_helpers.rs"]
mod helpers;
use helpers::{incident_edges, pick_vertex_with_valence};

fn fillet_hexagon6_variable() -> BlendResult {
    let pyramid = BRep::n_gonal_pyramid(6, 1.0, 0.8, Tolerance::default()).unwrap();
    let apex = pick_vertex_with_valence(&pyramid, 6).expect("apex must exist");
    let edges = incident_edges(&pyramid, apex);
    let radii = [0.02, 0.04, 0.06, 0.06, 0.04, 0.02];
    let specs: Vec<_> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    fillet_edges(
        &pyramid,
        &specs,
        &FilletOptions::default(),
        Tolerance::default(),
        None,
    )
    .expect("fillet should not panic")
}

fn describe_surface(geometry: &SurfaceGeometry) -> String {
    match geometry {
        SurfaceGeometry::Plane {
            origin,
            normal,
            u_axis,
            v_axis,
        } => format!(
            "Plane {{ origin=({:.4},{:.4},{:.4}), normal=({:.4},{:.4},{:.4}), u_axis=({:.4},{:.4},{:.4}), v_axis=({:.4},{:.4},{:.4}) }}",
            origin.x,
            origin.y,
            origin.z,
            normal.x,
            normal.y,
            normal.z,
            u_axis.x,
            u_axis.y,
            u_axis.z,
            v_axis.x,
            v_axis.y,
            v_axis.z,
        ),
        SurfaceGeometry::Cylinder {
            origin,
            axis,
            radius,
        } => format!(
            "Cylinder {{ origin=({:.4},{:.4},{:.4}), axis=({:.4},{:.4},{:.4}), radius={:.4} }}",
            origin.x, origin.y, origin.z, axis.x, axis.y, axis.z, radius
        ),
        SurfaceGeometry::Sphere { center, radius } => format!(
            "Sphere {{ center=({:.4},{:.4},{:.4}), radius={:.4} }}",
            center.x, center.y, center.z, radius
        ),
        SurfaceGeometry::Nurbs(_) => "Nurbs(...)".to_string(),
        SurfaceGeometry::Gregory(g) => match g {
            GregoryPatch::RectangularBiNonic(b) => format!(
                "Gregory::RectangularBiNonic grid[0][0]=({:.4},{:.4},{:.4}) grid[9][9]=({:.4},{:.4},{:.4})",
                b.grid[0][0].x,
                b.grid[0][0].y,
                b.grid[0][0].z,
                b.grid[9][9].x,
                b.grid[9][9].y,
                b.grid[9][9].z,
            ),
            _ => format!("Gregory({:?})", std::any::type_name_of_val(g)),
        },
        _ => format!("{:?}", geometry),
    }
}

#[test]
#[ignore = "probe — run with --ignored --nocapture"]
fn hexagon6_variable_post_rollback_cap_support_surface_inspect() {
    let result = fillet_hexagon6_variable();
    eprintln!("\n## Post-rollback BRep surface inspect — faces named by the decline diagnostics\n");

    let mut targets: Vec<(
        popkern_core::topology::FaceHandle,
        popkern_core::topology::FaceHandle,
        popkern_core::topology::EdgeHandle,
        &str,
    )> = Vec::new();
    for d in &result.repair_report.apex_diagnostics {
        if let ApexConstructionDiagnostic::CapSupportContinuity(c) = d {
            targets.push((c.cap_face, c.support_face, c.edge, &c.measurement_scope));
        }
    }

    for (cap_face, support_face, edge, scope) in &targets {
        eprintln!("--- scope={} ---", scope);
        eprintln!("  cap_face={:?}", cap_face);
        match result.brep.faces.get(*cap_face) {
            Some(f) => match &f.surface {
                Some(s) => eprintln!("    surface = {}", describe_surface(&s.geometry)),
                None => eprintln!("    surface = <None>"),
            },
            None => eprintln!("    (face handle is STALE in post-rollback brep)"),
        }
        eprintln!("  support_face={:?}", support_face);
        match result.brep.faces.get(*support_face) {
            Some(f) => match &f.surface {
                Some(s) => eprintln!("    surface = {}", describe_surface(&s.geometry)),
                None => eprintln!("    surface = <None>"),
            },
            None => eprintln!("    (face handle is STALE in post-rollback brep)"),
        }
        eprintln!("  edge={:?}", edge);
        match result.brep.edges.get(*edge) {
            Some(e) => {
                eprintln!("    edge.curve.is_some = {}", e.curve.is_some());
                if let Some(curve) = &e.curve {
                    let p0 = curve.geometry.point_at(curve.t_start);
                    let p1 = curve.geometry.point_at(curve.t_end);
                    eprintln!(
                        "    curve start=({:.4},{:.4},{:.4}) end=({:.4},{:.4},{:.4})",
                        p0.x, p0.y, p0.z, p1.x, p1.y, p1.z
                    );
                }
                let he_a = e.half_edge;
                let he_a_data = result.brep.half_edges.get(he_a);
                let he_b = he_a_data.map(|d| d.twin);
                eprintln!(
                    "    he_a={:?} pcurve.is_some={} ; he_b={:?} pcurve.is_some={}",
                    he_a,
                    he_a_data.and_then(|d| d.pcurve.as_ref()).is_some(),
                    he_b,
                    he_b.and_then(|h| result.brep.half_edges.get(h))
                        .and_then(|d| d.pcurve.as_ref())
                        .is_some(),
                );
                if let Some(ha) = he_a_data
                    && let Some(pcurve) = &ha.pcurve
                {
                    let (lo, hi) = (pcurve.t_start, pcurve.t_end);
                    let uv0 = pcurve.point_at(lo);
                    let uv1 = pcurve.point_at(hi);
                    eprintln!(
                        "    he_a pcurve [{:.4},{:.4}] -> uv start=({:.4},{:.4}) end=({:.4},{:.4})",
                        lo, hi, uv0.x, uv0.y, uv1.x, uv1.y
                    );
                }
                if let Some(hb) = he_b.and_then(|h| result.brep.half_edges.get(h))
                    && let Some(pcurve) = &hb.pcurve
                {
                    let (lo, hi) = (pcurve.t_start, pcurve.t_end);
                    let uv0 = pcurve.point_at(lo);
                    let uv1 = pcurve.point_at(hi);
                    eprintln!(
                        "    he_b pcurve [{:.4},{:.4}] -> uv start=({:.4},{:.4}) end=({:.4},{:.4})",
                        lo, hi, uv0.x, uv0.y, uv1.x, uv1.y
                    );
                }
                if let (Some(ha), Some(hb_data)) =
                    (he_a_data, he_b.and_then(|h| result.brep.half_edges.get(h)))
                {
                    let va = result.brep.vertices.get(ha.origin).map(|v| v.position);
                    let vb = result.brep.vertices.get(hb_data.origin).map(|v| v.position);
                    eprintln!(
                        "    he_a origin vertex pos={:?} ; he_b origin vertex pos={:?}",
                        va, vb
                    );
                }
            }
            None => eprintln!("    (edge handle is STALE)"),
        }
    }
}

use popkern_core::geometry::curve::CurveEval;
