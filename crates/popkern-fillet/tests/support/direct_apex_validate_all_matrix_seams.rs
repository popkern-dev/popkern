use super::*;

#[test]
#[ignore = "native Gregory seam-ownership gate; current repeated-apex direct-apex path is under containment fallback"]
fn direct_apex_native_cap_support_seam_ownership_inventory() {
    struct Case {
        name: &'static str,
        radii: &'static [f64],
        is_pyramid: bool,
    }
    let cases: &[Case] = &[
        Case {
            name: "cube_mixed",
            radii: &[0.08, 0.10, 0.12],
            is_pyramid: false,
        },
        Case {
            name: "cube_extreme",
            radii: &[0.02, 0.10, 0.20],
            is_pyramid: false,
        },
        Case {
            name: "pyramid3_uniform",
            radii: &[0.10, 0.10, 0.10],
            is_pyramid: true,
        },
        Case {
            name: "pyramid4_mixed",
            radii: &[0.04, 0.06, 0.05, 0.07],
            is_pyramid: true,
        },
    ];
    for case in cases {
        let result = build_under_direct_apex(case.is_pyramid, case.radii);
        let cap_faces = result
            .corner_patches
            .iter()
            .map(|p| p.face)
            .collect::<std::collections::BTreeSet<_>>();
        let cap_methods = result
            .corner_patches
            .iter()
            .map(|p| (p.face, p.method.as_str()))
            .collect::<std::collections::BTreeMap<_, _>>();
        let mut findings = Vec::new();
        for (edge, edge_data) in result.brep.edges.iter() {
            let he = edge_data.half_edge;
            let Some(ha) = result.brep.half_edges.get(he) else {
                continue;
            };
            let Some(hb) = result.brep.half_edges.get(ha.twin) else {
                continue;
            };
            let Some(loop_a) = result.brep.loops.get(ha.loop_) else {
                continue;
            };
            let Some(loop_b) = result.brep.loops.get(hb.loop_) else {
                continue;
            };
            let face_a = loop_a.face;
            let face_b = loop_b.face;
            if face_a == face_b {
                continue;
            }
            let a_is_cap = cap_faces.contains(&face_a);
            let b_is_cap = cap_faces.contains(&face_b);
            if a_is_cap == b_is_cap {
                continue;
            }
            let cap_he = if a_is_cap { he } else { ha.twin };
            let support_he = if a_is_cap { ha.twin } else { he };
            let support_face = if a_is_cap { face_b } else { face_a };
            let support_kind = result
                .brep
                .faces
                .get(support_face)
                .and_then(|f| f.surface.as_ref())
                .map(|s| match &s.geometry {
                    popkern_core::geometry::surface::SurfaceGeometry::Plane { .. } => "Plane",
                    popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. } => "Cylinder",
                    popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. } => "Sphere",
                    popkern_core::geometry::surface::SurfaceGeometry::Torus { .. } => "Torus",
                    popkern_core::geometry::surface::SurfaceGeometry::Cone { .. } => "Cone",
                    popkern_core::geometry::surface::SurfaceGeometry::Gregory(_) => "Gregory",
                    popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_) => "Nurbs",
                })
                .unwrap_or("None");
            let issue = match &edge_data.curve {
                None => Some("missing 3D curve"),
                Some(curve)
                    if matches!(curve.geometry, CurveGeometry::Nurbs(_))
                        && curve
                            .geometry
                            .as_line_segment_in_range((curve.t_start, curve.t_end), 1.0e-8)
                            .is_some() =>
                {
                    Some("straight NURBS")
                }
                Some(_) => None,
            };
            let missing_cap_pcurve = result.brep.half_edges[cap_he].pcurve.is_none();
            let missing_support_pcurve = result.brep.half_edges[support_he].pcurve.is_none();
            if issue.is_some() || missing_cap_pcurve || missing_support_pcurve {
                let start = result.brep.half_edges[cap_he].origin;
                let end = result.brep.half_edges[result.brep.half_edges[cap_he].twin].origin;
                let cap_face = if a_is_cap { face_a } else { face_b };
                let on_cap_outer = result
                    .brep
                    .faces
                    .get(cap_face)
                    .and_then(|f| result.brep.loops.get(f.outer_loop))
                    .is_some_and(|loop_data| {
                        let mut cursor = loop_data.half_edge;
                        for _ in 0..(result.brep.half_edges.len() + 1) {
                            if cursor == cap_he {
                                return true;
                            }
                            let Some(he_data) = result.brep.half_edges.get(cursor) else {
                                return false;
                            };
                            cursor = he_data.next;
                            if cursor == loop_data.half_edge {
                                return false;
                            }
                        }
                        false
                    });
                findings.push(format!(
                    "{edge:?} he_edge={:?} cap={cap_face:?}/{} support={support_face:?}/{support_kind} cap_he={cap_he:?} outer={on_cap_outer} support_he={support_he:?} verts={start:?}->{end:?} issue={} cap_pcurve={} support_pcurve={}",
                    result.brep.half_edges[he].edge,
                    cap_methods
                        .get(&cap_face)
                        .copied()
                        .unwrap_or("unknown"),
                    issue.unwrap_or("none"),
                    !missing_cap_pcurve,
                    !missing_support_pcurve
                ));
            }
        }
        eprintln!(
            "[direct-apex seam ownership] {}: {} finding(s){}",
            case.name,
            findings.len(),
            if findings.is_empty() {
                String::new()
            } else {
                format!(": {}", findings.join(", "))
            }
        );
        if !findings.is_empty() && std::env::var_os("POPKERN_DIRECT_APEX_TRACE").is_some() {
            dump_cap_loops(case.name, &result, &cap_faces);
        }
        assert!(
            findings.is_empty(),
            "[{}] direct-apex cap/support seam ownership findings: {}",
            case.name,
            findings.join(", ")
        );
    }
}

fn dump_cap_loops(
    case_name: &str,
    result: &popkern_fillet::BlendResult,
    cap_faces: &std::collections::BTreeSet<popkern_core::topology::FaceHandle>,
) {
    for cap_face in cap_faces {
        let Some(face_data) = result.brep.faces.get(*cap_face) else {
            continue;
        };
        let Some(loop_data) = result.brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let mut cursor = loop_data.half_edge;
        let mut entries = Vec::new();
        for _ in 0..(result.brep.half_edges.len() + 1) {
            let Some(he) = result.brep.half_edges.get(cursor) else {
                break;
            };
            entries.push(format!("{:?}:{:?}", cursor, he.edge));
            cursor = he.next;
            if cursor == loop_data.half_edge {
                break;
            }
        }
        eprintln!(
            "[direct-apex seam ownership loop] {} {cap_face:?}: {}",
            case_name,
            entries.join(" -> ")
        );
    }
}
