//! Geometry edits for blend-boundary seam re-anchoring: vertex
//! snap with incident-curve re-trim, seam edge curve rewrite, and
//! support-side pcurve refit. See the parent module docs for the
//! transactional model.

use super::*;

/// Move `v` to `target`, keeping every incident edge curve joined:
/// NURBS curves get an endpoint-CP drag; analytic curves a
/// closest-point re-trim (their endpoint must already lie on the
/// curve to within `tolerance` — the junction slides ALONG a wall
/// line or fillet arc, it never leaves it). Edges in `seam_edges`
/// are exempt from the residual check — they are rewritten with the
/// exact segment curve immediately after.
pub(super) fn snap_vertex_with_retrim(
    brep: &mut BRep,
    v: VertexHandle,
    target: Point3,
    tolerance: f64,
    seam_edges: &[EdgeHandle],
) -> Result<(), &'static str> {
    let Some(vd) = brep.vertices.get(v) else {
        return Err("snap_vertex_missing");
    };
    if vd.position.distance_to(&target) < 1e-12 {
        return Ok(());
    }
    let mut off_curve_probe: Option<EdgeHandle> = None;
    let halves: Vec<HalfEdgeHandle> = brep
        .half_edges
        .iter()
        .filter_map(|(h, he)| (he.origin == v).then_some(h))
        .collect();
    // PK-222: straight support-face chords cannot be re-trimmed onto
    // a junction they never passed through — their endpoint genuinely
    // moves with the ring vertex. When both carrier faces are planes
    // that contain the target, rebuild the chord between the target
    // and its far endpoint (the post-contraction rebuild contract);
    // curved carriers still take the residual-gated re-trim below.
    let mut rebuilt: Vec<EdgeHandle> = Vec::new();
    for &he_h in &halves {
        let edge_h = brep.half_edges[he_h].edge;
        if seam_edges.contains(&edge_h) || rebuilt.contains(&edge_h) {
            continue;
        }
        let Some(edge_data) = brep.edges.get(edge_h) else {
            continue;
        };
        let Some(coe) = edge_data.curve.as_ref() else {
            continue;
        };
        if !matches!(coe.geometry, CurveGeometry::Line { .. }) {
            continue;
        }
        let residual_ok = coe
            .geometry
            .closest_point(&target)
            .map(|t| coe.geometry.point_at(t).distance_to(&target) <= tolerance)
            .unwrap_or(false);
        if residual_ok {
            continue;
        }
        if !edge_side_faces_contain_target(brep, edge_h, &target, tolerance) {
            continue;
        }
        let he0 = edge_data.half_edge;
        let Some(primary) = brep.half_edges.get(he0) else {
            continue;
        };
        let Some(twin) = brep.half_edges.get(primary.twin) else {
            continue;
        };
        let (origin_v, dest_v) = (primary.origin, twin.origin);
        let other = if origin_v == v { dest_v } else { origin_v };
        let Some(other_pos) = brep.vertices.get(other).map(|vd| vd.position) else {
            continue;
        };
        let (a, b) = if origin_v == v {
            (target, other_pos)
        } else {
            (other_pos, target)
        };
        let d = b - a;
        let len = d.length();
        let Some(dir) = (len > 1e-12).then(|| d.normalized()).flatten() else {
            continue;
        };
        let (he0, twin_h) = (he0, primary.twin);
        if let Some(edge_data) = brep.edges.get_mut(edge_h) {
            edge_data.curve = Some(CurveOnEdge {
                geometry: CurveGeometry::Line {
                    origin: a,
                    direction: dir,
                },
                t_start: 0.0,
                t_end: len,
            });
        }
        // Stale pcurves would resurrect the old chord; downstream
        // reprojection re-authors from the rebuilt line.
        if let Some(he) = brep.half_edges.get_mut(he0) {
            he.pcurve = None;
        }
        if let Some(he) = brep.half_edges.get_mut(twin_h) {
            he.pcurve = None;
        }
        rebuilt.push(edge_h);
    }
    for he_h in halves {
        let edge_h = brep.half_edges[he_h].edge;
        if rebuilt.contains(&edge_h) {
            continue;
        }
        let is_primary = brep.edges.get(edge_h).map(|e| e.half_edge) == Some(he_h);
        let is_seam = seam_edges.contains(&edge_h);
        let side_kinds = if crate::direct_apex_trace::cap_g0_enabled() {
            edge_side_face_kinds(brep, edge_h)
        } else {
            String::new()
        };
        let Some(coe) = brep.edges.get_mut(edge_h).and_then(|e| e.curve.as_mut()) else {
            continue;
        };
        match &mut coe.geometry {
            CurveGeometry::Nurbs(nc) if !nc.control_points.is_empty() => {
                let n = nc.control_points.len();
                if is_primary {
                    nc.control_points[0] = target;
                } else {
                    nc.control_points[n - 1] = target;
                }
            }
            geometry => {
                // Analytic curve: slide the trim end to the closest
                // parameter. The junction must lie on the curve.
                let Some(t) = geometry.closest_point(&target) else {
                    if is_seam {
                        continue;
                    }
                    return Err("incident_retrim_projection_failed");
                };
                let residual = geometry.point_at(t).distance_to(&target);
                if residual > tolerance && !is_seam {
                    if crate::direct_apex_trace::cap_g0_enabled() {
                        eprintln!(
                            "[direct-apex-trace cap-g0] reanchor retrim off-curve: edge={edge_h:?} \
                             kind={} residual={residual:.3e} tol={tolerance:.3e} v={v:?} sides={side_kinds}",
                            curve_kind_name(&coe.geometry),
                        );
                        off_curve_probe = Some(edge_h);
                        continue;
                    }
                    return Err("incident_retrim_off_curve");
                }
                if is_primary {
                    coe.t_start = t;
                } else {
                    coe.t_end = t;
                }
            }
        }
    }
    if off_curve_probe.is_some() {
        return Err("incident_retrim_off_curve");
    }
    if let Some(vd) = brep.vertices.get_mut(v) {
        vd.position = target;
    }
    Ok(())
}

/// Write the exact boundary segment curve into the seam edge
/// (oriented to the edge's primary half-edge) and refit the
/// support-side pcurve on the twin face. The cap-side pcurve is
/// left for `reproject_loop_pcurves_onto_gregory`, which rebuilds
/// every cap loop pcurve from edge curves after surface assignment.
pub(super) fn rewire_seam_edge(
    brep: &mut BRep,
    boundary: &BlendBoundary,
    job: &ReanchorJob,
) -> Result<(), &'static str> {
    let segment = segment_curve(&boundary.segments[job.segment_index]);
    let our_is_primary = brep.edges.get(job.edge).map(|e| e.half_edge) == Some(job.half_edge);
    // `job.reversed` relates the segment direction to OUR half-edge;
    // the stored curve must run with the PRIMARY half-edge. Reverse
    // when the segment opposes the primary: segment-vs-ours XOR
    // ours-vs-primary, which reduces to the equality below.
    let curve_reversed = job.reversed == our_is_primary;
    let oriented = if curve_reversed {
        reverse_nurbs(segment)
    } else {
        segment.clone()
    };
    let (t0, t1) = oriented.parameter_range();
    let geometry = CurveGeometry::Nurbs(oriented.clone());
    let geometry = geometry
        .canonical_line_in_range((t0, t1), 1.0e-8)
        .unwrap_or(geometry);
    let Some(edge_data) = brep.edges.get_mut(job.edge) else {
        return Err("rewire_edge_missing");
    };
    edge_data.curve = Some(CurveOnEdge {
        geometry,
        t_start: t0,
        t_end: t1,
    });

    // Support-side pcurve refit.
    let twin = brep.half_edges[job.half_edge].twin;
    let twin_surface = brep
        .half_edges
        .get(twin)
        .and_then(|he| brep.loops.get(he.loop_))
        .and_then(|l| brep.faces.get(l.face))
        .and_then(|f| f.surface.as_ref())
        .map(|sof| sof.geometry.clone());
    let Some(twin_surface) = twin_surface else {
        return Err("rewire_twin_surface_missing");
    };
    let pcurve = support_pcurve_for_curve(&oriented, &twin_surface, t0, t1)
        .ok_or("rewire_support_pcurve_fit_failed")?;
    if let Some(he) = brep.half_edges.get_mut(twin) {
        he.pcurve = Some(pcurve);
    }
    Ok(())
}

/// True when both side faces of `edge` carry surfaces that contain
/// `target` (the straight-chord rebuild precondition — a rebuilt
/// chord's endpoints must stay on both carriers). Planes are gated
/// at the authoring tolerance; fitted carriers (strips, caps) at the
/// carve-chord class (1e-3) — the rebuilt chord is exactly as
/// approximate against them as the carve chord it replaces.
fn edge_side_faces_contain_target(
    brep: &BRep,
    edge: EdgeHandle,
    target: &Point3,
    gate: f64,
) -> bool {
    use popkern_core::geometry::surface::SurfaceEval;
    let Some(edge_data) = brep.edges.get(edge) else {
        return false;
    };
    let Some(primary) = brep.half_edges.get(edge_data.half_edge) else {
        return false;
    };
    for he in [edge_data.half_edge, primary.twin] {
        let contains = brep
            .half_edges
            .get(he)
            .and_then(|h| brep.loops.get(h.loop_))
            .and_then(|l| brep.faces.get(l.face))
            .and_then(|f| f.surface.as_ref())
            .is_some_and(|sof| {
                let side_gate = if matches!(sof.geometry, SurfaceGeometry::Plane { .. }) {
                    gate
                } else {
                    gate.max(1.0e-3)
                };
                let distance = sof
                    .geometry
                    .closest_point(target)
                    .map(|(u, v)| sof.geometry.point_at(u, v).distance_to(target));
                if crate::direct_apex_trace::cap_g0_enabled()
                    && distance.is_none_or(|d| d > side_gate)
                {
                    eprintln!(
                        "[direct-apex-trace cap-g0] chord-rebuild gate: edge={edge:?} side {} \
                         off target by {distance:?} (gate {side_gate:.1e})",
                        surface_kind_name(&sof.geometry)
                    );
                }
                distance.is_some_and(|d| d <= side_gate)
            });
        if !contains {
            return false;
        }
    }
    true
}

/// Analytic Linear pcurve on planes (exact), otherwise a sampled
/// degree-3 fit through UV projections.
fn support_pcurve_for_curve(
    curve: &NurbsCurve,
    surface: &SurfaceGeometry,
    t0: f64,
    t1: f64,
) -> Option<Pcurve> {
    if let SurfaceGeometry::Plane { .. } = surface {
        let uv0 = project_point_to_uv(surface, &curve.point_at(t0))?;
        let uv1 = project_point_to_uv(surface, &curve.point_at(t1))?;
        return Some(Pcurve::linear(uv0, uv1, t0, t1));
    }
    const N: usize = 12;
    let mut points = Vec::with_capacity(N);
    let mut params = Vec::with_capacity(N);
    for i in 0..N {
        let s = i as f64 / (N - 1) as f64;
        let t = t0 + s * (t1 - t0);
        let (u, v) = project_point_to_uv(surface, &curve.point_at(t))?;
        points.push(Point3::new(u, v, 0.0));
        params.push(t);
    }
    fitting::fit_curve(&points, &params, 3, 5.min(points.len()))
        .ok()
        .map(Into::into)
}

pub(super) fn edge_curve_endpoint_gap(brep: &BRep, job: &ReanchorJob) -> Option<f64> {
    let edge = brep.edges.get(job.edge)?;
    let coe = edge.curve.as_ref()?;
    let primary = brep.half_edges.get(edge.half_edge)?;
    let p_start = brep.vertices.get(primary.origin)?.position;
    let p_end = brep
        .vertices
        .get(brep.half_edges.get(primary.twin)?.origin)?
        .position;
    let c_start = coe.geometry.point_at(coe.t_start);
    let c_end = coe.geometry.point_at(coe.t_end);
    Some(p_start.distance_to(&c_start).max(p_end.distance_to(&c_end)))
}

/// Reverse a NURBS curve (CPs, weights, reflected knots).
fn reverse_nurbs(curve: &NurbsCurve) -> NurbsCurve {
    let (t_min, t_max) = curve.parameter_range();
    let mut knots: Vec<f64> = curve
        .knots
        .iter()
        .rev()
        .map(|&k| t_min + t_max - k)
        .collect();
    knots.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let mut control_points = curve.control_points.clone();
    control_points.reverse();
    let mut weights = curve.weights.clone();
    weights.reverse();
    NurbsCurve {
        degree: curve.degree,
        knots,
        control_points,
        weights,
    }
}

fn edge_side_face_kinds(brep: &BRep, edge: EdgeHandle) -> String {
    let Some(edge_data) = brep.edges.get(edge) else {
        return "?".into();
    };
    let Some(primary) = brep.half_edges.get(edge_data.half_edge) else {
        return "?".into();
    };
    [edge_data.half_edge, primary.twin]
        .into_iter()
        .map(|he| {
            brep.half_edges
                .get(he)
                .and_then(|h| brep.loops.get(h.loop_))
                .and_then(|l| brep.faces.get(l.face))
                .map(|f| {
                    let kind = f
                        .surface
                        .as_ref()
                        .map(|sof| surface_kind_name(&sof.geometry))
                        .unwrap_or("none");
                    kind.to_string()
                })
                .unwrap_or_else(|| "dead".into())
        })
        .collect::<Vec<_>>()
        .join("|")
}

fn surface_kind_name(geometry: &SurfaceGeometry) -> &'static str {
    match geometry {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cylinder",
        SurfaceGeometry::Sphere { .. } => "Sphere",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
        _ => "Other",
    }
}

fn curve_kind_name(geometry: &CurveGeometry) -> &'static str {
    match geometry {
        CurveGeometry::Line { .. } => "Line",
        CurveGeometry::Circle { .. } => "Circle",
        CurveGeometry::Nurbs(_) => "Nurbs",
        _ => "Other",
    }
}
