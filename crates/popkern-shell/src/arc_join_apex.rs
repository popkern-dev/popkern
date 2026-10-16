//! Dilated cone-apex sphere caps (PK-326): the Minkowski dilation of a
//! cone rounds its apex with the sphere of radius d centered at the
//! ORIGINAL apex, tangent to the offset cone along an exact circle.
//! Surgery: split the (recomputed) seam edge at the tangency point,
//! `mef` a closed tangency-circle edge from the split vertex to itself
//! (separating the tip), `kev` the apex spur, and swap the tip face's
//! surface for the sphere. The tangency circle is exact: center
//! `A' + â·d·cos²β/sinβ`, radius `d·cosβ`, and every circle point
//! sits at distance d from the original apex (|C−A|² + r² = d²).

use popkern_math::PortableMath;
use std::collections::{BTreeSet, HashMap};

use crate::arc_join::cap_orientation;
use crate::skin_map::{JoinKind, JoinRecord, SkinMap};
use popkern_core::geometry::analytic_curves::circle_basis;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle, VertexHandle};

/// Construct sphere caps on every served dilated cone apex. Returns
/// the join records and the OUTER apex vertices that were served (the
/// caller retires their `ConeApexSharp` warnings). Apexes the surgery
/// cannot serve (tangency outside the lateral span, non-seam apex
/// valence) are skipped — their warnings stand as containment.
pub(crate) fn construct_cone_apex_caps(
    skin: &mut BRep,
    offset_surfaces: &HashMap<FaceHandle, SurfaceGeometry>,
    skin_map: &SkinMap,
    outer: &BRep,
    tolerance: &Tolerance,
) -> Result<(Vec<JoinRecord>, BTreeSet<VertexHandle>), String> {
    let reverse_vertex = skin_map.reverse_vertex_map();
    let mut records = Vec::new();
    let mut served: BTreeSet<VertexHandle> = BTreeSet::new();

    // Dilated offset cones, sorted for determinism.
    let mut cones: Vec<(FaceHandle, Point3, Vec3, f64, Point3)> = Vec::new();
    let mut cone_faces: Vec<(&FaceHandle, &SurfaceGeometry)> = offset_surfaces
        .iter()
        .filter(|(_, s)| matches!(s, SurfaceGeometry::Cone { .. }))
        .collect();
    cone_faces.sort_by_key(|(fh, _)| **fh);
    for (outer_fh, off) in cone_faces {
        let Some(SurfaceGeometry::Cone {
            apex: orig_apex,
            axis,
            half_angle,
        }) = outer
            .faces
            .get(*outer_fh)
            .and_then(|f| f.surface.as_ref())
            .map(|s| &s.geometry)
        else {
            continue;
        };
        let SurfaceGeometry::Cone {
            apex: offset_apex, ..
        } = off
        else {
            continue;
        };
        let shift = *offset_apex - *orig_apex;
        if shift.dot(axis) >= 0.0 {
            continue; // erosion or unmoved: sharp apex is exact
        }
        cones.push((*outer_fh, *orig_apex, *axis, *half_angle, *offset_apex));
    }

    for (_outer_fh, orig_apex, axis, half_angle, offset_apex) in cones {
        let Some(a_hat) = axis.normalized() else {
            continue;
        };
        let s = offset_apex.distance_to(&orig_apex);
        let d = s * half_angle.psin();
        if d <= tolerance.linear() * 10.0 {
            continue;
        }

        // Skin vertices whose OUTER original sits at the apex.
        let mut apex_vertices: Vec<VertexHandle> = skin
            .vertices
            .keys()
            .filter(|vh| {
                reverse_vertex.get(vh).is_some_and(|ov| {
                    outer.vertices.get(*ov).is_some_and(|v| {
                        v.position.distance_to(&orig_apex) <= tolerance.linear() * 10.0
                    })
                })
            })
            .collect();
        apex_vertices.sort();

        for vh in apex_vertices {
            if let Some(tip_face) = construct_one_apex_cap(
                skin,
                vh,
                &orig_apex,
                &offset_apex,
                a_hat,
                half_angle,
                d,
                tolerance,
            )? {
                let outer_vh = reverse_vertex.get(&vh).copied().unwrap_or(vh);
                served.insert(outer_vh);
                records.push(JoinRecord {
                    kind: JoinKind::SphereFill {
                        source_vertex: Some(outer_vh),
                    },
                    face: tip_face,
                });
            }
        }
    }
    Ok((records, served))
}

/// Serve one apex vertex. Returns Ok(None) when the apex stays sharp
/// (containment, warning remains), Ok(Some(tip_face)) on success.
#[allow(clippy::too_many_arguments)]
fn construct_one_apex_cap(
    skin: &mut BRep,
    vh: VertexHandle,
    orig_apex: &Point3,
    offset_apex: &Point3,
    a_hat: Vec3,
    half_angle: f64,
    d: f64,
    tolerance: &Tolerance,
) -> Result<Option<FaceHandle>, String> {
    // The apex must hang on exactly one seam edge (valence 1).
    let anchor = skin
        .vertices
        .get(vh)
        .map(|v| v.outgoing)
        .ok_or("apex vertex missing")?;
    let Some(anchor_he) = skin.half_edges.get(anchor) else {
        return Ok(None);
    };
    if anchor_he.origin != vh {
        return Ok(None);
    }
    let valence = skin.half_edges.values().filter(|h| h.origin == vh).count();
    if valence != 1 {
        return Ok(None); // trimmed / non-seam apex: not this class
    }
    let seam_eh = anchor_he.edge;
    let seam = skin.edges.get(seam_eh).ok_or("seam edge missing")?;
    let curve = seam
        .curve
        .clone()
        .ok_or("apex seam edge carries no curve")?;

    // Tangency point on the seam: slant distance d/tanβ from the
    // OFFSET apex.
    let slant = d / half_angle.ptan();
    let p0 = curve.geometry.point_at(curve.t_start);
    let p1 = curve.geometry.point_at(curve.t_end);
    let (apex_end_t, far_end_t, far_p) =
        if p0.distance_to(offset_apex) <= p1.distance_to(offset_apex) {
            (curve.t_start, curve.t_end, p1)
        } else {
            (curve.t_end, curve.t_start, p0)
        };
    let seam_len = far_p.distance_to(offset_apex);
    if !(seam_len.is_finite()) || slant >= seam_len * 0.95 || slant <= seam_len * 0.05 {
        return Ok(None); // tangency outside the lateral span: stays sharp
    }
    // Fraction along the curve parameter (seam slant lines are
    // uniform-speed; verify the evaluated point matches the exact
    // tangency position).
    let frac_from_apex = slant / seam_len;
    let t_split = apex_end_t + (far_end_t - apex_end_t) * frac_from_apex;
    let (t_lo, t_hi) = (
        curve.t_start.min(curve.t_end),
        curve.t_start.max(curve.t_end),
    );
    let frac = (t_split - t_lo) / (t_hi - t_lo);
    let dir = (far_p - *offset_apex) * (1.0 / seam_len);
    let t_exact = *offset_apex + dir * slant;
    let gate = (tolerance.linear() * 100.0).max(d * 1.0e-6);
    if curve.geometry.point_at(t_split).distance_to(&t_exact) > gate {
        return Ok(None); // non-uniform seam parametrization: stays sharp
    }

    let trace = popkern_core::trace::trace_enabled("shell-join");
    if trace {
        eprintln!("[apex] splitting seam {seam_eh:?} at frac {frac:.4}");
    }
    let split = skin
        .split_edge(seam_eh, frac)
        .map_err(|e| format!("apex seam split failed: {e}"))?;
    let v_t = split.vertex;
    skin.vertices[v_t].position = t_exact;

    // The two half-edges at v_t inside the lateral loop: one climbs to
    // the apex (its edge's far vertex is the apex), one descends.
    let mut he_up = None;
    let mut he_down = None;
    let mut cursor = skin.vertices[v_t].outgoing;
    for _ in 0..skin.half_edges.len() {
        let Some(he) = skin.half_edges.get(cursor) else {
            break;
        };
        if he.origin != v_t {
            break;
        }
        let Some(twin) = skin.half_edges.get(he.twin) else {
            break;
        };
        if twin.origin == vh {
            he_up = Some(cursor);
        } else {
            he_down = Some(cursor);
        }
        let next = twin.next;
        if next == skin.vertices[v_t].outgoing {
            break;
        }
        cursor = next;
    }
    let (Some(he_up), Some(he_down)) = (he_up, he_down) else {
        return Err("apex split fan did not yield up/down half-edges".into());
    };

    if trace {
        eprintln!("[apex] split ok, v_t={v_t:?}; mef {he_up:?} {he_down:?}");
    }
    let mef = skin
        .mef(he_up, he_down)
        .map_err(|e| format!("apex circle mef failed: {e}"))?;

    // Identify the tip face (the one whose loop reaches the apex).
    let tip_face = [mef.face, skin.loops[skin.half_edges[he_up].loop_].face]
        .into_iter()
        .find(|f| face_touches_vertex(skin, *f, vh))
        .ok_or("apex tip face not found after mef")?;

    if trace {
        eprintln!("[apex] mef ok, tip_face={tip_face:?}; repositioning pole");
    }
    // Keep the seam remnant as the cap's pole meridian (the shape the
    // sphere-cap machinery serves natively — a bare circle boundary
    // gives the mesher a zero-v-extent trim and the dome collapses to
    // a flat ring, measured on the first probe). Move the apex vertex
    // to the sphere's pole and author the meridian as a great arc.
    let pole = *orig_apex + a_hat * (-d);
    skin.vertices[vh].position = pole;
    let spur = loop_spur_into(skin, tip_face, vh).ok_or("apex spur not found in tip face")?;
    let seam_remnant = skin.half_edges[spur].edge;
    let meridian = crate::arc_join_geometry::great_arc(*orig_apex, d, t_exact, pole, tolerance)
        .ok_or("apex meridian arc is degenerate")?;
    skin.edges[seam_remnant].curve = Some(meridian);

    // Exact tangency circle, authored from v_t's angle so the closed
    // parameter range starts at the vertex.
    let center = *offset_apex + a_hat * (slant * half_angle.pcos());
    let radius = d * half_angle.pcos();
    let (u, v) = circle_basis(a_hat);
    let rvec = t_exact - center;
    let theta = rvec.dot(&v).patan2(rvec.dot(&u));
    let circle = CurveOnEdge {
        geometry: CurveGeometry::Arc {
            center,
            normal: a_hat,
            radius,
            start_angle: theta,
            end_angle: theta + std::f64::consts::TAU,
        },
        t_start: theta,
        t_end: theta + std::f64::consts::TAU,
    };
    // Self-checks: the circle passes through the split vertex and lies
    // on the apex sphere.
    if circle.geometry.point_at(theta).distance_to(&t_exact) > gate {
        return Err("apex tangency circle misses its seam vertex".into());
    }
    let probe = circle
        .geometry
        .point_at(theta + std::f64::consts::FRAC_PI_2);
    if (probe.distance_to(orig_apex) - d).abs() > gate {
        return Err("apex tangency circle leaves the apex sphere".into());
    }
    if trace {
        eprintln!("[apex] spur retired; circle authored");
    }
    skin.edges[mef.edge].curve = Some(circle);

    let sphere = SurfaceGeometry::Sphere {
        center: *orig_apex,
        radius: d,
        axis: popkern_core::geometry::Vec3::Z,
    };
    let orientation = cap_orientation(skin, tip_face, &sphere).unwrap_or(true);
    skin.faces[tip_face].surface = Some(SurfaceOnFace {
        geometry: sphere,
        orientation,
    });
    Ok(Some(tip_face))
}

fn face_touches_vertex(skin: &BRep, face: FaceHandle, vh: VertexHandle) -> bool {
    let Some(f) = skin.faces.get(face) else {
        return false;
    };
    let Some(loop_data) = skin.loops.get(f.outer_loop) else {
        return false;
    };
    let start = loop_data.half_edge;
    if !skin.half_edges.contains_key(start) {
        return false;
    }
    let mut cursor = start;
    for _ in 0..skin.half_edges.len() {
        let he = &skin.half_edges[cursor];
        if he.origin == vh {
            return true;
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    false
}

/// The spur half-edge in `face`'s outer loop whose destination is `vh`
/// (he.next == he.twin and twin.origin == vh).
fn loop_spur_into(
    skin: &BRep,
    face: FaceHandle,
    vh: VertexHandle,
) -> Option<popkern_core::topology::HalfEdgeHandle> {
    let start = skin.loops.get(skin.faces.get(face)?.outer_loop)?.half_edge;
    let mut cursor = start;
    for _ in 0..skin.half_edges.len() {
        let he = skin.half_edges.get(cursor)?;
        if he.next == he.twin && skin.half_edges.get(he.twin)?.origin == vh {
            return Some(cursor);
        }
        cursor = he.next;
        if cursor == start {
            break;
        }
    }
    None
}

/// Retire `ConeApexSharp` warnings whose vertex was served by a cap.
pub(crate) fn retire_served_apex_warnings(
    warnings: &mut Vec<crate::error::ShellWarning>,
    served: &BTreeSet<VertexHandle>,
) {
    warnings.retain(|w| {
        !matches!(w, crate::error::ShellWarning::ConeApexSharp { vertex } if served.contains(vertex))
    });
}
