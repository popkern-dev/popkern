//! Edge analysis for fillet/chamfer operations.
//!
//! For each edge to fillet or chamfer, identifies the two adjacent faces,
//! classifies edge convexity, detects G1-continuous chains, and groups
//! vertices where multiple fillet edges meet.

mod chains;

pub use chains::{detect_chains, group_vertices};

use popkern_core::geometry::Vec3;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle, VertexHandle};

/// Analysis of a single edge for filleting.
#[derive(Debug)]
pub struct EdgeAnalysis {
    /// The edge handle in the BRep.
    pub edge: EdgeHandle,
    /// First adjacent face.
    pub face_1: FaceHandle,
    /// Second adjacent face.
    pub face_2: FaceHandle,
    /// Surface geometry of the first face.
    pub surface_1: SurfaceGeometry,
    /// Surface geometry of the second face.
    pub surface_2: SurfaceGeometry,
    /// Face orientation for surface 1. When `false`, the effective face
    /// normal is opposite to the surface geometry's normal.
    pub orientation_1: bool,
    /// Face orientation for surface 2.
    pub orientation_2: bool,
    /// Whether the shell containing these faces is globally
    /// OUTWARD-oriented (measured via signed volume). Primitives
    /// build inward shells; boolean/sweep output builds outward
    /// ones. True face outwardness = `orientation_i == shell_outward`.
    pub shell_outward: bool,
    /// Convexity classification at the edge.
    pub convexity: Convexity,
    /// Vertex at the start of the edge.
    pub vertex_start: VertexHandle,
    /// Vertex at the end of the edge.
    pub vertex_end: VertexHandle,
}

/// Convexity classification of an edge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Convexity {
    /// The edge is convex (exterior angle < 180 degrees).
    Convex,
    /// The edge is concave (exterior angle > 180 degrees).
    Concave,
    /// Mixed or indeterminate convexity (e.g., tangent surfaces).
    Mixed,
}

/// Analyze edges for filleting.
///
/// For each edge, finds the two adjacent faces, retrieves their surface
/// geometry, classifies convexity, and identifies the start/end vertices.
/// If `propagate_tangent` is true, additional tangent-continuous edges
/// may be included (future extension point).
pub fn analyze_edges(
    brep: &BRep,
    edges: &[EdgeHandle],
    _propagate_tangent: bool,
) -> Result<Vec<EdgeAnalysis>, String> {
    let mut analyses = Vec::with_capacity(edges.len());

    for &edge_handle in edges {
        let analysis = analyze_single_edge(brep, edge_handle)?;
        analyses.push(analysis);
    }

    Ok(analyses)
}

/// Analyze a single edge: find adjacent faces, classify convexity, get vertices.
fn analyze_single_edge(brep: &BRep, edge_handle: EdgeHandle) -> Result<EdgeAnalysis, String> {
    let edge = brep
        .edges
        .get(edge_handle)
        .ok_or_else(|| format!("edge {edge_handle:?} not found in BRep"))?;

    let he1_handle = edge.half_edge;
    let he1 = brep
        .half_edges
        .get(he1_handle)
        .ok_or_else(|| format!("half-edge {he1_handle:?} not found"))?;

    let he2_handle = he1.twin;
    let he2 = brep
        .half_edges
        .get(he2_handle)
        .ok_or_else(|| format!("twin half-edge {he2_handle:?} not found"))?;

    // Walk half-edge -> loop -> face.
    let face_1 = loop_to_face(brep, he1.loop_)?;
    let face_2 = loop_to_face(brep, he2.loop_)?;

    // Seam edges have the same face on both sides — cannot be filleted.
    if face_1 == face_2 {
        return Err(format!(
            "edge {edge_handle:?} is a seam edge (same face on both sides) and cannot be filleted"
        ));
    }

    let (surface_1, stored_orientation_1) = get_face_surface_with_orientation(brep, face_1)?;
    let (surface_2, stored_orientation_2) = get_face_surface_with_orientation(brep, face_2)?;
    // EFFECTIVE orientation: does the surface's natural normal point
    // OUT of the material? Derived from the outer loop's UV winding
    // (outer loops wind CCW in the chart iff the natural normal is
    // outward — the manifold convention), not from the stored flag:
    // face producers disagree on the flag (fillet strips carry
    // radial-out cylinder normals into concave grooves; sweep caps
    // and boolean output vary), and the rolling-ball side depends on
    // the PHYSICAL outward direction. Falls back to the stored flag
    // when pcurves are unavailable.
    // orientation_1/2 keep legacy STORED-flag semantics (a large
    // consumer surface depends on them); shell_outward carries the
    // measured global shell orientation so ball placement can
    // compute true outwardness (stored == shell_outward).
    let shell_outward = shell_outward_sign(brep) > 0.0;
    let orientation_1 = stored_orientation_1;
    let orientation_2 = stored_orientation_2;

    // Get the two vertices at edge endpoints.
    let vertex_start = he1.origin;
    let vertex_end = he2.origin;

    // Classify convexity at the edge midpoint.
    let convexity = classify_convexity(brep, edge_handle, &surface_1, &surface_2, face_1, face_2)?;

    Ok(EdgeAnalysis {
        edge: edge_handle,
        face_1,
        face_2,
        surface_1,
        surface_2,
        orientation_1,
        orientation_2,
        shell_outward,
        convexity,
        vertex_start,
        vertex_end,
    })
}

/// Find the face that owns a given loop.
fn loop_to_face(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Result<FaceHandle, String> {
    let loop_ = brep
        .loops
        .get(loop_handle)
        .ok_or_else(|| format!("loop {loop_handle:?} not found"))?;
    Ok(loop_.face)
}

/// Get surface geometry for a face.
#[allow(dead_code)]
fn get_face_surface(brep: &BRep, face_handle: FaceHandle) -> Result<SurfaceGeometry, String> {
    let face = brep
        .faces
        .get(face_handle)
        .ok_or_else(|| format!("face {face_handle:?} not found"))?;

    let surface_on_face = face
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face_handle:?} has no surface geometry"))?;

    Ok(surface_on_face.geometry.clone())
}

/// Get the surface geometry AND orientation flag from a face.
fn get_face_surface_with_orientation(
    brep: &BRep,
    face_handle: FaceHandle,
) -> Result<(SurfaceGeometry, bool), String> {
    let face = brep
        .faces
        .get(face_handle)
        .ok_or_else(|| format!("face {face_handle:?} not found"))?;
    let surface_on_face = face
        .surface
        .as_ref()
        .ok_or_else(|| format!("face {face_handle:?} has no surface geometry"))?;
    Ok((
        surface_on_face.geometry.clone(),
        surface_on_face.orientation,
    ))
}

/// Number of sample points along the edge for convexity classification.
const CONVEXITY_SAMPLES: usize = 7;

/// Classify the convexity of an edge between two surfaces.
///
/// Samples normals at [`CONVEXITY_SAMPLES`] points along the edge and
/// uses the sign of (N1 x T) . N2 at each sample:
///   - positive -> convex
///   - negative -> concave
///   - near zero -> indeterminate (skipped)
///
/// If all definite samples agree, returns [`Convexity::Convex`] or
/// [`Convexity::Concave`]. If samples disagree, returns [`Convexity::Mixed`].
fn classify_convexity(
    brep: &BRep,
    edge_handle: EdgeHandle,
    surface_1: &SurfaceGeometry,
    surface_2: &SurfaceGeometry,
    face_1: FaceHandle,
    face_2: FaceHandle,
) -> Result<Convexity, String> {
    let edge = brep
        .edges
        .get(edge_handle)
        .ok_or_else(|| format!("edge {edge_handle:?} not found"))?;

    let curve = edge
        .curve
        .as_ref()
        .ok_or_else(|| format!("edge {edge_handle:?} has no curve geometry"))?;

    const MIXED_THRESHOLD: f64 = 1e-6;

    let mut convex_count: usize = 0;
    let mut concave_count: usize = 0;

    // Winding-free convexity: OUTWARD normals (natural normal
    // corrected by the winding-derived effective orientation) plus
    // INTO-FACE tangent directions (toward each face's uv-interior).
    // For a convex edge the material wedge is < π and each face's
    // interior direction opposes the OTHER face's outward normal;
    // for a concave (reentrant) edge it agrees:
    //     dot(d1, n2) + dot(d2, n1)  < 0  ⇔ convex, > 0 ⇔ concave.
    // (Box top edge: d1=+y into the top face, n2=−y outward front
    // ⇒ −1. Pocket floor/wall: d1=+y, n2=+y ⇒ +1.) The previous
    // formula used the RAW curve tangent — untied to either face's
    // winding — and stored orientation flags, so its sign was
    // producer-dependent: every convex plate edge on boolean output
    // classified Concave, sending the rolling ball to the exterior
    // weld-bead mirror.
    // Stored orientation flags are the manifold truth: primitives
    // encode outwardness in the flag alone (their loop windings are
    // unenforced), and boolean/sweep output is winding-flag
    // consistent since the PK-099 source fixes (validated: zero
    // mismatches across the flagship union). The former shoelace
    // override compensated for a since-fixed producer defect.
    let stored = |face: FaceHandle| -> bool {
        brep.faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
            .map(|sof| sof.orientation)
            .unwrap_or(true)
    };
    // GLOBAL shell orientation: primitive constructors build
    // INWARD-oriented shells (extrude keeps the base face's
    // winding), sweep/boolean output builds OUTWARD ones — the flag
    // is only honest relative to its shell. One signed-volume
    // measurement disambiguates.
    let shell = shell_outward_sign(brep);

    // Into-face direction at a boundary sample: walking a face's own
    // half-edge with the face's OUTWARD normal up, the interior is on
    // the LEFT — d = n_outward × t̂, with t̂ the half-edge traversal
    // direction from its pcurve. Holds for outer AND inner loops
    // (that is why inner loops wind opposite), unlike any
    // uv-box-center heuristic (the plate top's box center sits in
    // its own hole).
    let into_face_dir = |face: FaceHandle,
                         point: &popkern_core::geometry::Point3,
                         outward: &Vec3|
     -> Option<Vec3> {
        let face_data = brep.faces.get(face)?;
        let surface = &face_data.surface.as_ref()?.geometry;
        // The half-edge of THIS edge whose loop belongs to `face`.
        let edge_data = brep.edges.get(edge_handle)?;
        let he_a = edge_data.half_edge;
        let he_b = brep.half_edges.get(he_a)?.twin;
        let mut he_face = None;
        for cand in [he_a, he_b] {
            let d = brep.half_edges.get(cand)?;
            if brep.loops.get(d.loop_)?.face == face {
                he_face = Some(cand);
                break;
            }
        }
        let d = brep.half_edges.get(he_face?)?;
        let pcurve = d.pcurve.as_ref()?;
        let (t0, t1) = pcurve.parameter_range();
        // Pcurve parameter nearest the sample point (in 3D).
        let mut best = (f64::INFINITY, t0);
        for k in 0..=16 {
            let t = t0 + (t1 - t0) * k as f64 / 16.0;
            let uv = pcurve.point_at(t);
            let q = surface.point_at(uv.x, uv.y);
            let dist = q.distance_to(point);
            if dist < best.0 {
                best = (dist, t);
            }
        }
        let t = best.1;
        let h = (t1 - t0) * 1.0e-4;
        let (ta, tb) = if t + h <= t1 { (t, t + h) } else { (t - h, t) };
        let uv_a = pcurve.point_at(ta);
        let uv_b = pcurve.point_at(tb);
        let pa = surface.point_at(uv_a.x, uv_a.y);
        let pb = surface.point_at(uv_b.x, uv_b.y);
        let mut tangent = (pb - pa).normalized()?;
        // Direction verification (no repo-wide pcurve direction
        // contract; see effective_outward_orientation).
        let origin = brep.vertices.get(d.origin)?.position;
        let dest = brep
            .vertices
            .get(brep.half_edges.get(d.twin)?.origin)?
            .position;
        let uv0 = pcurve.point_at(t0);
        let q0 = surface.point_at(uv0.x, uv0.y);
        if q0.distance_to(&dest) + 1.0e-12 < q0.distance_to(&origin) {
            tangent = tangent * -1.0;
        }
        (outward.cross(&tangent)).normalized()
    };

    for i in 0..CONVEXITY_SAMPLES {
        let frac = (i as f64) / ((CONVEXITY_SAMPLES - 1).max(1) as f64);
        let t = curve.t_start + (curve.t_end - curve.t_start) * frac;
        let point = curve.geometry.point_at(t);

        // Outward normals fold in the shell sign; the interior-left
        // rule below is a property of the face's OWN manifold
        // normal (flag-corrected natural, shell-independent) — on an
        // inward-oriented shell the two differ and mixing them
        // mirrors the into-face direction.
        let nat_1 = surface_normal_at_point_natural(surface_1, &point);
        let nat_2 = surface_normal_at_point_natural(surface_2, &point);
        let flag_sign = |f: bool| if f { 1.0 } else { -1.0 };
        let m1 = nat_1 * flag_sign(stored(face_1));
        let m2 = nat_2 * flag_sign(stored(face_2));
        let n1 = m1 * shell;
        let n2 = m2 * shell;

        let Some(d1) = into_face_dir(face_1, &point, &m1) else {
            continue;
        };
        let Some(d2) = into_face_dir(face_2, &point, &m2) else {
            continue;
        };
        let dot = d1.dot(&n2) + d2.dot(&n1);
        if std::env::var_os("PK_CONVEXITY_TRACE").is_some() && i == 0 {
            eprintln!(
                "  [cvx] {edge_handle:?} p({:.2},{:.2},{:.2}) n1({:.2},{:.2},{:.2}) d1({:.2},{:.2},{:.2}) n2({:.2},{:.2},{:.2}) d2({:.2},{:.2},{:.2}) dot {dot:.3}",
                point.x,
                point.y,
                point.z,
                n1.x,
                n1.y,
                n1.z,
                d1.x,
                d1.y,
                d1.z,
                n2.x,
                n2.y,
                n2.z,
                d2.x,
                d2.y,
                d2.z
            );
        }

        if dot < -MIXED_THRESHOLD {
            convex_count += 1;
        } else if dot > MIXED_THRESHOLD {
            concave_count += 1;
        }
        // Near-zero samples are indeterminate — skipped.
    }

    if convex_count > 0 && concave_count > 0 {
        Ok(Convexity::Mixed)
    } else if convex_count > 0 {
        Ok(Convexity::Convex)
    } else if concave_count > 0 {
        Ok(Convexity::Concave)
    } else {
        // All samples indeterminate.
        Ok(Convexity::Mixed)
    }
}

/// Evaluate the surface normal at a 3D point, accounting for face orientation.
fn surface_normal_at_point_natural(
    surface: &SurfaceGeometry,
    point: &popkern_core::geometry::Point3,
) -> Vec3 {
    let (u, v) = find_uv_for_point(surface, point);
    surface.normal_at(u, v)
}

#[allow(dead_code)]
fn surface_normal_at_point(
    surface: &SurfaceGeometry,
    point: &popkern_core::geometry::Point3,
    brep: &BRep,
    face_handle: FaceHandle,
) -> Result<Vec3, String> {
    // Find UV coordinates by closest-point projection.
    let (u, v) = find_uv_for_point(surface, point);
    let mut normal = surface.normal_at(u, v);

    // Account for face orientation.
    let face = brep
        .faces
        .get(face_handle)
        .ok_or_else(|| format!("face {face_handle:?} not found"))?;
    if let Some(ref sof) = face.surface
        && !sof.orientation
    {
        normal = Vec3::new(-normal.x, -normal.y, -normal.z);
    }

    Ok(normal)
}

/// Find (u, v) parameters on a surface closest to a 3D point.
///
/// Uses analytic projection for planes and a coarse grid search for others.
fn find_uv_for_point(
    surface: &SurfaceGeometry,
    point: &popkern_core::geometry::Point3,
) -> (f64, f64) {
    // Try analytic projection first.
    if let Some(uv) = popkern_core::geometry::uv_projection::project_point_to_uv(surface, point) {
        return uv;
    }

    // Fallback: grid search over the parameter range.
    let ((u_min, u_max), (v_min, v_max)) = surface.parameter_range();
    let u_lo = if u_min.is_finite() { u_min } else { -10.0 };
    let u_hi = if u_max.is_finite() { u_max } else { 10.0 };
    let v_lo = if v_min.is_finite() { v_min } else { -10.0 };
    let v_hi = if v_max.is_finite() { v_max } else { 10.0 };

    let n = 20;
    let mut best_u = (u_lo + u_hi) * 0.5;
    let mut best_v = (v_lo + v_hi) * 0.5;
    let mut best_dist = f64::INFINITY;

    for iu in 0..n {
        let u = u_lo + (u_hi - u_lo) * (iu as f64) / ((n - 1).max(1) as f64);
        for iv in 0..n {
            let v = v_lo + (v_hi - v_lo) * (iv as f64) / ((n - 1).max(1) as f64);
            let sp = surface.point_at(u, v);
            let d = sp.distance_to(point);
            if d < best_dist {
                best_dist = d;
                best_u = u;
                best_v = v;
            }
        }
    }

    (best_u, best_v)
}

/// Whether `face`'s natural surface normal points out of the
/// material, from the outer loop's UV winding (shoelace over sampled
/// pcurves in loop order). Returns `None` when any loop half-edge
/// lacks a pcurve or the area is degenerate.
/// +1.0 when the shell's loop windings + orientation flags describe
/// an OUTWARD-oriented boundary, −1.0 for inward, via the divergence
/// theorem over outer-loop polygon fans (sign-exact for sane closed
/// shells). Primitive constructors currently build inward shells
/// (extrude keeps the base winding); boolean/sweep output builds
/// outward ones — consumers needing true outwardness must fold this
/// in until constructors are normalized (PK-100).
fn shell_outward_sign(brep: &BRep) -> f64 {
    use popkern_core::geometry::surface::SurfaceEval;
    // Divergence theorem per face using the FLAG-CORRECTED NATURAL
    // normal — deliberately NOT the loop winding: producers disagree
    // on winding-vs-flag consistency (the cube winds inward-wound
    // loops, the square pyramid outward ones, both with flag=true),
    // so the loop fan supplies only area and centroid, while the
    // direction comes from what the flags claim. V > 0 ⇔ the
    // flag-corrected normals describe an outward shell.
    let mut vol = 0.0;
    for (fh, face) in brep.faces.iter() {
        let Some(sof) = face.surface.as_ref() else {
            continue;
        };
        let flag_sign = if sof.orientation { 1.0 } else { -1.0 };
        let pts = brep.face_outer_loop_positions(fh);
        if pts.len() < 3 {
            continue;
        }
        let m = pts.len() as f64;
        let centroid = popkern_core::geometry::Point3::new(
            pts.iter().map(|p| p.x).sum::<f64>() / m,
            pts.iter().map(|p| p.y).sum::<f64>() / m,
            pts.iter().map(|p| p.z).sum::<f64>() / m,
        );
        let mut area_vec = Vec3::new(0.0, 0.0, 0.0);
        for i in 1..pts.len() - 1 {
            let v1 = pts[i] - pts[0];
            let v2 = pts[i + 1] - pts[0];
            area_vec = area_vec + v1.cross(&v2) * 0.5;
        }
        let area = area_vec.length();
        if area < 1e-16 {
            continue;
        }
        let (u, v) = find_uv_for_point(&sof.geometry, &centroid);
        let n = sof.geometry.normal_at(u, v) * flag_sign;
        let x = Vec3::new(centroid.x, centroid.y, centroid.z);
        vol += n.dot(&x) * area;
    }
    if std::env::var_os("PK_CONVEXITY_TRACE").is_some() {
        eprintln!("    [shell] signed 6V = {vol:.3e}");
    }
    if vol >= 0.0 { 1.0 } else { -1.0 }
}

fn effective_orientation_pcurves_missing(brep: &BRep, face: FaceHandle) -> bool {
    let Some(face_data) = brep.faces.get(face) else {
        return true;
    };
    let Some(l) = brep.loops.get(face_data.outer_loop) else {
        return true;
    };
    let start = l.half_edge;
    let mut he = start;
    for _ in 0..10_000 {
        let Some(data) = brep.half_edges.get(he) else {
            return true;
        };
        if data.pcurve.is_none() {
            return true;
        }
        he = data.next;
        if he == start {
            return false;
        }
    }
    true
}

fn effective_outward_orientation_centroid(brep: &BRep, face: FaceHandle) -> Option<bool> {
    use popkern_core::geometry::surface::SurfaceEval;
    let face_data = brep.faces.get(face)?;
    let surface = &face_data.surface.as_ref()?.geometry;
    let positions = brep.face_outer_loop_positions(face);
    if positions.len() < 3 {
        return None;
    }
    let m = positions.len() as f64;
    let fc = popkern_core::geometry::Point3::new(
        positions.iter().map(|p| p.x).sum::<f64>() / m,
        positions.iter().map(|p| p.y).sum::<f64>() / m,
        positions.iter().map(|p| p.z).sum::<f64>() / m,
    );
    let n = brep.vertices.len() as f64;
    if n < 4.0 {
        return None;
    }
    let mut bc = popkern_core::geometry::Point3::new(0.0, 0.0, 0.0);
    for (_, v) in brep.vertices.iter() {
        bc.x += v.position.x / n;
        bc.y += v.position.y / n;
        bc.z += v.position.z / n;
    }
    let (u, v) = find_uv_for_point(surface, &fc);
    let normal = surface.normal_at(u, v);
    let out = fc - bc;
    let d = normal.dot(&out);
    if d.abs() < 1e-12 {
        return None;
    }
    Some(d > 0.0)
}

fn effective_outward_orientation(brep: &BRep, face: FaceHandle) -> Option<bool> {
    let face_data = brep.faces.get(face)?;
    let dbg = std::env::var_os("PK_CONVEXITY_TRACE").is_some();
    if effective_orientation_pcurves_missing(brep, face) {
        if dbg {
            eprintln!(
                "    [eff] face {face:?}: pcurve-less -> centroid fallback = {:?}",
                effective_outward_orientation_centroid(brep, face)
            );
        }
        // Pcurve-less faces (primitive constructors): outwardness
        // from geometry — the natural normal at the face centroid
        // vs the body centroid. Exact for convex bodies, which is
        // the primitive-constructor population; boolean/sweep
        // output always carries pcurves and never reaches this.
        return effective_outward_orientation_centroid(brep, face);
    }
    let start = brep.loops.get(face_data.outer_loop)?.half_edge;
    let surface = &face_data.surface.as_ref()?.geometry;
    let mut pts: Vec<(f64, f64)> = Vec::new();
    let mut he = start;
    for _ in 0..10_000 {
        let data = brep.half_edges.get(he)?;
        let pcurve = data.pcurve.as_ref()?;
        let (t0, t1) = pcurve.parameter_range();
        // No repo-wide pcurve direction contract exists (primitive
        // constructors store edge-aligned pcurves; recomputed ones
        // are traversal-aligned). Orient each run by mapping the
        // pcurve start through the surface and comparing against
        // the half-edge's own endpoints; ties (closed edges) keep
        // the stored direction.
        let origin = brep.vertices.get(data.origin)?.position;
        let dest = brep
            .vertices
            .get(brep.half_edges.get(data.twin)?.origin)?
            .position;
        let uv0 = pcurve.point_at(t0);
        let q0 = surface.point_at(uv0.x, uv0.y);
        let flip = q0.distance_to(&dest) + 1.0e-12 < q0.distance_to(&origin);
        for k in 0..8 {
            let t = if flip {
                t1 + (t0 - t1) * k as f64 / 8.0
            } else {
                t0 + (t1 - t0) * k as f64 / 8.0
            };
            let p = pcurve.point_at(t);
            pts.push((p.x, p.y));
        }
        he = data.next;
        if he == start {
            break;
        }
    }
    if pts.len() < 3 {
        return None;
    }
    let mut area = 0.0;
    for i in 0..pts.len() {
        let j = (i + 1) % pts.len();
        area += pts[i].0 * pts[j].1 - pts[j].0 * pts[i].1;
    }
    if area.abs() < 1.0e-14 {
        return None;
    }
    if std::env::var_os("PK_CONVEXITY_TRACE").is_some() {
        eprintln!("    [eff] face {face:?}: shoelace area {area:.3e}");
    }
    Some(area > 0.0)
}

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
