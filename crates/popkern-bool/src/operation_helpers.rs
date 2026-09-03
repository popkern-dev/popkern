//! Helper functions for the boolean pipeline orchestrator: the BRep
//! copy/combine utilities, empty-selection handling, containment
//! verification and record tracing. The no-intersection dispatch
//! lives in `operation_no_intersection.rs` (500-line cap).

use std::collections::HashMap;

use popkern_core::topology::{BRep, FaceHandle};

use crate::error::{BoolError, BoolOperation, Classification};
use crate::split::SplitResult;

/// Handle the case where no intersection curves were found between the solids.
///
/// Determines the spatial relationship (disjoint, containment, or coincident)
/// by classifying a sample point from A against B, then returns the
/// appropriate result for the given operation.
/// Combine two disjoint BReps into a single BRep by copying all topology
/// from `b` into `a`.
///
/// The two inputs should not share geometry — this is a structural
/// concatenation that copies handles from `b` into fresh slotmap slots
/// in `a` without welding at intersection points. Used by the sheet
/// dispatcher to implement non-overlapping `SheetSheetUnion`.
pub(crate) fn combine_breps(a: BRep, b: BRep) -> Result<BRep, BoolError> {
    copy_brep_into(a, b, false)
}

/// Handle the case where the select stage produces 0 kept faces.
///
/// This happens when all fragments from both operands classify the same way
/// (e.g. all Inside). Infers containment from the classification maps and
/// returns the appropriate operand directly. When containment cannot be
/// determined, returns an empty BRep.
pub(super) fn handle_empty_selection(
    op: BoolOperation,
    split_result: SplitResult,
    classifications_a: &HashMap<FaceHandle, Classification>,
    classifications_b: &HashMap<FaceHandle, Classification>,
) -> BRep {
    let all_a_inside = !classifications_a.is_empty()
        && classifications_a
            .values()
            .all(|c| *c == Classification::Inside);
    let all_b_inside = !classifications_b.is_empty()
        && classifications_b
            .values()
            .all(|c| *c == Classification::Inside);

    if all_a_inside {
        // A is entirely inside B.
        return match op {
            BoolOperation::Union => split_result.b,
            BoolOperation::Intersection => split_result.a,
            BoolOperation::Subtraction => BRep::new(split_result.a.tolerance),
        };
    }

    if all_b_inside {
        // B is entirely inside A.
        return match op {
            BoolOperation::Union => split_result.a,
            BoolOperation::Intersection => split_result.b,
            // A minus B when B is inside A — full hollow solid requires
            // combining A with flipped B. For now, return A (the outer shell).
            BoolOperation::Subtraction => split_result.a,
        };
    }

    // Cannot determine containment — return empty result.
    BRep::new(split_result.a.tolerance)
}

/// Copy all topology from `b` into `a`, optionally flipping B's face normals.
///
/// When `flip_b_normals` is true, the surface orientation on every copied face
/// is inverted, producing an inner void shell for hollow-solid subtraction.
pub(crate) fn copy_brep_into(
    mut a: BRep,
    b: BRep,
    flip_b_normals: bool,
) -> Result<BRep, BoolError> {
    use popkern_core::topology::{Edge, Face, HalfEdge, Loop, Shell, Solid, Vertex};

    // Build handle maps for B → new handles in A.
    let mut v_map: HashMap<
        popkern_core::topology::VertexHandle,
        popkern_core::topology::VertexHandle,
    > = HashMap::new();
    let mut he_map: HashMap<
        popkern_core::topology::HalfEdgeHandle,
        popkern_core::topology::HalfEdgeHandle,
    > = HashMap::new();
    let mut e_map: HashMap<popkern_core::topology::EdgeHandle, popkern_core::topology::EdgeHandle> =
        HashMap::new();
    let mut l_map: HashMap<popkern_core::topology::LoopHandle, popkern_core::topology::LoopHandle> =
        HashMap::new();
    let mut f_map: HashMap<popkern_core::topology::FaceHandle, popkern_core::topology::FaceHandle> =
        HashMap::new();
    let mut sh_map: HashMap<
        popkern_core::topology::ShellHandle,
        popkern_core::topology::ShellHandle,
    > = HashMap::new();
    let mut so_map: HashMap<
        popkern_core::topology::SolidHandle,
        popkern_core::topology::SolidHandle,
    > = HashMap::new();

    // Copy vertices.
    for (old_h, v) in &b.vertices {
        let new_h = a.vertices.insert(Vertex {
            position: v.position,
            outgoing: popkern_core::topology::HalfEdgeHandle::default(),
            name: v.name,
            tolerance: None,
        });
        v_map.insert(old_h, new_h);
    }

    // Copy edges (with placeholder half_edge).
    for (old_h, e) in &b.edges {
        let new_h = a.insert_edge_tracked(Edge {
            half_edge: popkern_core::topology::HalfEdgeHandle::default(),
            curve: e.curve.clone(),
            name: e.name,
            tolerance: None,
            canonical: e.canonical,
            role: e.role,
        });
        e_map.insert(old_h, new_h);
    }

    // Copy half-edges (with placeholders).
    for (old_h, _he) in &b.half_edges {
        let new_h = a.half_edges.insert(HalfEdge {
            origin: popkern_core::topology::VertexHandle::default(),
            twin: popkern_core::topology::HalfEdgeHandle::default(),
            next: popkern_core::topology::HalfEdgeHandle::default(),
            prev: popkern_core::topology::HalfEdgeHandle::default(),
            edge: popkern_core::topology::EdgeHandle::default(),
            loop_: popkern_core::topology::LoopHandle::default(),
            pcurve: None,
        });
        he_map.insert(old_h, new_h);
    }

    // Copy loops (with placeholders).
    for (old_h, _l) in &b.loops {
        let new_h = a.loops.insert(Loop {
            half_edge: popkern_core::topology::HalfEdgeHandle::default(),
            face: popkern_core::topology::FaceHandle::default(),
        });
        l_map.insert(old_h, new_h);
    }

    // Copy solids (with placeholders).
    for (old_h, _s) in &b.solids {
        let new_h = a.solids.insert(Solid {
            shell: popkern_core::topology::ShellHandle::default(),
        });
        so_map.insert(old_h, new_h);
    }

    // Copy shells (with placeholders).
    for (old_h, _s) in &b.shells {
        let new_h = a.shells.insert(Shell {
            face: popkern_core::topology::FaceHandle::default(),
            solid: popkern_core::topology::SolidHandle::default(),
        });
        sh_map.insert(old_h, new_h);
    }

    // Copy faces, optionally flipping surface orientation.
    for (old_h, f) in &b.faces {
        let surface = if flip_b_normals {
            f.surface
                .as_ref()
                .map(|s| popkern_core::geometry::SurfaceOnFace {
                    geometry: s.geometry.clone(),
                    orientation: !s.orientation,
                })
        } else {
            f.surface.clone()
        };
        let new_h = a.faces.insert(Face {
            outer_loop: popkern_core::topology::LoopHandle::default(),
            inner_loops: Vec::new(),
            shell: popkern_core::topology::ShellHandle::default(),
            surface,
            name: f.name,
            tolerance: None,
            role: Default::default(),
        });
        f_map.insert(old_h, new_h);
    }

    // Now fix up all the cross-references.

    // Fix half-edges. Source BReps that came out of `BRepBuilder`
    // without a full loop wiring can have default (unset) handles
    // in next/prev/twin — look those up via `.get()` so we pass
    // them through as default in the target instead of panicking
    // on an unknown key. The source's brokenness is preserved, but
    // the copy itself doesn't panic.
    for (old_h, he) in &b.half_edges {
        let new_h = he_map[&old_h];
        a.half_edges[new_h].origin = v_map.get(&he.origin).copied().unwrap_or_default();
        a.half_edges[new_h].twin = he_map.get(&he.twin).copied().unwrap_or_default();
        a.half_edges[new_h].next = he_map.get(&he.next).copied().unwrap_or_default();
        a.half_edges[new_h].prev = he_map.get(&he.prev).copied().unwrap_or_default();
        a.half_edges[new_h].edge = e_map.get(&he.edge).copied().unwrap_or_default();
        a.half_edges[new_h].loop_ = l_map.get(&he.loop_).copied().unwrap_or_default();
        a.half_edges[new_h].pcurve = he.pcurve.clone();
    }

    // Fix vertices (outgoing). Same tolerance for broken inputs as
    // the half-edge fixup above.
    for (old_h, v) in &b.vertices {
        let new_h = v_map[&old_h];
        if b.half_edges.contains_key(v.outgoing) {
            a.vertices[new_h].outgoing = he_map.get(&v.outgoing).copied().unwrap_or_default();
        }
    }

    // Fix edges (half_edge).
    for (old_h, e) in &b.edges {
        let new_h = e_map[&old_h];
        a.edges[new_h].half_edge = he_map.get(&e.half_edge).copied().unwrap_or_default();
    }

    // Fix loops.
    for (old_h, l) in &b.loops {
        let new_h = l_map[&old_h];
        a.loops[new_h].half_edge = he_map.get(&l.half_edge).copied().unwrap_or_default();
        a.loops[new_h].face = f_map.get(&l.face).copied().unwrap_or_default();
    }

    // Fix faces.
    for (old_h, f) in &b.faces {
        let new_h = f_map[&old_h];
        a.faces[new_h].outer_loop = l_map.get(&f.outer_loop).copied().unwrap_or_default();
        a.faces[new_h].shell = sh_map.get(&f.shell).copied().unwrap_or_default();
        a.faces[new_h].inner_loops = f
            .inner_loops
            .iter()
            .map(|il| l_map.get(il).copied().unwrap_or_default())
            .collect();
    }

    // Fix shells.
    for (old_h, s) in &b.shells {
        let new_h = sh_map[&old_h];
        a.shells[new_h].face = f_map.get(&s.face).copied().unwrap_or_default();
        a.shells[new_h].solid = so_map.get(&s.solid).copied().unwrap_or_default();
    }

    // Fix solids.
    for (old_h, s) in &b.solids {
        let new_h = so_map[&old_h];
        a.solids[new_h].shell = sh_map.get(&s.shell).copied().unwrap_or_default();
    }

    Ok(a)
}

/// the boolean trace topic helper: endpoint/chord dump of every intersection
/// record (a fragmented closed curve shows up as short-chord copies).
pub(crate) fn trace_records(records: &[crate::intersect::IntersectionRecord], n_overlaps: usize) {
    eprintln!("  [exec] {} records, {n_overlaps} overlaps", records.len());
    for (i, r) in records.iter().enumerate() {
        use popkern_core::geometry::curve::CurveEval as _;
        let c = &r.curve.curve_3d;
        let (t0, t1) = c.parameter_range();
        let s = c.point_at(t0);
        let e = c.point_at(t1);
        let uv_span = |c: &popkern_core::geometry::NurbsCurve| -> (f64, f64) {
            let (t0, t1) = c.parameter_range();
            let (mut lo, mut hi) = (f64::MAX, f64::MIN);
            for k in 0..=16 {
                let u = c.point_at(t0 + (t1 - t0) * k as f64 / 16.0).x;
                lo = lo.min(u);
                hi = hi.max(u);
            }
            (lo, hi)
        };
        let (u1lo, u1hi) = uv_span(&r.curve.curve_uv1);
        let (u2lo, u2hi) = uv_span(&r.curve.curve_uv2);
        eprintln!(
            "    [rec {i}] ({:?},{:?}) ({:.6},{:.6},{:.6})->({:.6},{:.6},{:.6}) chord {:.4} u1[{u1lo:.4},{u1hi:.4}] u2[{u2lo:.4},{u2hi:.4}]",
            r.face_a,
            r.face_b,
            s.x,
            s.y,
            s.z,
            e.x,
            e.y,
            e.z,
            s.distance_to(&e)
        );
    }
}

/// Verify a claimed whole-body containment directly against the
/// ORIGINAL operands: every vertex of `inner` must classify Inside
/// or OnBoundary of `outer`. The empty-selection shortcut's claim
/// comes from seed-propagated face classifications, and a silently
/// under-split operand poisons every seed (measured: a rim-straddling
/// key subtract left the disc's top/bottom faces unsplit, one seed
/// wiped each side, and the "contained" flange came back as an EMPTY
/// body that passed every downstream gate vacuously — PK-173).
pub(super) fn verify_whole_body_containment(inner: &BRep, outer: &BRep) -> bool {
    use crate::classify::point_in_solid_final;
    use crate::error::PointClassification;
    use crate::intersect::compute_face_aabbs;
    let aabbs = compute_face_aabbs(outer);
    let tolerance = outer.tolerance;
    // Vertices alone under-cover analytic bodies (a plain cylinder
    // has only its two seam vertices — both can sit inside a partial
    // tool while most of the body is outside), so every face's
    // boundary is curve-sampled as well.
    let mut probes: Vec<popkern_core::geometry::Point3> =
        inner.vertices.values().map(|v| v.position).collect();
    for face in inner.faces.keys() {
        probes.extend(inner.face_outer_loop_sample_positions(face, 8));
    }
    if probes.is_empty() {
        return false;
    }
    for probe in &probes {
        match point_in_solid_final(probe, outer, &aabbs, &tolerance) {
            Ok(PointClassification::Inside) | Ok(PointClassification::OnBoundary) => {}
            _ => return false,
        }
    }
    true
}
