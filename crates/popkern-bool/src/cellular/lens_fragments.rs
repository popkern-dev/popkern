//! A lens fragment takes its side from its neighbours (W3 V cells,
//! tangency design §4.8).
//!
//! A fragment no point of which clears the other operand past the
//! linear reach (a lens between two records, every point of it within
//! the tolerance of the other face) reads its side by noise: its sample
//! is as near the other operand as the tolerance. Its records cut the
//! face, and crossing a record flips the side of the other operand, so
//! a decided neighbour across any of them says its side; where every
//! decided neighbour agrees, that side is the fragment's (measured,
//! CC-V 64: two arms of a tangent crossing running 3e-4 to B's rim
//! bounded a lens 1e-7 deep on both cones; A's read Outside beside an
//! Outside neighbour across an arm, B's read Outside across from Inside,
//! and the union kept both sides of the lens and broke its loops). The
//! propagation planner leaves these out: a crossing's arms are too
//! near tangent for its transversality gate.

use popkern_math::hash::FastHashMap;

use crate::error::PointClassification;
use popkern_core::geometry::Point3;

use super::fragment_classify_one::ClassifyCtx;
use crate::merge::graph::{EdgeSource, FaceRef};
use crate::merge::walk::Fragment;

/// One fragment's verdict: sample, class, co-orientation.
type Verdict = (Point3, PointClassification, Option<bool>);

fn flipped(class: PointClassification) -> Option<PointClassification> {
    match class {
        PointClassification::Inside => Some(PointClassification::Outside),
        PointClassification::Outside => Some(PointClassification::Inside),
        _ => None,
    }
}

/// Whether record half-edge `h` crosses: its two faces' normals at its
/// middle part by more than the angular tolerance. A tangent record (a
/// curve the faces touch along) bounds the same side on both of its
/// sides, and nothing flips across it (a dimple's rim touching a face
/// at milli scale, `tangency_corpus` cone-dimple-touching-edge).
fn transversal(ctx: &ClassifyCtx<'_>, h: usize) -> bool {
    use popkern_core::geometry::curve::CurveEval as _;
    use popkern_core::geometry::surface::SurfaceEval as _;
    let graph = ctx.graph;
    let he = &graph.half_edges[h];
    let EdgeSource::Shared {
        face_a,
        face_b,
        operands,
    } = he.source
    else {
        return false;
    };
    let Some(curve) = he
        .curve
        .as_ref()
        .or_else(|| graph.half_edges[he.twin].curve.as_ref())
    else {
        return false;
    };
    let mid = curve.geometry.point_at(0.5 * (curve.t_start + curve.t_end));
    let normal = |operand: crate::merge::graph::Operand, face| {
        let brep = if operand == crate::merge::graph::Operand::A {
            ctx.brep_a
        } else {
            ctx.brep_b
        };
        let s = &brep.faces.get(face)?.surface.as_ref()?.geometry;
        let (u, v) = s.closest_point(&mid)?;
        s.normal_at(u, v).normalized()
    };
    normal(operands[0], face_a)
        .zip(normal(operands[1], face_b))
        .is_some_and(|(na, nb)| na.cross(&nb).length() > ctx.tolerance.angular())
}

/// Give each lens fragment the side its decided neighbours across its
/// records agree on (see the module doc).
///
/// Only fragments sampled on their own (not deferred to a propagation
/// component), without holes and bounded by records alone are tested:
/// the test projects onto the other operand's faces, which costs the
/// most on the multi-wrap walls the planner exists to spare.
pub(super) fn reclassify_lenses(
    ctx: &ClassifyCtx<'_>,
    fragments: &[Fragment],
    (degenerate, deferred): (&[bool], &[Option<usize>]),
    results: &mut [Option<Verdict>],
) {
    let graph = ctx.graph;
    let bounded_by_records = |f: &Fragment| {
        f.inner_loops.is_empty()
            && f.outer_loop.len() >= 2
            && f.outer_loop
                .iter()
                .all(|&h| matches!(graph.half_edges[h].source, EdgeSource::Shared { .. }))
    };
    let lens: Vec<bool> = fragments
        .iter()
        .enumerate()
        .map(|(i, f)| {
            !degenerate[i]
                && deferred[i].is_none()
                && bounded_by_records(f)
                && results[i].is_some_and(|r| {
                    flipped(r.1).is_some() && super::fragment_sample_deep::is_lens(ctx, f, r.0)
                })
        })
        .collect();
    if !lens.iter().any(|l| *l) {
        return;
    }
    let mut owner: FastHashMap<(FaceRef, usize), usize> = FastHashMap::default();
    for (i, f) in fragments.iter().enumerate() {
        for &h in f.outer_loop.iter().chain(f.inner_loops.iter().flatten()) {
            owner.insert((f.face, h), i);
        }
    }
    for (i, f) in fragments.iter().enumerate() {
        if !lens[i] {
            continue;
        }
        let mut implied: Option<PointClassification> = None;
        let mut agree = true;
        for &h in &f.outer_loop {
            let he = &graph.half_edges[h];
            if !transversal(ctx, h) {
                continue;
            }
            let Some(&j) = owner.get(&(f.face, he.twin)) else {
                continue;
            };
            if j == i || lens[j] {
                continue;
            }
            let Some(side) = results[j].and_then(|r| flipped(r.1)) else {
                continue;
            };
            match implied {
                None => implied = Some(side),
                Some(s) if s != side => agree = false,
                _ => {}
            }
        }
        let (Some(side), true) = (implied, agree) else {
            continue;
        };
        if let Some(r) = results[i].as_mut()
            && r.1 != side
        {
            if popkern_core::trace::trace_enabled("cellular") {
                eprintln!(
                    "  [lens] {:?} outer={} read {:?}, its neighbours across its records say {side:?}",
                    f.face,
                    f.outer_loop.len(),
                    r.1
                );
            }
            r.1 = side;
        }
    }
}
