//! A sweep along a path of several edges (brep-bench §sweep): the
//! edges chained end to end from the path edge, each run swept with
//! one frame carried across the joints (the rotation-minimising frame
//! of each run, its start the end of the last turned by the least
//! rotation from one tangent to the next), so the section keeps its
//! attitude round every corner.
//!
//! Where runs meet smoothly their end sections coincide and the runs
//! join. At a sharp corner (the tangent turning by more than 1e-9
//! rad): with `Mitre` each run is carried straight on past the corner
//! and cut by the plane bisecting it, `(p - c) · (a + b) = 0` for the
//! incoming tangent `a` and outgoing `b`, so the runs meet in it; with
//! `Round` the section turns about the corner (the axis through it
//! along `a × b`) from one run to the next. The pieces are united by
//! the boolean (they share faces flush). A U-turn has no bisecting
//! plane and fails; a twist is not distributed over several runs and
//! fails. A closed chain turns at its start too. Faces are named `Generated i` in face order (as a sweep's).

use popkern_core::geometry::{NurbsCurve, Point3, Vec3};
use popkern_core::history::TrackingContext;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;
use popkern_sweep::{FrameMethod, SweepOptions, SweepProfile};

use crate::curve_chain::ChainPiece;
use crate::error::ReplayError;
use crate::feature_sweep_path::SweepCorner;

fn fail(detail: impl Into<String>) -> ReplayError {
    ReplayError::FeatureFailed {
        label: "sweep".into(),
        error: crate::failure::FeatureFailure::sweep("sweep", &detail.into()).into(),
    }
}

/// `v` turned by the least rotation taking unit `from` to unit `to`.
fn carry(v: Vec3, from: Vec3, to: Vec3) -> Vec3 {
    let axis = from.cross(&to);
    let s = axis.length();
    let c = from.dot(&to);
    if s < 1e-15 {
        return v;
    }
    let k = axis * (1.0 / s);
    let angle = s.patan2(c);
    v * angle.pcos() + k.cross(&v) * angle.psin() + k * (k.dot(&v) * (1.0 - angle.pcos()))
}

/// The straight path from `a` to `b`.
fn segment(a: Point3, b: Point3) -> NurbsCurve {
    NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![a, b],
        weights: vec![1.0, 1.0],
    }
}

fn sweep_one(
    profile: &SweepProfile,
    path: &NurbsCurve,
    frame: FrameMethod,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let options = SweepOptions {
        frame,
        ..SweepOptions::default()
    };
    popkern_sweep::sweep(profile, path, &options, tolerance, None)
        .map(|r| r.brep)
        .map_err(|e| fail(format!("sweeping a run of the path: {e}")))
}

/// The run's start normal: Bishop's own seed for the first run, the
/// carried normal after.
fn start_frame(normal: Option<Vec3>) -> FrameMethod {
    match normal {
        None => FrameMethod::Bishop,
        Some(n) => FrameMethod::BishopFrom(n),
    }
}

/// The normal at the end of a run swept from `frame`.
fn end_normal(path: &NurbsCurve, frame: &FrameMethod) -> Result<Vec3, ReplayError> {
    let table = popkern_sweep::arc_length::ArcLengthTable::build(path, 200);
    let (seq, _) = popkern_sweep::frame::compute_frames(path, &table, frame, 200, 0.0, false)
        .map_err(|e| fail(format!("carrying the frame: {e}")))?;
    seq.frames
        .last()
        .map(|f| f.normal)
        .ok_or_else(|| fail("a run of the path has no frame"))
}

/// The half-space `(p - c) · n >= 0` as a box reaching `reach` from `c`.
fn half_space(c: Point3, n: Vec3, reach: f64, tolerance: Tolerance) -> Result<BRep, ReplayError> {
    let n = n.normalized().ok_or_else(|| fail("a degenerate mitre"))?;
    let seed = if n.x.abs() < 0.9 {
        Vec3::new(1.0, 0.0, 0.0)
    } else {
        Vec3::new(0.0, 1.0, 0.0)
    };
    let u = (seed - n * seed.dot(&n))
        .normalized()
        .ok_or_else(|| fail("a degenerate mitre"))?;
    let v = n.cross(&u);
    let square = |base: Point3| {
        [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)]
            .iter()
            .map(|(a, b)| base + u * (a * reach) + v * (b * reach))
            .collect::<Vec<_>>()
    };
    BRep::polygon_frustum(&square(c), &square(c + n * reach), tolerance)
        .map_err(|e| fail(format!("the mitre's cutter: {e}")))
}

fn boolean(
    op: popkern_bool::BoolOperation,
    a: BRep,
    b: BRep,
    tolerance: Tolerance,
    what: &str,
) -> Result<BRep, ReplayError> {
    popkern_bool::boolean(op, a, b, tolerance, None).map_err(|e| fail(format!("{what}: {e}")))
}

/// The profile's farthest reach from its centroid (bounds the mitre's
/// run-on past a corner).
fn profile_reach(profile: &SweepProfile, tolerance: &Tolerance) -> Result<f64, ReplayError> {
    let analyzed = popkern_sweep::profile::analyze_profile(profile, tolerance)
        .map_err(|e| fail(format!("the profile: {e}")))?;
    let mut reach: f64 = 0.0;
    for span in &analyzed.spans {
        for p in &span.control_points {
            reach = reach.max((p.x * p.x + p.y * p.y).sqrt());
        }
    }
    Ok(reach)
}

/// Sweep `profile` along the chained `pieces`.
pub(super) fn sweep_chain(
    profile: &SweepProfile,
    pieces: &[ChainPiece],
    corner: SweepCorner,
    tracking: &mut TrackingContext,
    tolerance: Tolerance,
) -> Result<BRep, ReplayError> {
    let paths: Vec<NurbsCurve> = pieces.iter().map(ChainPiece::nurbs).collect();
    let n = pieces.len();
    // The tangents in and out of each joint, and whether it is sharp.
    let joints: Vec<(Point3, Vec3, Vec3)> = (0..n - 1)
        .map(|i| {
            (
                pieces[i].ends().1,
                pieces[i].tangent(1.0),
                pieces[i + 1].tangent(0.0),
            )
        })
        .collect();
    let sharp = |(_, a, b): &(Point3, Vec3, Vec3)| a.cross(b).length() > 1e-9 || a.dot(b) < 0.0;
    // A closed chain also turns where it ends: the joint from the last
    // run back into the first.
    let closed = n > 1
        && pieces[0].ends().0.distance_to(&pieces[n - 1].ends().1) <= tolerance.linear() * 10.0;
    let wrap = closed.then(|| {
        (
            pieces[0].ends().0,
            pieces[n - 1].tangent(1.0),
            pieces[0].tangent(0.0),
        )
    });
    let before = |i: usize| if i > 0 { Some(joints[i - 1]) } else { wrap };
    let after = |i: usize| if i + 1 < n { Some(joints[i]) } else { wrap };
    for (c, a, b) in joints.iter().chain(wrap.iter()) {
        if (*a + *b).length() < 1e-9 {
            return Err(fail(format!(
                "the path turns straight back at ({:.6}, {:.6}, {:.6}): no plane bisects a U-turn",
                c.x, c.y, c.z
            )));
        }
    }
    let reach = profile_reach(profile, &tolerance)?;
    let size = paths
        .iter()
        .flat_map(|p| p.control_points.iter())
        .fold(0.0f64, |m, p| {
            m.max((*p - Point3::new(0.0, 0.0, 0.0)).length())
        })
        + reach;
    let cutter_reach = 4.0 * size + 1.0;

    // Frames: each run starts with the last run's end normal carried
    // across the joint.
    let mut starts: Vec<Option<Vec3>> = vec![None];
    let mut ends: Vec<Vec3> = Vec::with_capacity(n);
    for i in 0..n {
        let end = end_normal(&paths[i], &start_frame(starts[i]))?;
        ends.push(end);
        if i + 1 < n {
            starts.push(Some(carry(end, joints[i].1, joints[i].2)));
        }
    }

    let mut result: Option<BRep> = None;
    for i in 0..n {
        // Carry the run straight on past each sharp corner, then cut it
        // back to its side of the bisecting plane. A straight run is
        // swept along its extended line; a curved one is joined to a
        // straight run-on (two coaxial runs joined end to end made a
        // tangent seam the union declined).
        let run_on = 2.0 * reach + 10.0 * tolerance.linear();
        let mitre = corner == SweepCorner::Mitre;
        let lead_joint = before(i).filter(|j| mitre && sharp(j));
        let tail_joint = after(i).filter(|j| mitre && sharp(j));
        let straight = matches!(
            pieces[i].curve.geometry,
            popkern_core::geometry::curve::CurveGeometry::Line { .. }
        );
        let mut run = if straight && (lead_joint.is_some() || tail_joint.is_some()) {
            let (p0, p1) = pieces[i].ends();
            let d = pieces[i].tangent(0.0);
            let from = if lead_joint.is_some() {
                p0 - d * run_on
            } else {
                p0
            };
            let to = if tail_joint.is_some() {
                p1 + d * run_on
            } else {
                p1
            };
            let normal_here = start_normal(&paths[i], starts[i])?;
            sweep_one(
                profile,
                &segment(from, to),
                FrameMethod::Fixed(normal_here),
                tolerance,
            )?
        } else {
            sweep_one(profile, &paths[i], start_frame(starts[i]), tolerance)?
        };
        if mitre {
            if let Some(joint) = lead_joint {
                let (c, a, b) = joint;
                if !straight {
                    let normal_here = start_normal(&paths[i], starts[i])?;
                    let lead = sweep_one(
                        profile,
                        &segment(c - b * run_on, c),
                        FrameMethod::Fixed(normal_here),
                        tolerance,
                    )?;
                    run = boolean(
                        popkern_bool::BoolOperation::Union,
                        run,
                        lead,
                        tolerance,
                        &format!("carrying run {i} back past its corner"),
                    )?;
                }
                let keep_out = half_space(c, -(a + b), cutter_reach, tolerance)?;
                run = boolean(
                    popkern_bool::BoolOperation::Subtraction,
                    run,
                    keep_out,
                    tolerance,
                    &format!("mitring run {i}"),
                )?;
            }
            if let Some(joint) = tail_joint {
                let (c, a, b) = joint;
                if !straight {
                    let tail = sweep_one(
                        profile,
                        &segment(c, c + a * run_on),
                        FrameMethod::Fixed(ends[i]),
                        tolerance,
                    )?;
                    run = boolean(
                        popkern_bool::BoolOperation::Union,
                        run,
                        tail,
                        tolerance,
                        &format!("carrying run {i} on past its corner"),
                    )?;
                }
                let keep_out = half_space(c, a + b, cutter_reach, tolerance)?;
                run = boolean(
                    popkern_bool::BoolOperation::Subtraction,
                    run,
                    keep_out,
                    tolerance,
                    &format!("mitring run {i}"),
                )?;
            }
        }
        result = Some(match result {
            None => run,
            Some(body) => boolean(
                popkern_bool::BoolOperation::Union,
                body,
                run,
                tolerance,
                &format!("joining run {i} to the runs before it"),
            )?,
        });
        if let Some(joint) = after(i).filter(|j| corner == SweepCorner::Round && sharp(j)) {
            let (c, a, b) = joint;
            let knuckle =
                super::replay_sweep_corner::knuckle(profile, c, (a, b), ends[i], tolerance)?;
            let body = result.take().ok_or_else(|| fail("no run"))?;
            result = Some(boolean(
                popkern_bool::BoolOperation::Union,
                body,
                knuckle,
                tolerance,
                &format!("rounding the corner after run {i}"),
            )?);
        }
    }
    let mut body = result.ok_or_else(|| fail("the path has no runs"))?;
    name_faces(&mut body, tracking);
    Ok(body)
}

/// The start normal of a run swept from `seed`.
fn start_normal(path: &NurbsCurve, seed: Option<Vec3>) -> Result<Vec3, ReplayError> {
    let table = popkern_sweep::arc_length::ArcLengthTable::build(path, 8);
    let (seq, _) =
        popkern_sweep::frame::compute_frames(path, &table, &start_frame(seed), 8, 0.0, false)
            .map_err(|e| fail(format!("carrying the frame: {e}")))?;
    seq.frames
        .first()
        .map(|f| f.normal)
        .ok_or_else(|| fail("a run of the path has no frame"))
}

/// `Generated i` names over the result's faces, edges and vertices.
fn name_faces(brep: &mut BRep, ctx: &mut TrackingContext) {
    use popkern_core::naming::{EntityName, EntityRole, RoleCategory};
    use popkern_core::naming_criteria::{
        GeometricCriteria, capture_edge, capture_face, capture_vertex,
    };
    let mut mint = |index: usize| {
        let name = EntityName {
            operation: ctx.operation,
            role: EntityRole {
                category: RoleCategory::Generated,
                index: index as u32,
            },
            ancestors: Vec::new(),
        };
        ctx.record_generated(name.clone());
        name
    };
    let faces: Vec<_> = brep.faces.keys().collect();
    for (i, fh) in faces.into_iter().enumerate() {
        let criteria = capture_face(brep, fh).map(GeometricCriteria::Face);
        let id = brep.name_pool.intern(mint(i), criteria);
        brep.faces[fh].name = Some(id);
    }
    let edges: Vec<_> = brep.edges.keys().collect();
    for (i, eh) in edges.into_iter().enumerate() {
        let criteria = capture_edge(brep, eh).map(GeometricCriteria::Edge);
        let id = brep.name_pool.intern(mint(i), criteria);
        brep.edges[eh].name = Some(id);
    }
    let vertices: Vec<_> = brep.vertices.keys().collect();
    for (i, vh) in vertices.into_iter().enumerate() {
        let criteria = capture_vertex(brep, vh).map(GeometricCriteria::Vertex);
        let id = brep.name_pool.intern(mint(i), criteria);
        brep.vertices[vh].name = Some(id);
    }
}
