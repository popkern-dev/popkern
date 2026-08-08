//! PK-142 stage 2 — TrimBoth attach for face–face blends: trim both
//! support walls back to the marched contact curves, remove the
//! enclosed gap region, and seam the strip in as a single face.
//!
//! The surgery is transactional: it operates on a private copy of the
//! input body and either returns a validated result or declines with
//! a structured [`BlendError::AttachDeclined`] naming the stage and
//! the measured evidence. Stages:
//!
//! 1. Clip the marched chain to both walls' real trim polygons and
//!    resample it uniformly along the centre polyline — uniform
//!    sections keep the degree-1 skin parameterization consistent
//!    with the trim-curve parameterization authored on the walls.
//! 2. Split each wall face along its contact polyline, and each strip
//!    end's supporting face along the section end arc.
//! 3. Flood-collect the doomed gap region (bounded by exactly the
//!    four split edges), merge it to a single face with
//!    `merge_faces_across_edge`, then dissolve leftover interior spur
//!    edges with `dissolve_same_loop_bridge`. Maximal merging leaves
//!    an acyclic spur forest whose trees touch the region boundary at
//!    most once (a leftover boundary-to-boundary path would separate
//!    the region, contradicting the merge spanning tree), so tip
//!    dissolution always terminates.
//! 4. Re-surface the merged face with the strip skin and author its
//!    chart-border pcurves — one geometric authority from
//!    construction (the PK-148 lesson).

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::contact::ContactSolution;
use crate::error::BlendError;
use crate::face_blend::{FaceBlendAttach, FaceBlendResult, FaceBlendSpec};
use crate::face_blend_attach_geometry::{decline, face_uv_polygon, point_in_polygon};
use crate::face_blend_attach_surgery::surgery;
use crate::face_blend_sheet::strip_sheet_body;

/// Re-solves one blend section with the ball centre pinned near a
/// target point, warm-started from a neighbouring solution.
pub(crate) type SectionResolver<'a> =
    dyn Fn(Point3, &ContactSolution) -> Option<ContactSolution> + 'a;

/// Result of a completed TrimBoth attach.
pub(crate) struct AttachOutcome {
    /// The modified body (walls trimmed, strip seamed in).
    pub body: BRep,
    /// The strip face inside `body`.
    pub strip_face: FaceHandle,
    /// The clipped, uniformly resampled section chain.
    pub chain: Vec<ContactSolution>,
}

/// Shared tail for both march paths: emit the sheet and, for
/// TrimBoth, run the attach surgery.
pub(crate) fn finish_face_blend(
    brep: &BRep,
    spec: &FaceBlendSpec,
    chain: Vec<ContactSolution>,
    resolve: &SectionResolver,
    bulge: Option<Vec3>,
    tolerance: Tolerance,
) -> Result<FaceBlendResult, BlendError> {
    match spec.attach {
        FaceBlendAttach::SheetOnly => {
            let contacts_a = chain.iter().map(|s| s.contact0).collect();
            let contacts_b = chain.iter().map(|s| s.contact1).collect();
            let sheet = strip_sheet_body(&chain, spec.radius, bulge, tolerance)?;
            Ok(FaceBlendResult {
                sheet,
                attached: None,
                strip_face: None,
                contacts_a,
                contacts_b,
            })
        }
        FaceBlendAttach::TrimBoth => {
            let mut outcome = attach_trim_both(brep, spec, &chain, resolve, bulge, tolerance)?;
            // Serving-boundary flag reconcile (PK-244/PK-245): the
            // orientation flag is derived state, and the cap-extension
            // surgery can author a face whose flag contradicts its
            // committed loop winding (measured: through-notch and
            // deep-channel cap extensions each shipped one CW loop
            // with flag true). Derive the honest value from the final
            // winding, as the fillet/chamfer/boolean tails do.
            outcome.body.reconcile_orientation_flags();
            let contacts_a = outcome.chain.iter().map(|s| s.contact0).collect();
            let contacts_b = outcome.chain.iter().map(|s| s.contact1).collect();
            let sheet = strip_sheet_body(&outcome.chain, spec.radius, bulge, tolerance)?;
            Ok(FaceBlendResult {
                sheet,
                attached: Some(outcome.body),
                strip_face: Some(outcome.strip_face),
                contacts_a,
                contacts_b,
            })
        }
    }
}

fn attach_trim_both(
    brep: &BRep,
    spec: &FaceBlendSpec,
    chain: &[ContactSolution],
    resolve: &SectionResolver,
    bulge: Option<Vec3>,
    tolerance: Tolerance,
) -> Result<AttachOutcome, BlendError> {
    let poly_a = face_uv_polygon(brep, spec.face_a)?;
    let poly_b = face_uv_polygon(brep, spec.face_b)?;
    let clipped = clip_and_resample(chain, &poly_a, &poly_b, spec.radius, resolve, tolerance)?;
    surgery(brep, spec, &clipped, bulge, tolerance)
}

fn inside_both(sol: &ContactSolution, poly_a: &[(f64, f64)], poly_b: &[(f64, f64)]) -> bool {
    point_in_polygon(sol.uv0, poly_a) && point_in_polygon(sol.uv1, poly_b)
}

pub(crate) fn lerp(a: Point3, b: Point3, s: f64) -> Point3 {
    Point3::new(
        a.x + (b.x - a.x) * s,
        a.y + (b.y - a.y) * s,
        a.z + (b.z - a.z) * s,
    )
}

/// Bisect the section chain between an inside and an outside sample
/// until the section sits on the trim boundary.
fn bisect_to_boundary(
    inside: &ContactSolution,
    outside: &ContactSolution,
    poly_a: &[(f64, f64)],
    poly_b: &[(f64, f64)],
    resolve: &SectionResolver,
) -> ContactSolution {
    let mut best = inside.clone();
    let (base_in, base_out) = (inside.center, outside.center);
    let (mut lo, mut hi) = (0.0_f64, 1.0_f64);
    for _ in 0..48 {
        let s = 0.5 * (lo + hi);
        let target = lerp(base_in, base_out, s);
        match resolve(target, &best) {
            Some(sol) if inside_both(&sol, poly_a, poly_b) => {
                best = sol;
                lo = s;
            }
            _ => hi = s,
        }
    }
    best
}

/// March past the last inside section until the trim border is
/// crossed, then bisect back onto it. Falls back to the last inside
/// section when the solver cannot advance (the endpoint gate declines
/// honestly downstream).
fn extend_to_boundary(
    from: &ContactSolution,
    outward: Option<Vec3>,
    step: f64,
    poly_a: &[(f64, f64)],
    poly_b: &[(f64, f64)],
    resolve: &SectionResolver,
) -> ContactSolution {
    let Some(dir) = outward else {
        return from.clone();
    };
    let mut cur = from.clone();
    let mut probe = step;
    for _ in 0..40 {
        let target = cur.center + dir * probe;
        match resolve(target, &cur) {
            Some(sol) if inside_both(&sol, poly_a, poly_b) => {
                if sol.center.distance_to(&cur.center) < probe * 1.0e-3 {
                    break;
                }
                cur = sol;
            }
            Some(sol) => {
                return bisect_to_boundary(&cur, &sol, poly_a, poly_b, resolve);
            }
            None => {
                probe *= 0.5;
                if probe < 1.0e-9 {
                    break;
                }
            }
        }
    }
    cur
}

fn clip_and_resample(
    chain: &[ContactSolution],
    poly_a: &[(f64, f64)],
    poly_b: &[(f64, f64)],
    radius: f64,
    resolve: &SectionResolver,
    tolerance: Tolerance,
) -> Result<Vec<ContactSolution>, BlendError> {
    // Longest contiguous run of sections inside both trims.
    let flags: Vec<bool> = chain
        .iter()
        .map(|s| inside_both(s, poly_a, poly_b))
        .collect();
    let (mut best_start, mut best_len) = (0_usize, 0_usize);
    let mut run_start = 0_usize;
    let mut run_len = 0_usize;
    for (i, &f) in flags.iter().enumerate() {
        if f {
            if run_len == 0 {
                run_start = i;
            }
            run_len += 1;
            if run_len > best_len {
                best_start = run_start;
                best_len = run_len;
            }
        } else {
            run_len = 0;
        }
    }
    if best_len == 0 {
        return Err(decline(
            "trim clip",
            "no marched section lies inside both wall trims".into(),
        ));
    }
    let (i0, i1) = (best_start, best_start + best_len - 1);
    // Boundary sections: bisect against an outside neighbour when the
    // march produced one; otherwise the march stopped short of the
    // trim (solver-range grounds, measured 0.15 short on the
    // trapezoid pocket) — extend with shrinking resolve steps until
    // the trim border is crossed, then bisect.
    let step = (radius * 0.5).max(tolerance.linear() * 100.0);
    let outward_start = (chain[i0].center - chain[(i0 + 1).min(chain.len() - 1)].center)
        .normalized()
        .or_else(|| (chain[i0.saturating_sub(1)].center - chain[i0].center).normalized());
    let start_section = if i0 > 0 {
        bisect_to_boundary(&chain[i0], &chain[i0 - 1], poly_a, poly_b, resolve)
    } else {
        extend_to_boundary(&chain[i0], outward_start, step, poly_a, poly_b, resolve)
    };
    let outward_end = (chain[i1].center - chain[i1.saturating_sub(1)].center)
        .normalized()
        .or_else(|| (chain[(i1 + 1).min(chain.len() - 1)].center - chain[i1].center).normalized());
    let end_section = if i1 + 1 < chain.len() {
        bisect_to_boundary(&chain[i1], &chain[i1 + 1], poly_a, poly_b, resolve)
    } else {
        extend_to_boundary(&chain[i1], outward_end, step, poly_a, poly_b, resolve)
    };

    // Uniform resample along the centre polyline.
    let mut sections: Vec<ContactSolution> = Vec::with_capacity(best_len + 2);
    sections.push(start_section);
    sections.extend(chain[i0..=i1].iter().cloned());
    sections.push(end_section);
    let centres: Vec<Point3> = sections.iter().map(|s| s.center).collect();
    let mut cum = vec![0.0_f64; centres.len()];
    for i in 1..centres.len() {
        cum[i] = cum[i - 1] + centres[i].distance_to(&centres[i - 1]);
    }
    let total = cum[centres.len() - 1];
    if total < tolerance.linear().max(1.0e-9) * 10.0 {
        return Err(decline(
            "trim clip",
            format!("clipped strip span is degenerate: length {total:.3e}"),
        ));
    }
    let step = (radius * 0.25).max(tolerance.linear() * 100.0);
    let m = ((total / step).ceil() as usize + 1).clamp(2, 800);
    let mut out: Vec<ContactSolution> = Vec::with_capacity(m);
    out.push(sections[0].clone());
    for k in 1..m - 1 {
        let sk = total * k as f64 / (m - 1) as f64;
        let j = cum.iter().rposition(|&c| c <= sk).unwrap_or(0);
        let j = j.min(centres.len() - 2);
        let seg = (cum[j + 1] - cum[j]).max(1.0e-30);
        let target = lerp(centres[j], centres[j + 1], (sk - cum[j]) / seg);
        let sol = resolve(target, &sections[j]).ok_or_else(|| {
            decline(
                "resample",
                format!("section solve failed at station {k}/{m}"),
            )
        })?;
        out.push(sol);
    }
    out.push(sections[sections.len() - 1].clone());
    Ok(out)
}
