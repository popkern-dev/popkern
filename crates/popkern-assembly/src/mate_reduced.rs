//! Reduced-coordinate rescue for revolute chains (PK-763).
//!
//! The cluster solve's unknowns are a 6-DoF increment per free root, and
//! a revolute pair (concentric + axial offset between two axes) is two
//! constraint mates on those increments. On a serial chain that is a
//! poor parametrisation: a target posing jump moves the flange by
//! dragging twelve-dimensional placements back onto their joint axes,
//! and Gauss–Newton and Levenberg–Marquardt both stall where a
//! joint-space solve walks straight through (six-axis-arm §11 E4: four
//! of the ten misses from home were reachable by damped least squares
//! on the joint angles).
//!
//! The rescue, tried only after both full-coordinate solves miss:
//! - the revolute pairs among the cluster's mates span a tree from the
//!   instances the cluster does not move; each free root it reaches is
//!   parametrised by ONE angle about its joint axis (as posed at the
//!   start), carried by its parent's motion; every other free root keeps
//!   its 6-DoF increment;
//! - Levenberg–Marquardt runs on those unknowns against every mate of
//!   the cluster (the pairs hold by construction), from the current
//!   pose first and then, if that misses, from each joint turned by a
//!   half turn (the chain's other IK branches); among the solutions the
//!   one nearest the current pose (least total joint turn) wins, so a
//!   reachable continuous solution is never traded for a flip, unless
//!   it leaves a released joint past its limits: solutions are ranked by
//!   how far past them they lie first (PK-902);
//! - the solution maps back to the cluster's own increments, so the
//!   caller's residual test, diagnostics and placements are unchanged.
//!
//! Nothing here can serve a placement the caller would not: the result
//! is an increment the cluster problem itself verifies.

use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_math::PortableMath;

use crate::error::AssemblyError;
use crate::mate_linalg::Svd;
use crate::mate_solve::{MATE_TOLERANCE, Problem};
use crate::mates::{MateKind, Primitive};
use crate::placement::Placement;

/// The residual rows of the reduced unknowns.
type ResidualFn<'a> = dyn Fn(&[f64]) -> Result<Vec<f64>, AssemblyError> + 'a;

/// Branch seeds are tried for chains up to this many joints (2^n
/// solves); longer chains take the from-current solve only.
const MAX_BRANCH_JOINTS: usize = 8;

/// A request moving the cluster's mates by less than this fraction of
/// the model size (its largest start residual row; rows are
/// dimensionless) is a continuous edit: it must stay on its IK branch,
/// and is refused where that branch ends rather than served flipped (a
/// 5 mm slider step on the arm turned a joint 112° on another branch).
/// A larger jump, target posing from home, may land on another branch.
pub(crate) const BRANCH_JUMP: f64 = 0.02;

/// One free root reached through a revolute pair.
struct Link {
    /// Index into `problem.free`.
    k: usize,
    /// Its parent's index into `problem.free` (`None`: a root the cluster
    /// does not move).
    parent: Option<usize>,
    /// The joint axis as posed at the start (world).
    origin: Point3,
    direction: Vec3,
}

/// The tree of revolute links (parents before children) and the free
/// roots outside it.
fn joint_tree(problem: &Problem<'_>) -> (Vec<Link>, Vec<usize>) {
    let index = |root: usize| problem.free.iter().position(|&f| f == root);
    // (root a, root b, world axis at the start) per revolute pair.
    let mut pairs = Vec::new();
    for (mate, a, b, _) in problem.resolved {
        if !matches!(mate.kind, MateKind::Concentric) || a.root == b.root {
            continue;
        }
        let seated = problem.resolved.iter().any(|(m, ma, mb, _)| {
            matches!(m.kind, MateKind::AxialOffset { .. })
                && ((ma.root == a.root && mb.root == b.root)
                    || (ma.root == b.root && mb.root == a.root))
        });
        if let (true, Primitive::Axis(o, d)) = (seated, a.local.placed(&problem.base[a.root])) {
            pairs.push((a.root, b.root, o, d));
        }
    }
    let mut placed: Vec<bool> = problem.free.iter().map(|_| false).collect();
    let mut links = Vec::new();
    let mut grew = true;
    while grew {
        grew = false;
        for &(ra, rb, o, d) in &pairs {
            for (from, to) in [(ra, rb), (rb, ra)] {
                let Some(kt) = index(to) else { continue };
                if placed[kt] {
                    continue;
                }
                let parent = match index(from) {
                    None => None,
                    Some(kf) if placed[kf] => Some(kf),
                    Some(_) => continue,
                };
                placed[kt] = true;
                links.push(Link {
                    k: kt,
                    parent,
                    origin: o,
                    direction: d,
                });
                grew = true;
            }
        }
    }
    let rest = (0..problem.free.len()).filter(|&k| !placed[k]).collect();
    (links, rest)
}

/// The placements of the free roots for reduced unknowns `y` (one angle
/// per link, then six per remaining root).
fn placements(
    problem: &Problem<'_>,
    links: &[Link],
    rest: &[usize],
    y: &[f64],
) -> Result<Vec<Placement>, AssemblyError> {
    let n = problem.free.len();
    let mut motion: Vec<Transform3> = vec![Transform3::identity(); n];
    let mut out: Vec<Placement> = problem.free.iter().map(|&r| problem.base[r]).collect();
    for (j, link) in links.iter().enumerate() {
        let turn =
            Transform3::rotation_around_axis(link.origin, link.direction, y[j]).map_err(|e| {
                AssemblyError::Mate {
                    mate: "revolute chain".into(),
                    detail: format!("joint axis: {e}"),
                }
            })?;
        let carried = match link.parent {
            Some(p) => motion[p].compose(&turn),
            None => turn,
        };
        motion[link.k] = carried;
        out[link.k] =
            Placement::from_raw(carried.compose(problem.base[problem.free[link.k]].transform()));
    }
    for (i, &k) in rest.iter().enumerate() {
        let x6 = &y[links.len() + 6 * i..links.len() + 6 * i + 6];
        out[k] = problem.placement(x6, 0, problem.free[k]);
    }
    Ok(out)
}

/// Every cluster mate's rows with the free roots at `placed` (indexed
/// as `problem.free`), the others where they start.
fn residual_at(problem: &Problem<'_>, placed: &[Placement]) -> Result<Vec<f64>, AssemblyError> {
    let world = |r: &crate::mates::Resolved| match problem.free.iter().position(|&f| f == r.root) {
        Some(k) => r.local.placed(&placed[k]),
        None => r.local.placed(&problem.base[r.root]),
    };
    let mut out = Vec::new();
    for (mate, a, b, about) in problem.resolved {
        out.extend(crate::mate_solve::mate_rows(
            mate,
            &world(a),
            &world(b),
            about.as_ref().map(world),
            problem.size,
        )?);
    }
    Ok(out)
}

/// The cluster's own increment (ω, t/size) taking the start placement
/// `base` to `to`: `to = T(t) ∘ R(o, ω) ∘ base` with `o` base's origin.
fn increment(base: &Placement, to: &Placement, size: f64) -> Result<[f64; 6], AssemblyError> {
    let inv = base
        .transform()
        .inverse()
        .map_err(|e| AssemblyError::Mate {
            mate: "revolute chain".into(),
            detail: format!("placement: {e}"),
        })?;
    let q = to.transform().compose(&inv).matrix;
    let w = rotation_log([
        [q[0][0], q[0][1], q[0][2]],
        [q[1][0], q[1][1], q[1][2]],
        [q[2][0], q[2][1], q[2][2]],
    ]);
    let o = base.origin();
    let ro = Vec3::new(
        q[0][0] * o.x + q[0][1] * o.y + q[0][2] * o.z,
        q[1][0] * o.x + q[1][1] * o.y + q[1][2] * o.z,
        q[2][0] * o.x + q[2][1] * o.y + q[2][2] * o.z,
    );
    let t = Vec3::new(q[0][3], q[1][3], q[2][3]) + ro - Vec3::new(o.x, o.y, o.z);
    Ok([w.x, w.y, w.z, t.x / size, t.y / size, t.z / size])
}

/// The rotation vector (axis × angle, angle in [0, π]) of `r`, through
/// its unit quaternion (Shepperd's largest-component extraction, exact
/// to rounding at every angle; the trace formula keeps only √ε near a
/// half turn).
fn rotation_log(r: [[f64; 3]; 3]) -> Vec3 {
    let trace = r[0][0] + r[1][1] + r[2][2];
    let (w, x, y, z) = if trace > r[0][0].max(r[1][1]).max(r[2][2]) {
        let s = (1.0 + trace).sqrt() * 2.0;
        (
            0.25 * s,
            (r[2][1] - r[1][2]) / s,
            (r[0][2] - r[2][0]) / s,
            (r[1][0] - r[0][1]) / s,
        )
    } else if r[0][0] >= r[1][1] && r[0][0] >= r[2][2] {
        let s = (1.0 + r[0][0] - r[1][1] - r[2][2]).sqrt() * 2.0;
        (
            (r[2][1] - r[1][2]) / s,
            0.25 * s,
            (r[0][1] + r[1][0]) / s,
            (r[0][2] + r[2][0]) / s,
        )
    } else if r[1][1] >= r[2][2] {
        let s = (1.0 + r[1][1] - r[0][0] - r[2][2]).sqrt() * 2.0;
        (
            (r[0][2] - r[2][0]) / s,
            (r[0][1] + r[1][0]) / s,
            0.25 * s,
            (r[1][2] + r[2][1]) / s,
        )
    } else {
        let s = (1.0 + r[2][2] - r[0][0] - r[1][1]).sqrt() * 2.0;
        (
            (r[1][0] - r[0][1]) / s,
            (r[0][2] + r[2][0]) / s,
            (r[1][2] + r[2][1]) / s,
            0.25 * s,
        )
    };
    // The short way round: w >= 0.
    let (w, v) = if w < 0.0 {
        (-w, Vec3::new(-x, -y, -z))
    } else {
        (w, Vec3::new(x, y, z))
    };
    let sin_half = v.length();
    if sin_half < 1e-300 {
        return Vec3::new(0.0, 0.0, 0.0);
    }
    let angle = 2.0 * sin_half.patan2(w);
    v * (angle / sin_half)
}

/// Gauss–Newton polish of the cluster's own increments from `x` (the
/// mapped joint-space solution), a few truncated-SVD steps.
fn polish(problem: &Problem<'_>, mut x: Vec<f64>) -> Result<Vec<f64>, AssemblyError> {
    let max = |v: &[f64]| v.iter().fold(0.0_f64, |m, e| m.max(e.abs()));
    let mut r = problem.residual(&x)?;
    for _ in 0..8 {
        if max(&r) <= MATE_TOLERANCE * 1e-2 {
            break;
        }
        let svd = Svd::new(&problem.jacobian(&x)?, x.len());
        let minus_r: Vec<f64> = r.iter().map(|v| -v).collect();
        let step = svd.solve(&minus_r, crate::mate_solve::RANK_REL);
        let trial: Vec<f64> = x.iter().zip(&step).map(|(a, s)| a + s).collect();
        let tr = problem.residual(&trial)?;
        if max(&tr) >= max(&r) {
            break;
        }
        x = trial;
        r = tr;
    }
    Ok(x)
}

/// Levenberg–Marquardt on `y` (central-differenced Jacobian): the end
/// point and its largest residual row.
fn solve_from(
    residual: &ResidualFn<'_>,
    mut y: Vec<f64>,
) -> Result<(Vec<f64>, f64), AssemblyError> {
    let n = y.len();
    let norm = |v: &[f64]| v.iter().map(|e| e * e).sum::<f64>();
    let max = |v: &[f64]| v.iter().fold(0.0_f64, |m, e| m.max(e.abs()));
    let mut r = residual(&y)?;
    let mut lambda: Option<f64> = None;
    for _ in 0..300 {
        if max(&r) <= MATE_TOLERANCE * 1e-2 {
            break;
        }
        let h = 1e-7;
        let mut columns = Vec::with_capacity(n);
        for j in 0..n {
            let (mut p, mut m) = (y.clone(), y.clone());
            p[j] += h;
            m[j] -= h;
            let (rp, rm) = (residual(&p)?, residual(&m)?);
            columns.push(
                rp.iter()
                    .zip(&rm)
                    .map(|(a, b)| (a - b) / (2.0 * h))
                    .collect::<Vec<_>>(),
            );
        }
        let rows = columns.first().map_or(0, Vec::len);
        let jac: Vec<Vec<f64>> = (0..rows)
            .map(|i| columns.iter().map(|c| c[i]).collect())
            .collect();
        let svd = Svd::new(&jac, n);
        let top = svd.max_sigma();
        if top == 0.0 {
            break;
        }
        let minus_r: Vec<f64> = r.iter().map(|v| -v).collect();
        let mut damping = lambda.unwrap_or(1e-3 * top);
        let mut accepted = false;
        while damping < 1e8 * top {
            let step = svd.damped_solve(&minus_r, damping);
            let trial: Vec<f64> = y.iter().zip(&step).map(|(a, s)| a + s).collect();
            let tr = residual(&trial)?;
            if norm(&tr) < norm(&r) {
                y = trial;
                r = tr;
                accepted = true;
                break;
            }
            damping *= 4.0;
        }
        if !accepted {
            break;
        }
        lambda = Some(damping / 3.0);
    }
    let worst = max(&r);
    Ok((y, worst))
}

/// The cluster's increment from a joint-space solve, when one reaches
/// the tolerance: from the current pose, else (when `allow_branches`:
/// a jump, see [`BRANCH_JUMP`]) the branch seed nearest it. `None` when the cluster has no revolute links or no seed solves.
pub(crate) fn rescue(
    problem: &Problem<'_>,
    allow_branches: bool,
) -> Result<Option<Vec<f64>>, AssemblyError> {
    let (links, rest) = joint_tree(problem);
    if links.is_empty() {
        return Ok(None);
    }
    let unknowns = links.len() + 6 * rest.len();
    let residual = |y: &[f64]| -> Result<Vec<f64>, AssemblyError> {
        let placed = placements(problem, &links, &rest, y)?;
        residual_at(problem, &placed)
    };
    let mut seeds = vec![vec![0.0; unknowns]];
    if allow_branches && links.len() <= MAX_BRANCH_JOINTS {
        for mask in 1..(1usize << links.len()) {
            let mut seed = vec![0.0; unknowns];
            for (j, slot) in seed.iter_mut().take(links.len()).enumerate() {
                if mask & (1 << j) != 0 {
                    *slot = std::f64::consts::PI;
                }
            }
            seeds.push(seed);
        }
    }
    let wrap = |a: f64| {
        (a + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU) - std::f64::consts::PI
    };
    let mut best: Option<((f64, f64), Vec<f64>)> = None;
    for (s, seed) in seeds.into_iter().enumerate() {
        let (y, worst) = solve_from(&residual, seed)?;
        if worst > MATE_TOLERANCE {
            continue;
        }
        let placed = placements(problem, &links, &rest, &y)?;
        let excess = crate::mate_solve::limit_excess(problem.limits, |r| {
            match problem.free.iter().position(|&f| f == r.root) {
                Some(k) => r.local.placed(&placed[k]),
                None => r.local.placed(&problem.base[r.root]),
            }
        });
        let turn: f64 = y.iter().take(links.len()).map(|a| wrap(*a).abs()).sum();
        if best
            .as_ref()
            .is_none_or(|((e, t), _)| (excess, turn) < (*e, *t))
        {
            best = Some(((excess, turn), y));
        }
        // The from-current solution is the continuous one: keep it,
        // unless it leaves a joint past its limits.
        if s == 0 && excess == 0.0 {
            break;
        }
    }
    let Some((_, y)) = best else {
        return Ok(None);
    };
    let placed = placements(problem, &links, &rest, &y)?;
    let mut x = Vec::with_capacity(6 * problem.free.len());
    for (k, &root) in problem.free.iter().enumerate() {
        x.extend(increment(&problem.base[root], &placed[k], problem.size)?);
    }
    let x = polish(problem, x)?;
    Ok(Some(x))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn matrix(t: &Transform3) -> [[f64; 3]; 3] {
        let m = t.matrix;
        [
            [m[0][0], m[0][1], m[0][2]],
            [m[1][0], m[1][1], m[1][2]],
            [m[2][0], m[2][1], m[2][2]],
        ]
    }

    #[test]
    fn rotation_log_is_exact_through_a_half_turn() {
        let axis = Vec3::new(0.3, -0.5, 0.81).normalized().expect("axis");
        for angle in [
            0.0,
            1e-9,
            0.7,
            2.5,
            std::f64::consts::PI - 1e-7,
            std::f64::consts::PI,
        ] {
            let t = Transform3::rotation_around_axis(Point3::new(0.0, 0.0, 0.0), axis, angle)
                .expect("rotation");
            let w = rotation_log(matrix(&t));
            let back = if w.length() > 0.0 {
                Transform3::rotation_around_axis(
                    Point3::new(0.0, 0.0, 0.0),
                    w * (1.0 / w.length()),
                    w.length(),
                )
                .expect("rotation")
            } else {
                Transform3::identity()
            };
            assert!(
                t.approx_eq(&back, 1e-15 * 8.0),
                "angle {angle}: log {w:?} does not round-trip"
            );
        }
    }

    #[test]
    fn increments_round_trip_placements() {
        let base = Placement::from_raw(
            Transform3::rotation_around_axis(
                Point3::new(5.0, -2.0, 1.0),
                Vec3::new(0.0, 0.0, 1.0),
                0.4,
            )
            .expect("rotation")
            .compose(&Transform3::translation(30.0, 10.0, -4.0)),
        );
        let to = Placement::from_raw(
            Transform3::rotation_around_axis(
                Point3::new(-7.0, 3.0, 2.0),
                Vec3::new(1.0, 1.0, 0.0).normalized().expect("axis"),
                3.1,
            )
            .expect("rotation")
            .compose(base.transform()),
        );
        let size = 120.0;
        let x = increment(&base, &to, size).expect("increment");
        let w = Vec3::new(x[0], x[1], x[2]);
        let t = Vec3::new(x[3], x[4], x[5]) * size;
        let turn =
            Transform3::rotation_around_axis(base.origin(), w * (1.0 / w.length()), w.length())
                .expect("rotation");
        let rebuilt =
            Transform3::translation(t.x, t.y, t.z).compose(&turn.compose(base.transform()));
        assert!(
            rebuilt.approx_eq(to.transform(), 1e-12),
            "{rebuilt:?} vs {:?}",
            to.transform()
        );
    }
}
