use super::*;

/// A bicubic patch over `[0, 1]²` with points `f(u, v)` at its
/// control stations (a polynomial surface of degree at most one in each
/// parameter is reproduced exactly by them).
fn patch(f: impl Fn(f64, f64) -> [f64; 3], w: impl Fn(usize, usize) -> f64) -> Patch {
    let mut net = Vec::new();
    for a in 0..4 {
        for b in 0..4 {
            let p = f(a as f64 / 3.0, b as f64 / 3.0);
            let wt = w(a, b);
            net.push([p[0] * wt, p[1] * wt, p[2] * wt, wt]);
        }
    }
    Patch::new(4, 4, net).expect("patch")
}

/// A tilted plane's normal is never across an axis it is not
/// perpendicular to, and always across one lying in the plane.
#[test]
fn a_plane_patch_is_decided_against_an_axis() {
    let tilted = patch(|u, v| [u, v, 0.5 * u], |_, _| 1.0);
    assert!(normals_never_across(
        std::slice::from_ref(&tilted),
        Vec3::new(0.0, 0.0, 1.0)
    ));
    // The axis (1, 0, 0.5) lies in the plane: every normal is across it.
    let in_plane = Vec3::new(1.0, 0.0, 0.5).normalized().expect("axis");
    assert!(!normals_never_across(&[tilted], in_plane));
}

/// A cone-like flank (normals leaning off the axis everywhere) is never
/// across its axis; a cylinder-like wall (radial normals) always is; a
/// saddle whose normal turns through the axis's perpendicular is not
/// decided, nor is a rational patch.
#[test]
fn flanks_walls_saddles_and_rational_patches() {
    let z = Vec3::new(0.0, 0.0, 1.0);
    // A ruled strip from radius 1 at z 0 to radius 2 at z 1 over a
    // sixth of a turn, its rulings at 45° to the axis (a cone's flank,
    // drawn by its chords).
    let flank = patch(
        |u, v| {
            let t = u * std::f64::consts::FRAC_PI_3;
            let r = 1.0 + v;
            [r * (1.0 - t * t / 2.0), r * (t - t * t * t / 6.0), v]
        },
        |_, _| 1.0,
    );
    assert!(normals_never_across(std::slice::from_ref(&flank), z));
    let wall = patch(
        |u, v| {
            let t = u * std::f64::consts::FRAC_PI_3;
            [1.0 - t * t / 2.0, t - t * t * t / 6.0, v]
        },
        |_, _| 1.0,
    );
    assert!(!normals_never_across(&[wall], z));
    let saddle = patch(
        |u, v| [u - 0.5, v - 0.5, (u - 0.5) * (v - 0.5) * 4.0],
        |_, _| 1.0,
    );
    assert!(!normals_never_across(&[saddle], Vec3::new(1.0, 0.0, 0.0)));
    let rational = patch(|u, v| [u, v, 0.5 * u], |a, _| 1.0 + 0.1 * a as f64);
    assert!(!normals_never_across(&[rational], z));
}
