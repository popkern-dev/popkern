use super::*;
use popkern_core::geometry::surface::SurfaceEval;

/// A rational surface of degree 2 x 3 over unequal spans, a double
/// interior knot in u, weights 0.6 to 1.8, its net a waved sheet.
fn wavy(offset: Vec3) -> NurbsSurface {
    let (degree_u, degree_v) = (2, 3);
    let knots_u = vec![0.0, 0.0, 0.0, 0.3, 0.3, 0.55, 1.0, 1.0, 1.0];
    let knots_v = vec![0.0, 0.0, 0.0, 0.0, 0.2, 0.7, 1.0, 1.0, 1.0, 1.0];
    let (num_u, num_v) = (knots_u.len() - degree_u - 1, knots_v.len() - degree_v - 1);
    let mut control_points = Vec::new();
    let mut weights = Vec::new();
    for i in 0..num_u {
        for j in 0..num_v {
            let (x, y) = (i as f64, 0.8 * j as f64);
            let z = 0.4 * ((i * 7 + j * 3) % 5) as f64 - 0.6;
            control_points.push(Point3::new(x, y, z) + offset);
            weights.push(0.6 + 0.2 * ((i * 5 + j * 2) % 7) as f64);
        }
    }
    NurbsSurface {
        degree_u,
        degree_v,
        knots_u,
        knots_v,
        control_points,
        num_u,
        num_v,
        weights,
    }
}

fn whole(s: &NurbsSurface) -> ((f64, f64), (f64, f64)) {
    (
        (s.knots_u[s.degree_u], s.knots_u[s.num_u]),
        (s.knots_v[s.degree_v], s.knots_v[s.num_v]),
    )
}

/// A patch evaluated at local `(a, b)` in `[0, 1]^2` (de Casteljau over
/// its homogeneous net).
fn patch_at(p: &Patch, a: f64, b: f64) -> Point3 {
    let lerp = |row: &mut Vec<[f64; 4]>, t: f64| {
        while row.len() > 1 {
            for m in 0..row.len() - 1 {
                row[m] = std::array::from_fn(|c| (1.0 - t) * row[m][c] + t * row[m + 1][c]);
            }
            row.pop();
        }
        row[0]
    };
    let mut col: Vec<[f64; 4]> = (0..p.nu)
        .map(|i| lerp(&mut (0..p.nv).map(|j| p.net[i * p.nv + j]).collect(), b))
        .collect();
    let h = lerp(&mut col, a);
    Point3::new(h[0] / h[3], h[1] / h[3], h[2] / h[3])
}

/// The nonempty spans of `knots` (degree `p`, `count` points).
fn spans(knots: &[f64], p: usize, count: usize) -> Vec<(f64, f64)> {
    (p..count)
        .filter(|&i| knots[i] < knots[i + 1])
        .map(|i| (knots[i], knots[i + 1]))
        .collect()
}

/// Each span's Bézier patch is the surface over that span, to 1e-12 at a
/// lattice of points (the double knot's spans and the unequal ones too),
/// and its control box holds those points.
#[test]
fn each_span_patch_is_the_surface_over_its_span() {
    let s = wavy(Vec3::new(0.0, 0.0, 0.0));
    let patches = window_patches(&s, &whole(&s)).expect("patches");
    let (su, sv) = (
        spans(&s.knots_u, s.degree_u, s.num_u),
        spans(&s.knots_v, s.degree_v, s.num_v),
    );
    assert_eq!(patches.len(), su.len() * sv.len());
    let surface = SurfaceGeometry::Nurbs(s.clone());
    for (k, p) in patches.iter().enumerate() {
        let ((u0, u1), (v0, v1)) = (su[k / sv.len()], sv[k % sv.len()]);
        for i in 0..=8 {
            for j in 0..=8 {
                let (a, b) = (f64::from(i) / 8.0, f64::from(j) / 8.0);
                let want = surface.point_at(u0 + a * (u1 - u0), v0 + b * (v1 - v0));
                let got = patch_at(p, a, b);
                assert!(
                    got.distance_to(&want) < 1e-12,
                    "span {k} at ({a}, {b}): {got:?} vs {want:?}"
                );
                let inside = |lo: f64, x: f64, hi: f64| x >= lo - 1e-12 && x <= hi + 1e-12;
                assert!(
                    inside(p.bx.min.x, want.x, p.bx.max.x)
                        && inside(p.bx.min.y, want.y, p.bx.max.y)
                        && inside(p.bx.min.z, want.z, p.bx.max.z),
                    "span {k}: {want:?} outside {:?}",
                    p.bx
                );
            }
        }
    }
}

/// A window keeps only the spans it meets.
#[test]
fn a_window_keeps_the_spans_it_meets() {
    let s = wavy(Vec3::new(0.0, 0.0, 0.0));
    // u spans [0, .3] [.3, .55] [.55, 1]; v spans [0, .2] [.2, .7] [.7, 1].
    let some = window_patches(&s, &((0.35, 0.5), (0.1, 0.3))).expect("patches");
    assert_eq!(some.len(), 2);
}

/// Halving a patch along either direction gives the two halves of the
/// surface.
#[test]
fn halves_are_the_patch_halved() {
    let s = wavy(Vec3::new(0.0, 0.0, 0.0));
    let p = &window_patches(&s, &whole(&s)).expect("patches")[4];
    for in_u in [true, false] {
        let (lo, hi) = p.halves(in_u).expect("halves");
        for i in 0..=6 {
            for j in 0..=6 {
                let (a, b) = (f64::from(i) / 6.0, f64::from(j) / 6.0);
                let (whole_lo, whole_hi) = if in_u {
                    (patch_at(p, 0.5 * a, b), patch_at(p, 0.5 + 0.5 * a, b))
                } else {
                    (patch_at(p, a, 0.5 * b), patch_at(p, a, 0.5 + 0.5 * b))
                };
                assert!(patch_at(&lo, a, b).distance_to(&whole_lo) < 1e-12);
                assert!(patch_at(&hi, a, b).distance_to(&whole_hi) < 1e-12);
            }
        }
    }
}

/// Whatever the patch test clears stays farther than the pad: the waved
/// sheet against copies of itself moved over a sweep of offsets, every
/// cleared pair's sampled distance over the pad, and pairs both cleared
/// and kept.
#[test]
fn a_cleared_pair_never_comes_within_the_pad() {
    let a = wavy(Vec3::new(0.0, 0.0, 0.0));
    let pa = window_patches(&a, &whole(&a)).expect("patches");
    let sample = |s: &NurbsSurface| -> Vec<Point3> {
        let g = SurfaceGeometry::Nurbs(s.clone());
        (0..=40)
            .flat_map(|i| (0..=40).map(move |j| (f64::from(i) / 40.0, f64::from(j) / 40.0)))
            .map(|(u, v)| g.point_at(u, v))
            .collect()
    };
    let at = sample(&a);
    let (mut cleared, mut kept) = (0, 0);
    for k in 0..24 {
        let angle = f64::from(k) * 0.7;
        let length = 0.05 + 0.15 * f64::from(k % 6);
        let dir = Vec3::new(angle.cos() * 0.3, angle.sin() * 0.3, 1.0)
            .normalized()
            .expect("dir");
        let b = wavy(dir * length);
        let pb = window_patches(&b, &whole(&b)).expect("patches");
        let pad = 0.02;
        if !patches_apart(&pa, &Other::Patches(&pb), pad) {
            kept += 1;
            continue;
        }
        cleared += 1;
        let near = sample(&b)
            .iter()
            .flat_map(|q| at.iter().map(move |p| p.distance_to(q)))
            .fold(f64::MAX, f64::min);
        assert!(
            near > pad,
            "offset {length} at {angle}: cleared but {near} apart"
        );
    }
    assert!(cleared > 3 && kept > 3, "cleared {cleared}, kept {kept}");
}

/// Two long flanks either side of a crest 0.05 wide, 40 spans each (the
/// cycloidal disc's walls), are cleared at a tolerance's pad; touching
/// at the crest, they are not.
#[test]
fn flanks_either_side_of_a_crest_are_cleared() {
    let flank = |side: f64, gap: f64| {
        let num_v = 43;
        let mut knots_v = vec![0.0; 4];
        knots_v.extend((1..40).map(|k| f64::from(k) / 40.0));
        knots_v.extend(vec![1.0; 4]);
        let mut control_points = Vec::new();
        for i in 0..2 {
            for j in 0..num_v {
                // A helical strip: x along the crest with a rise in z,
                // falling away from the crest across u.
                let x = 10.0 * j as f64 / (num_v - 1) as f64;
                let rise = 0.3 * (0.6 * x).sin();
                let (y, z) = if i == 0 {
                    (side * 0.5 * gap, rise)
                } else {
                    (side * (0.5 * gap + 0.4), rise - 0.7)
                };
                control_points.push(Point3::new(x, y, z));
            }
        }
        NurbsSurface {
            degree_u: 1,
            degree_v: 3,
            knots_u: vec![0.0, 0.0, 1.0, 1.0],
            knots_v,
            control_points,
            num_u: 2,
            num_v,
            weights: vec![1.0; 2 * num_v],
        }
    };
    let pad = popkern_core::tolerance::Tolerance::default().linear_reach();
    for (gap, apart) in [(0.05, true), (0.0, false)] {
        let (a, b) = (flank(1.0, gap), flank(-1.0, gap));
        let (pa, pb) = (
            window_patches(&a, &whole(&a)).expect("patches"),
            window_patches(&b, &whole(&b)).expect("patches"),
        );
        assert_eq!(
            patches_apart(&pa, &Other::Patches(&pb), pad),
            apart,
            "gap {gap}"
        );
    }
}

/// The piece over a window is the surface there, to 1e-12, with the same
/// parameters, and covers the window; a window meeting every span has
/// no piece.
#[test]
fn the_piece_over_a_window_is_the_surface_there() {
    let s = wavy(Vec3::new(0.0, 0.0, 0.0));
    assert!(piece_over(&s, whole(&s)).is_none());
    let window = ((0.35, 0.5), (0.1, 0.3));
    let piece = piece_over(&s, window).expect("piece");
    let ((pu0, pu1), (pv0, pv1)) = whole(&piece);
    assert!(pu0 <= 0.35 && pu1 >= 0.5 && pv0 <= 0.1 && pv1 >= 0.3);
    let (a, b) = (
        SurfaceGeometry::Nurbs(s.clone()),
        SurfaceGeometry::Nurbs(piece),
    );
    for i in 0..=10 {
        for j in 0..=10 {
            let (u, v) = (
                pu0 + (pu1 - pu0) * f64::from(i) / 10.0,
                pv0 + (pv1 - pv0) * f64::from(j) / 10.0,
            );
            assert!(
                a.point_at(u, v).distance_to(&b.point_at(u, v)) < 1e-12,
                "({u}, {v})"
            );
        }
    }
}

/// A window cutting across spans keeps only the part of each span
/// inside it: each patch is the surface over its part to 1e-12 and its
/// box holds those points; and two windows either side of a cut through
/// one span (u 0.3 to 0.55) are apart at the tolerance, where with
/// whole spans both held the span and were never apart (PK-1172: a
/// thread flank either side of a slot).
#[test]
fn a_window_across_a_span_keeps_its_part_of_it() {
    let s = wavy(Vec3::new(0.0, 0.0, 0.0));
    let surface = SurfaceGeometry::Nurbs(s.clone());
    let window = ((0.1, 0.42), (0.25, 0.6));
    let patches = window_patches(&s, &window).expect("patches");
    // u parts [0.1, 0.3] and [0.3, 0.42], v part [0.25, 0.6] of [0.2, 0.7].
    assert_eq!(patches.len(), 2);
    let parts = [((0.1, 0.3), (0.25, 0.6)), ((0.3, 0.42), (0.25, 0.6))];
    for (p, ((u0, u1), (v0, v1))) in patches.iter().zip(parts) {
        for i in 0..=6 {
            for j in 0..=6 {
                let (a, b) = (f64::from(i) / 6.0, f64::from(j) / 6.0);
                let want = surface.point_at(u0 + a * (u1 - u0), v0 + b * (v1 - v0));
                assert!(
                    patch_at(p, a, b).distance_to(&want) < 1e-12,
                    "at ({a}, {b})"
                );
                let inside = |lo: f64, x: f64, hi: f64| x >= lo - 1e-12 && x <= hi + 1e-12;
                assert!(
                    inside(p.bx.min.x, want.x, p.bx.max.x)
                        && inside(p.bx.min.y, want.y, p.bx.max.y)
                        && inside(p.bx.min.z, want.z, p.bx.max.z)
                );
            }
        }
    }
    let left = window_patches(&s, &((0.0, 0.4), (0.0, 1.0))).expect("left");
    let right = window_patches(&s, &((0.45, 1.0), (0.0, 1.0))).expect("right");
    assert!(patches_apart(&left, &Other::Patches(&right), 1e-6));
    // Overlapping windows are not apart.
    let over = window_patches(&s, &((0.35, 1.0), (0.0, 1.0))).expect("over");
    assert!(!patches_apart(&left, &Other::Patches(&over), 1e-6));
}
