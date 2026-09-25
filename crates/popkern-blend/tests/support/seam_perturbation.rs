//! PK-498 negative-control perturbations for bi-nonic cap patches.
//!
//! A continuity gate is only evidence if it fails on the defect it
//! claims to gate. These helpers perturb one side of a bi-nonic
//! (degree-9 Bézier) cap piece so that only one order changes there:
//! - **Row 2** (along the patch normal) changes the transverse
//!   curvature. Position and normal along the side are exact, and so
//!   is the along-seam normal curvature (Meusnier).
//! - **Row 1** changes the cross tangent, so the normal tilts while
//!   position stays exact.
//!
//! The moved control points are the interior along-side indices 3..=6
//! only. That leaves the other three sides' position, tangent and
//! curvature rows untouched, and never touches the dual-jet
//! corner-class entries (1,2)/(2,1)/(2,2) and their mirrors, whose
//! evaluation is replaced by claims.
//!
//! Callers self-check with [`side_change`] before trusting a verdict.
#![allow(dead_code)]

use popkern_core::geometry::Vec3;
use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

/// A side of the unit parameter square.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PatchSide {
    /// u = 0.
    U0,
    /// u = 1.
    U1,
    /// v = 0.
    V0,
    /// v = 1.
    V1,
}

impl PatchSide {
    /// The parameter point at fraction `t` along the side.
    pub fn at(self, t: f64) -> (f64, f64) {
        match self {
            PatchSide::U0 => (0.0, t),
            PatchSide::U1 => (1.0, t),
            PatchSide::V0 => (t, 0.0),
            PatchSide::V1 => (t, 1.0),
        }
    }

    /// Grid index `(i, j)` (i ↔ u, j ↔ v) at depth `row` from this side
    /// and position `along` on it.
    fn index(self, row: usize, along: usize) -> (usize, usize) {
        match self {
            PatchSide::U0 => (row, along),
            PatchSide::U1 => (9 - row, along),
            PatchSide::V0 => (along, row),
            PatchSide::V1 => (along, 9 - row),
        }
    }
}

fn patch_normal(patch: &GregoryPatch, (u, v): (f64, f64)) -> Option<Vec3> {
    patch
        .partial_u(u, v)
        .cross(&patch.partial_v(u, v))
        .normalized()
}

/// Mean cross leg `|row 1 − row 0|` along `side`: the scale offsets
/// are chosen against.
pub fn side_leg(patch: &GregoryPatch, side: PatchSide) -> Option<f64> {
    let grid = patch.binonic_grid()?;
    let legs: Vec<f64> = (0..10)
        .map(|along| {
            let (i0, j0) = side.index(0, along);
            let (i1, j1) = side.index(1, along);
            grid[i0][j0].distance_to(&grid[i1][j1])
        })
        .collect();
    Some(legs.iter().sum::<f64>() / legs.len() as f64)
}

/// Move the interior control points at depth `row` from `side` by
/// `offset` along the patch normal at the side's midpoint.
pub fn perturb_side_row(
    patch: &mut GregoryPatch,
    side: PatchSide,
    row: usize,
    offset: f64,
) -> Result<(), String> {
    let normal = patch_normal(patch, side.at(0.5)).ok_or("degenerate side normal")?;
    let grid = patch
        .binonic_grid_mut()
        .ok_or("not a bi-nonic-family patch")?;
    for along in 3..=6 {
        let (i, j) = side.index(row, along);
        grid[i][j] = grid[i][j] + normal * offset;
    }
    Ok(())
}

/// Worst position and normal change along `side` between two versions
/// of a patch: the self-check that a perturbation kept the orders it
/// claims to keep.
pub fn side_change(before: &GregoryPatch, after: &GregoryPatch, side: PatchSide) -> (f64, f64) {
    let (mut position, mut normal) = (0.0_f64, 0.0_f64);
    for k in 0..=64 {
        let uv = side.at(0.05 + 0.9 * k as f64 / 64.0);
        position = position.max(
            before
                .point_at(uv.0, uv.1)
                .distance_to(&after.point_at(uv.0, uv.1)),
        );
        if let (Some(a), Some(b)) = (patch_normal(before, uv), patch_normal(after, uv)) {
            normal = normal.max(a.cross(&b).length());
        }
    }
    (position, normal)
}

/// The side of `face`'s patch that `edge` runs along, read from the
/// face's pcurve for that edge (both ends on one side of the unit
/// square).
pub fn seam_side_of(brep: &BRep, face: FaceHandle, edge: EdgeHandle) -> Option<PatchSide> {
    let he = brep.half_edges.iter().find_map(|(h, d)| {
        (d.edge == edge && brep.loops.get(d.loop_).map(|l| l.face) == Some(face)).then_some(h)
    })?;
    let pcurve = brep.half_edges.get(he)?.pcurve.as_ref()?;
    let (a, b) = (
        pcurve.point_at(pcurve.t_start),
        pcurve.point_at(pcurve.t_end),
    );
    let near = |x: f64, target: f64| (x - target).abs() < 1e-6;
    [
        (PatchSide::U0, near(a.x, 0.0) && near(b.x, 0.0)),
        (PatchSide::U1, near(a.x, 1.0) && near(b.x, 1.0)),
        (PatchSide::V0, near(a.y, 0.0) && near(b.y, 0.0)),
        (PatchSide::V1, near(a.y, 1.0) && near(b.y, 1.0)),
    ]
    .into_iter()
    .find_map(|(side, on)| on.then_some(side))
}

/// The Gregory patch on `face`.
pub fn patch_of_face(brep: &BRep, face: FaceHandle) -> Option<GregoryPatch> {
    match &brep.faces.get(face)?.surface.as_ref()?.geometry {
        SurfaceGeometry::Gregory(patch) => Some(patch.clone()),
        _ => None,
    }
}

/// `brep` with `face`'s surface replaced by `patch`.
pub fn with_patch(brep: &BRep, face: FaceHandle, patch: GregoryPatch) -> BRep {
    let mut out = brep.clone();
    if let Some(f) = out.faces.get_mut(face)
        && let Some(on_face) = f.surface.as_mut()
    {
        on_face.geometry = SurfaceGeometry::Gregory(patch);
    }
    out
}

/// A seam shared by two of `faces` along which the first face's patch
/// has a whole side, with that face and side.
pub fn cap_cap_seam(
    brep: &BRep,
    faces: &[FaceHandle],
) -> Option<(EdgeHandle, FaceHandle, PatchSide)> {
    for (edge, data) in brep.edges.iter() {
        let he = brep.half_edges.get(data.half_edge)?;
        let a = brep.loops.get(he.loop_)?.face;
        let b = brep.loops.get(brep.half_edges.get(he.twin)?.loop_)?.face;
        if a == b || !faces.contains(&a) || !faces.contains(&b) {
            continue;
        }
        for face in [a, b] {
            let is_binonic = patch_of_face(brep, face).is_some_and(|p| p.binonic_grid().is_some());
            if let (true, Some(side)) = (is_binonic, seam_side_of(brep, face, edge)) {
                return Some((edge, face, side));
            }
        }
    }
    None
}

/// A ring of `n` bi-nonic patches on the paraboloid `z = x² + y²`, in
/// the `measure_ring` layout: patch k's u = 1 side is patch k + 1's
/// v = 0 side reversed, and every patch's (1, 0) corner is the centre.
/// One smooth surface cut into quads, so every ring seam is exact to
/// all orders: the positive fixture a ring negative control needs.
/// Patch j is the quad (q_j, O, q_{j+1}, e_j) mapped bilinearly and
/// lifted to the paraboloid. That is degree (2, 2), so the degree-9
/// grid interpolating it at 10 × 10 Chebyshev nodes reproduces it
/// exactly.
pub fn paraboloid_ring(n: usize) -> Vec<GregoryPatch> {
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
    let radial = |angle: f64, radius: f64| (radius * angle.cos(), radius * angle.sin());
    let step = std::f64::consts::TAU / n as f64;
    let nodes: Vec<f64> = (0..10)
        .map(|a| 0.5 - 0.5 * ((2 * a + 1) as f64 * std::f64::consts::PI / 20.0).cos())
        .collect();
    let inverse = bernstein_collocation_inverse(&nodes);
    (0..n)
        .map(|j| {
            let q = radial(j as f64 * step, 0.5);
            let q_next = radial((j + 1) as f64 * step, 0.5);
            let e = radial((j as f64 + 0.5) * step, 0.9);
            let surface = |u: f64, v: f64| {
                let x = (1.0 - u) * (1.0 - v) * q.0 + u * v * q_next.0 + (1.0 - u) * v * e.0;
                let y = (1.0 - u) * (1.0 - v) * q.1 + u * v * q_next.1 + (1.0 - u) * v * e.1;
                Point3::new(x, y, x * x + y * y)
            };
            let mut grid = [[Point3::new(0.0, 0.0, 0.0); 10]; 10];
            for i in 0..10 {
                for jj in 0..10 {
                    let (mut x, mut y, mut z) = (0.0, 0.0, 0.0);
                    for a in 0..10 {
                        for b in 0..10 {
                            let w = inverse[i][a] * inverse[jj][b];
                            let p = surface(nodes[a], nodes[b]);
                            x += w * p.x;
                            y += w * p.y;
                            z += w * p.z;
                        }
                    }
                    grid[i][jj] = Point3::new(x, y, z);
                }
            }
            GregoryPatch::RectangularBiNonic(Box::new(RectangularGregoryBiNonic { grid }))
        })
        .collect()
}

/// Inverse of the degree-9 Bernstein collocation matrix `M[a][i] =
/// B_i(nodes[a])`, returned as `inv[i][a]`, so control point i is
/// `Σ_a inv[i][a] · value(nodes[a])`.
fn bernstein_collocation_inverse(nodes: &[f64]) -> [[f64; 10]; 10] {
    let binomial =
        |k: usize| -> f64 { (0..k).fold(1.0, |acc, m| acc * (9 - m) as f64 / (m + 1) as f64) };
    let mut m = [[0.0_f64; 20]; 10];
    for (a, &t) in nodes.iter().enumerate() {
        for (i, entry) in m[a].iter_mut().take(10).enumerate() {
            *entry = binomial(i) * t.powi(i as i32) * (1.0 - t).powi(9 - i as i32);
        }
        m[a][10 + a] = 1.0;
    }
    // Gauss–Jordan with partial pivoting on [M | I].
    for col in 0..10 {
        let pivot = (col..10)
            .max_by(|&x, &y| m[x][col].abs().total_cmp(&m[y][col].abs()))
            .unwrap_or(col);
        m.swap(col, pivot);
        let scale = m[col][col];
        for entry in m[col].iter_mut() {
            *entry /= scale;
        }
        let pivot_row = m[col];
        for (row, r) in m.iter_mut().enumerate() {
            if row != col {
                let factor = r[col];
                for (x, p) in r.iter_mut().zip(pivot_row.iter()) {
                    *x -= factor * p;
                }
            }
        }
    }
    // [I | M⁻¹]: M⁻¹[i][a] maps node values to control points.
    let mut inv = [[0.0_f64; 10]; 10];
    for (row, source) in inv.iter_mut().zip(m.iter()) {
        row.copy_from_slice(&source[10..20]);
    }
    inv
}
