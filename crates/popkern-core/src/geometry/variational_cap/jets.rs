//! Neighbour jets at rim samples: the foot point, the oriented unit
//! normal and the oriented second fundamental form, from exact span
//! polynomials (the certifier's jet source) rather than finite
//! differences.

use crate::geometry::certify::interval::Iv;
use crate::geometry::certify::surface_jet::SurfaceJetSource;
use crate::geometry::surface::SurfaceGeometry;
use crate::geometry::uv_projection::{project_point_to_uv, project_point_to_uv_near};
use crate::geometry::{Point3, Vec3};

/// The neighbour's local geometry at one rim sample.
#[derive(Debug, Clone, Copy)]
pub(super) struct NeighbourJet {
    /// Foot point on the neighbour.
    pub foot: Point3,
    /// Unit normal, oriented as the neighbour face's outward normal.
    pub normal: Vec3,
    /// The oriented second fundamental form as a symmetric 3×3 form on
    /// tangent vectors: `II(a, b) = aᵀ Q b` for `a, b` in the tangent
    /// plane.
    pub form: [[f64; 3]; 3],
}

impl NeighbourJet {
    /// `II(a, b)`.
    pub fn ii(&self, a: Vec3, b: Vec3) -> f64 {
        let q = &self.form;
        let qb = [
            q[0][0] * b.x + q[0][1] * b.y + q[0][2] * b.z,
            q[1][0] * b.x + q[1][1] * b.y + q[1][2] * b.z,
            q[2][0] * b.x + q[2][1] * b.y + q[2][2] * b.z,
        ];
        a.x * qb[0] + a.y * qb[1] + a.z * qb[2]
    }

    /// `Q a` (the gradient of `II(a, a)` is `2 Q a`).
    pub fn apply(&self, a: Vec3) -> Vec3 {
        let q = &self.form;
        Vec3::new(
            q[0][0] * a.x + q[0][1] * a.y + q[0][2] * a.z,
            q[1][0] * a.x + q[1][1] * a.y + q[1][2] * a.z,
            q[2][0] * a.x + q[2][1] * a.y + q[2][2] * a.z,
        )
    }
}

/// Walks a rim, projecting consecutive samples with continuation.
pub(super) struct NeighbourSampler<'a> {
    surface: &'a SurfaceGeometry,
    source: SurfaceJetSource,
    reversed: bool,
    hint: Option<(f64, f64)>,
}

impl<'a> NeighbourSampler<'a> {
    /// A sampler over `surface`; `reversed` flips the natural normal.
    pub fn new(surface: &'a SurfaceGeometry, reversed: bool) -> Result<Self, String> {
        let source = SurfaceJetSource::new(surface).map_err(|e| format!("{e:?}"))?;
        Ok(NeighbourSampler {
            surface,
            source,
            reversed,
            hint: None,
        })
    }

    /// The jet at the neighbour point nearest `p`.
    pub fn jet_at(&mut self, p: &Point3) -> Result<NeighbourJet, String> {
        let uv = match self.hint {
            Some((u, v)) => project_point_to_uv_near(self.surface, p, u, v)
                .or_else(|| project_point_to_uv(self.surface, p)),
            None => project_point_to_uv(self.surface, p),
        }
        .ok_or("projection failed")?;
        self.hint = Some(uv);
        let (iu, iv) = (Iv::point(uv.0), Iv::point(uv.1));
        let span = self
            .source
            .spans_touching(iu, iv)
            .map_err(|e| format!("{e:?}"))?
            .first()
            .copied()
            .ok_or("no span at the foot point")?;
        let jet = self
            .source
            .jet_in(iu, iv, span)
            .map_err(|e| format!("{e:?}"))?;
        let s = jet.partial(0, 0).mid();
        let su = jet.partial(1, 0).mid();
        let sv = jet.partial(0, 1).mid();
        let suu = jet.partial(2, 0).mid();
        let suv = jet.partial(1, 1).mid();
        let svv = jet.partial(0, 2).mid();
        let mut n = su.cross(&sv).normalized().ok_or("singular chart")?;
        if self.reversed {
            n = n * -1.0;
        }
        let form = form_from_partials(su, sv, suu, suv, svv, n)
            .ok_or("degenerate first fundamental form")?;
        Ok(NeighbourJet {
            foot: Point3::new(s.x, s.y, s.z),
            normal: n,
            form,
        })
    }
}

/// The second fundamental form w.r.t. unit normal `n` as a symmetric
/// 3×3 form on tangent vectors (`II(a, b) = aᵀ Q b`), from the chart
/// partials; `None` on a degenerate chart.
pub(super) fn form_from_partials(
    su: Vec3,
    sv: Vec3,
    suu: Vec3,
    suv: Vec3,
    svv: Vec3,
    n: Vec3,
) -> Option<[[f64; 3]; 3]> {
    let (e, f, g) = (su.dot(&su), su.dot(&sv), sv.dot(&sv));
    let det = e * g - f * f;
    if det.is_nan() || det <= 1e-24 * e * g {
        return None;
    }
    // P = (AᵀA)⁻¹ Aᵀ maps a tangent vector to chart coordinates.
    let rows = [
        su * (g / det) + sv * (-f / det),
        su * (-f / det) + sv * (e / det),
    ];
    let b = [[suu.dot(&n), suv.dot(&n)], [suv.dot(&n), svv.dot(&n)]];
    let comp = |v: Vec3, k: usize| match k {
        0 => v.x,
        1 => v.y,
        _ => v.z,
    };
    let mut form = [[0.0; 3]; 3];
    for (i, row) in form.iter_mut().enumerate() {
        for (j, slot) in row.iter_mut().enumerate() {
            let mut acc = 0.0;
            for (a, ra) in rows.iter().enumerate() {
                for (c, rc) in rows.iter().enumerate() {
                    acc += comp(*ra, i) * b[a][c] * comp(*rc, j);
                }
            }
            *slot = acc;
        }
    }
    Some(form)
}

/// Spectral norm of the difference of two forms restricted to the
/// tangent plane spanned by the orthonormal pair `(e1, e2)`.
pub(super) fn form_jump(a: &[[f64; 3]; 3], b: &[[f64; 3]; 3], e1: Vec3, e2: Vec3) -> f64 {
    let ev = |q: &[[f64; 3]; 3], x: Vec3, y: Vec3| {
        let qy = [
            q[0][0] * y.x + q[0][1] * y.y + q[0][2] * y.z,
            q[1][0] * y.x + q[1][1] * y.y + q[1][2] * y.z,
            q[2][0] * y.x + q[2][1] * y.y + q[2][2] * y.z,
        ];
        x.x * qy[0] + x.y * qy[1] + x.z * qy[2]
    };
    let d11 = ev(a, e1, e1) - ev(b, e1, e1);
    let d12 = ev(a, e1, e2) - ev(b, e1, e2);
    let d22 = ev(a, e2, e2) - ev(b, e2, e2);
    0.5 * (d11 + d22).abs() + (0.25 * (d11 - d22) * (d11 - d22) + d12 * d12).sqrt()
}
