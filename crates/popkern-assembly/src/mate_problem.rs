//! The mate solve's least-squares problem and each mate's residual
//! rows (split from `mate_solve.rs` under the 500-line budget).

use super::*;

pub(crate) struct Problem<'a> {
    pub(crate) resolved: &'a [(Mate, Resolved, Resolved, Option<Resolved>)],
    pub(crate) free: &'a [usize],
    pub(crate) base: &'a [Placement],
    pub(crate) size: f64,
}

impl Problem<'_> {
    pub(crate) fn placement(&self, x: &[f64], k: usize, root: usize) -> Placement {
        let base = self.base[root];
        let w = Vec3::new(x[6 * k], x[6 * k + 1], x[6 * k + 2]);
        let t = Vec3::new(x[6 * k + 3], x[6 * k + 4], x[6 * k + 5]) * self.size;
        let angle = w.length();
        let mut out = *base.transform();
        if angle > 0.0
            && let Ok(turn) =
                Transform3::rotation_around_axis(base.origin(), w * (1.0 / angle), angle)
        {
            out = turn.compose(&out);
        }
        Placement::from_raw(Transform3::translation(t.x, t.y, t.z).compose(&out))
    }

    fn world(&self, x: &[f64], r: &Resolved) -> Primitive {
        match self.free.iter().position(|&f| f == r.root) {
            Some(k) => r.local.placed(&self.placement(x, k, r.root)),
            None => r.local.placed(&self.base[r.root]),
        }
    }

    /// Row ranges per mate (same layout as [`Problem::residual`]).
    pub(crate) fn rows(&self, x: &[f64]) -> Result<Vec<std::ops::Range<usize>>, AssemblyError> {
        let mut out = Vec::new();
        let mut at = 0;
        for (mate, a, b, about) in self.resolved {
            let n = mate_rows(
                mate,
                &self.world(x, a),
                &self.world(x, b),
                about.as_ref().map(|r| self.world(x, r)),
                self.size,
            )?
            .len();
            out.push(at..at + n);
            at += n;
        }
        Ok(out)
    }

    pub(crate) fn residual(&self, x: &[f64]) -> Result<Vec<f64>, AssemblyError> {
        let mut out = Vec::new();
        for (mate, a, b, about) in self.resolved {
            out.extend(mate_rows(
                mate,
                &self.world(x, a),
                &self.world(x, b),
                about.as_ref().map(|r| self.world(x, r)),
                self.size,
            )?);
        }
        Ok(out)
    }

    pub(crate) fn jacobian(&self, x: &[f64]) -> Result<Vec<Vec<f64>>, AssemblyError> {
        let h = 1e-7;
        let mut columns = Vec::with_capacity(x.len());
        for j in 0..x.len() {
            let mut plus = x.to_vec();
            let mut minus = x.to_vec();
            plus[j] += h;
            minus[j] -= h;
            let (rp, rm) = (self.residual(&plus)?, self.residual(&minus)?);
            columns.push(
                rp.iter()
                    .zip(&rm)
                    .map(|(p, m)| (p - m) / (2.0 * h))
                    .collect::<Vec<f64>>(),
            );
        }
        let m = columns.first().map_or(0, Vec::len);
        Ok((0..m)
            .map(|i| columns.iter().map(|c| c[i]).collect())
            .collect())
    }
}

/// The residual rows of one mate between world primitives.
pub(crate) fn mate_rows(
    mate: &Mate,
    a: &Primitive,
    b: &Primitive,
    about: Option<Primitive>,
    size: f64,
) -> Result<Vec<f64>, AssemblyError> {
    let refuse = |detail: &str| AssemblyError::Mate {
        mate: mate.describe(),
        detail: detail.to_string(),
    };
    let v3 = |v: Vec3| [v.x, v.y, v.z];
    let lin = |v: Vec3| [v.x / size, v.y / size, v.z / size];
    Ok(match (&mate.kind, a, b) {
        (MateKind::Coincident { aligned }, Primitive::Plane(pa, na), Primitive::Plane(pb, nb))
        | (
            MateKind::Distance { aligned, .. },
            Primitive::Plane(pa, na),
            Primitive::Plane(pb, nb),
        ) => {
            let s = if *aligned { -1.0 } else { 1.0 };
            let offset = match mate.kind {
                MateKind::Distance { distance, .. } => distance,
                _ => 0.0,
            };
            let mut rows = v3(*na + *nb * s).to_vec();
            rows.push(((*pb - *pa).dot(na) - offset) / size);
            rows
        }
        (
            MateKind::Coincident { .. } | MateKind::Concentric,
            Primitive::Axis(pa, da),
            Primitive::Axis(pb, db),
        ) => {
            let mut rows = v3(da.cross(db)).to_vec();
            rows.extend(lin((*pb - *pa).cross(da)));
            rows
        }
        (MateKind::AxialOffset { distance }, Primitive::Axis(pa, da), Primitive::Axis(pb, _)) => {
            vec![((*pb - *pa).dot(da) - distance) / size]
        }
        (MateKind::AxialOffset { .. }, _, _) => {
            return Err(refuse("an axial offset needs two axes"));
        }
        (MateKind::Coincident { .. }, Primitive::Point(pa), Primitive::Point(pb)) => {
            lin(*pb - *pa).to_vec()
        }
        (MateKind::Coincident { .. }, Primitive::Point(p), Primitive::Plane(o, n))
        | (MateKind::Coincident { .. }, Primitive::Plane(o, n), Primitive::Point(p)) => {
            vec![(*p - *o).dot(n) / size]
        }
        (MateKind::Distance { distance, .. }, Primitive::Point(p), Primitive::Plane(o, n))
        | (MateKind::Distance { distance, .. }, Primitive::Plane(o, n), Primitive::Point(p)) => {
            vec![((*p - *o).dot(n) - distance) / size]
        }
        (
            MateKind::Coincident { .. } | MateKind::Concentric,
            Primitive::Point(p),
            Primitive::Axis(o, d),
        )
        | (
            MateKind::Coincident { .. } | MateKind::Concentric,
            Primitive::Axis(o, d),
            Primitive::Point(p),
        ) => lin((*p - *o).cross(d)).to_vec(),
        (MateKind::Coincident { .. }, Primitive::Plane(o, n), Primitive::Axis(p, d))
        | (MateKind::Coincident { .. }, Primitive::Axis(p, d), Primitive::Plane(o, n)) => {
            vec![d.dot(n), (*p - *o).dot(n) / size]
        }
        (MateKind::Angle { degrees, .. }, _, _) => {
            let (Some(ua), Some(ub)) = (a.direction(), b.direction()) else {
                return Err(refuse(
                    "an angle needs two directions (plane normals or axes)",
                ));
            };
            let theta = degrees.to_radians();
            match about.and_then(|k| k.direction()) {
                Some(k) => {
                    let signed = ua.cross(&ub).dot(&k).patan2(ua.dot(&ub));
                    let mut d = signed - theta;
                    d = (d + std::f64::consts::PI).rem_euclid(std::f64::consts::TAU)
                        - std::f64::consts::PI;
                    vec![d]
                }
                None => vec![ua.dot(&ub) - theta.pcos()],
            }
        }
        (MateKind::Concentric, _, _) => {
            return Err(refuse(
                "concentric needs two axes (cylinder, cone or torus faces, or datum axes) or a point and an axis",
            ));
        }
        _ => return Err(refuse("these two entities cannot take this mate")),
    })
}
