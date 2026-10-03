//! B-splines whose knot vectors are implicit in their type (PK-636):
//! `QUASI_UNIFORM_*`, `UNIFORM_*` and `BEZIER_*` curves and surfaces,
//! simple or in the complex rational form. ISO 10303-42 fixes the
//! knots from the degree d and the control-point count n:
//!
//! - quasi-uniform: 0 and n − d with multiplicity d + 1, every integer
//!   between once (clamped, uniform interior);
//! - uniform: −d, −d + 1, …, n, each once;
//! - Bézier: 0 … s with multiplicity d + 1 at the ends and d inside,
//!   s = (n − 1) / d segments.
//!
//! Measured on the CAx-IF corpus: as1-ec, as1-oc, io1-ca and piston
//! write their edges and faces this way, and every one was dropped as
//! an unsupported type (io1-ca meshed with 3424 open edges).

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::surface::SurfaceGeometry;

use crate::step::diagnostics::{DiagnosticLog, Phase};
use crate::step::parser::types::{EntityId, Parameter};

use super::{GeometryCache, expand_knots};

/// The implicit-knot B-spline subtypes.
pub(super) const IMPLICIT_CURVES: [&str; 3] =
    ["QUASI_UNIFORM_CURVE", "UNIFORM_CURVE", "BEZIER_CURVE"];
/// The implicit-knot B-spline surface subtypes.
pub(super) const IMPLICIT_SURFACES: [&str; 3] =
    ["QUASI_UNIFORM_SURFACE", "UNIFORM_SURFACE", "BEZIER_SURFACE"];

/// Knot values and multiplicities of `kind` for degree `d` over `n`
/// control points, or `None` when the count does not fit the type.
pub(super) fn implicit_knots(kind: &str, d: usize, n: usize) -> Option<(Vec<f64>, Vec<i64>)> {
    if d == 0 || n < d + 1 {
        return None;
    }
    let (values, mults): (Vec<f64>, Vec<i64>) = if kind.starts_with("QUASI_UNIFORM") {
        let last = n - d;
        let values = (0..=last).map(|k| k as f64).collect();
        let mults = (0..=last)
            .map(|k| if k == 0 || k == last { d as i64 + 1 } else { 1 })
            .collect();
        (values, mults)
    } else if kind.starts_with("UNIFORM") {
        (
            (0..=n + d).map(|k| k as f64 - d as f64).collect(),
            vec![1; n + d + 1],
        )
    } else if kind.starts_with("BEZIER") {
        if !(n - 1).is_multiple_of(d) {
            return None;
        }
        let s = (n - 1) / d;
        let values = (0..=s).map(|k| k as f64).collect();
        let mults = (0..=s)
            .map(|k| {
                if k == 0 || k == s {
                    d as i64 + 1
                } else {
                    d as i64
                }
            })
            .collect();
        (values, mults)
    } else {
        return None;
    };
    Some((values, mults))
}

impl GeometryCache<'_> {
    /// An implicit-knot B-spline curve, simple (`QUASI_UNIFORM_CURVE(name,
    /// d, points, form, closed, self)`) or complex (`B_SPLINE_CURVE(d,
    /// points, …)` beside the subtype and `RATIONAL_B_SPLINE_CURVE`).
    pub(super) fn resolve_implicit_knot_curve(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<CurveOnEdge> {
        let entity = self.graph.get(id)?.clone();
        let (kind, spline): (&str, Vec<Parameter>) =
            if IMPLICIT_CURVES.contains(&entity.type_name.as_str()) {
                (
                    IMPLICIT_CURVES.iter().find(|k| **k == entity.type_name)?,
                    entity.params.iter().skip(1).cloned().collect(),
                )
            } else {
                let kind = IMPLICIT_CURVES
                    .iter()
                    .find(|k| entity.complex_parts.iter().any(|(n, _)| n == *k))?;
                let (_, p) = entity
                    .complex_parts
                    .iter()
                    .find(|(n, _)| n == "B_SPLINE_CURVE")?;
                (kind, p.clone())
            };
        let Some(Parameter::Int(degree)) = spline.first() else {
            log.error(
                Phase::MapGeometry,
                Some(id),
                format!("{kind}: missing degree"),
            );
            return None;
        };
        let degree = *degree as usize;
        let cp_ids = crate::step::schema::as_ref_list(spline.get(1)?, kind)
            .map_err(|e| log.error(Phase::MapGeometry, Some(id), format!("{e}")))
            .ok()?;
        let mut control_points = Vec::with_capacity(cp_ids.len());
        for cp in cp_ids {
            control_points.push(self.resolve_point(cp, log)?);
        }
        let Some((values, mults)) = implicit_knots(kind, degree, control_points.len()) else {
            log.error(
                Phase::MapGeometry,
                Some(id),
                format!(
                    "{kind}: {} control points do not fit degree {degree}",
                    control_points.len()
                ),
            );
            return None;
        };
        let knots = expand_knots(&values, &mults);
        let weights = super::curves_complex::resolve_rational_bspline_curve_weights(&entity, log)
            .filter(|w| w.len() == control_points.len())
            .unwrap_or_else(|| vec![1.0; control_points.len()]);
        let (t_start, t_end) = (knots[degree], knots[knots.len() - 1 - degree]);
        Some(CurveOnEdge {
            geometry: CurveGeometry::Nurbs(NurbsCurve {
                degree,
                knots,
                control_points,
                weights,
            }),
            t_start,
            t_end,
        })
    }

    /// An implicit-knot B-spline surface, simple or complex: the knots
    /// are appended to the `B_SPLINE_SURFACE` attributes and the
    /// explicit-knot reader takes it from there (weights included).
    pub(super) fn resolve_implicit_knot_surface(
        &mut self,
        id: EntityId,
        log: &mut DiagnosticLog,
    ) -> Option<SurfaceGeometry> {
        let entity = self.graph.get(id)?.clone();
        let (kind, mut params): (&str, Vec<Parameter>) =
            if IMPLICIT_SURFACES.contains(&entity.type_name.as_str()) {
                (
                    IMPLICIT_SURFACES.iter().find(|k| **k == entity.type_name)?,
                    entity.params.clone(),
                )
            } else {
                let kind = IMPLICIT_SURFACES
                    .iter()
                    .find(|k| entity.complex_parts.iter().any(|(n, _)| n == *k))?;
                let (_, p) = entity
                    .complex_parts
                    .iter()
                    .find(|(n, _)| n == "B_SPLINE_SURFACE")?;
                let mut params = vec![Parameter::String(String::new())];
                params.extend(p.iter().cloned());
                (kind, params)
            };
        let (Some(Parameter::Int(du)), Some(Parameter::Int(dv)), Some(Parameter::List(rows))) =
            (params.get(1), params.get(2), params.get(3))
        else {
            log.error(
                Phase::MapGeometry,
                Some(id),
                format!("{kind}: malformed B_SPLINE_SURFACE attributes"),
            );
            return None;
        };
        let (du, dv) = (*du as usize, *dv as usize);
        let nu = rows.len();
        let nv = match rows.first() {
            Some(Parameter::List(r)) => r.len(),
            _ => 0,
        };
        let (Some((ku, mu)), Some((kv, mv))) =
            (implicit_knots(kind, du, nu), implicit_knots(kind, dv, nv))
        else {
            log.error(
                Phase::MapGeometry,
                Some(id),
                format!("{kind}: a {nu}×{nv} grid does not fit degrees ({du}, {dv})"),
            );
            return None;
        };
        let ints = |m: Vec<i64>| Parameter::List(m.into_iter().map(Parameter::Int).collect());
        let reals = |k: Vec<f64>| Parameter::List(k.into_iter().map(Parameter::Real).collect());
        params.truncate(8);
        params.extend([
            ints(mu),
            ints(mv),
            reals(ku),
            reals(kv),
            Parameter::Enum("UNSPECIFIED".into()),
        ]);
        self.resolve_bspline_surface_from(id, &params, log)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn implicit_knots_follow_iso_10303_42() {
        let (v, m) = implicit_knots("QUASI_UNIFORM_CURVE", 3, 6).expect("fits");
        assert_eq!(
            expand_knots(&v, &m),
            vec![0., 0., 0., 0., 1., 2., 3., 3., 3., 3.]
        );
        let (v, m) = implicit_knots("UNIFORM_CURVE", 2, 4).expect("fits");
        assert_eq!(expand_knots(&v, &m), vec![-2., -1., 0., 1., 2., 3., 4.]);
        let (v, m) = implicit_knots("BEZIER_CURVE", 3, 7).expect("fits");
        assert_eq!(
            expand_knots(&v, &m),
            vec![0., 0., 0., 0., 1., 1., 1., 2., 2., 2., 2.]
        );
        assert!(implicit_knots("BEZIER_CURVE", 3, 6).is_none());
        assert!(implicit_knots("QUASI_UNIFORM_CURVE", 3, 3).is_none());
        // Linear quasi-uniform: the io1-ca edges.
        let (v, m) = implicit_knots("QUASI_UNIFORM_CURVE", 1, 2).expect("fits");
        assert_eq!(expand_knots(&v, &m), vec![0., 0., 1., 1.]);
    }
}
