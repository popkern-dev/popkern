//! NURBS curve and surface export to STEP B-spline entities.

use popkern_core::geometry::Vec3;
use popkern_core::geometry::nurbs_curve::NurbsCurve;
use popkern_core::geometry::nurbs_surface::NurbsSurface;

use crate::step::export::{ExportContext, bool_, enum_, int, list, real, ref_, str_param};
use crate::step::parser::types::{EntityId, Parameter};

impl ExportContext {
    /// Export a NURBS curve as B_SPLINE_CURVE_WITH_KNOTS.
    pub(crate) fn export_nurbs_curve(&mut self, nurbs: &NurbsCurve) -> EntityId {
        let ctrl_pts: Vec<Parameter> = nurbs
            .control_points
            .iter()
            .map(|p| ref_(self.export_point(p)))
            .collect();

        let (knot_values, knot_multiplicities) = compress_knots(&nurbs.knots);

        let knot_vals: Vec<Parameter> = knot_values.iter().map(|&k| real(k)).collect();
        let knot_mults: Vec<Parameter> =
            knot_multiplicities.iter().map(|&m| int(m as i64)).collect();

        let is_rational = nurbs.weights.iter().any(|&w| (w - 1.0).abs() > 1e-15);

        if is_rational {
            let weights: Vec<Parameter> = nurbs.weights.iter().map(|&w| real(w)).collect();
            let id = self.alloc_id();
            self.entities
                .push(crate::step::parser::types::EntityInstance {
                    id,
                    type_name: "B_SPLINE_CURVE_WITH_KNOTS".into(),
                    params: vec![
                        str_param(""),
                        int(nurbs.degree as i64),
                        list(ctrl_pts),
                        enum_("UNSPECIFIED"),
                        bool_(false), // closed
                        bool_(false), // self-intersect
                        list(knot_mults),
                        list(knot_vals),
                        enum_("UNSPECIFIED"),
                    ],
                    complex_parts: vec![("RATIONAL_B_SPLINE_CURVE".into(), vec![list(weights)])],
                });
            id
        } else {
            self.emit(
                "B_SPLINE_CURVE_WITH_KNOTS",
                vec![
                    str_param(""),
                    int(nurbs.degree as i64),
                    list(ctrl_pts),
                    enum_("UNSPECIFIED"),
                    bool_(false), // closed
                    bool_(false), // self-intersect
                    list(knot_mults),
                    list(knot_vals),
                    enum_("UNSPECIFIED"),
                ],
            )
        }
    }

    /// Export a NURBS surface as B_SPLINE_SURFACE_WITH_KNOTS.
    pub(crate) fn export_nurbs_surface(&mut self, nurbs: &NurbsSurface) -> EntityId {
        let mut rows: Vec<Parameter> = Vec::new();
        for i in 0..nurbs.num_u {
            let row: Vec<Parameter> = (0..nurbs.num_v)
                .map(|j| {
                    let idx = i * nurbs.num_v + j;
                    ref_(self.export_point(&nurbs.control_points[idx]))
                })
                .collect();
            rows.push(list(row));
        }

        let (knot_u_vals, knot_u_mults) = compress_knots(&nurbs.knots_u);
        let (knot_v_vals, knot_v_mults) = compress_knots(&nurbs.knots_v);

        let ku_vals: Vec<Parameter> = knot_u_vals.iter().map(|&k| real(k)).collect();
        let ku_mults: Vec<Parameter> = knot_u_mults.iter().map(|&m| int(m as i64)).collect();
        let kv_vals: Vec<Parameter> = knot_v_vals.iter().map(|&k| real(k)).collect();
        let kv_mults: Vec<Parameter> = knot_v_mults.iter().map(|&m| int(m as i64)).collect();

        let is_rational = nurbs.weights.iter().any(|&w| (w - 1.0).abs() > 1e-15);

        if is_rational {
            let mut weight_rows: Vec<Parameter> = Vec::new();
            for i in 0..nurbs.num_u {
                let row: Vec<Parameter> = (0..nurbs.num_v)
                    .map(|j| real(nurbs.weights[i * nurbs.num_v + j]))
                    .collect();
                weight_rows.push(list(row));
            }

            let id = self.alloc_id();
            self.entities
                .push(crate::step::parser::types::EntityInstance {
                    id,
                    type_name: "B_SPLINE_SURFACE_WITH_KNOTS".into(),
                    params: vec![
                        str_param(""),
                        int(nurbs.degree_u as i64),
                        int(nurbs.degree_v as i64),
                        list(rows),
                        enum_("UNSPECIFIED"),
                        bool_(false),
                        bool_(false),
                        bool_(false),
                        list(ku_mults),
                        list(kv_mults),
                        list(ku_vals),
                        list(kv_vals),
                        enum_("UNSPECIFIED"),
                    ],
                    complex_parts: vec![(
                        "RATIONAL_B_SPLINE_SURFACE".into(),
                        vec![list(weight_rows)],
                    )],
                });
            id
        } else {
            self.emit(
                "B_SPLINE_SURFACE_WITH_KNOTS",
                vec![
                    str_param(""),
                    int(nurbs.degree_u as i64),
                    int(nurbs.degree_v as i64),
                    list(rows),
                    enum_("UNSPECIFIED"),
                    bool_(false),
                    bool_(false),
                    bool_(false),
                    list(ku_mults),
                    list(kv_mults),
                    list(ku_vals),
                    list(kv_vals),
                    enum_("UNSPECIFIED"),
                ],
            )
        }
    }
}

/// Compute an arbitrary perpendicular vector to a given normal.
///
/// Returns a unit vector perpendicular to `n`. Used when the STEP entity
/// requires an explicit reference direction that the popkern model does not
/// store (e.g., the X-axis for a cylinder placement).
pub(crate) fn compute_perpendicular(n: &Vec3) -> Vec3 {
    let abs_x = n.x.abs();
    let abs_y = n.y.abs();
    let abs_z = n.z.abs();
    let seed = if abs_x <= abs_y && abs_x <= abs_z {
        Vec3::new(1.0, 0.0, 0.0)
    } else if abs_y <= abs_z {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let perp = n.cross(&seed);
    perp.normalized().unwrap_or(Vec3::new(1.0, 0.0, 0.0))
}

/// Compress a knot vector into unique values and their multiplicities.
pub(crate) fn compress_knots(knots: &[f64]) -> (Vec<f64>, Vec<usize>) {
    let mut values: Vec<f64> = Vec::new();
    let mut multiplicities: Vec<usize> = Vec::new();
    for &k in knots {
        if let Some(last) = values.last()
            && (k - *last).abs() < 1e-15
        {
            if let Some(m) = multiplicities.last_mut() {
                *m += 1;
            }
            continue;
        }
        values.push(k);
        multiplicities.push(1);
    }
    (values, multiplicities)
}

#[cfg(test)]
mod tests {
    use super::*;
    use popkern_core::geometry::Point3;
    use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};

    #[test]
    fn export_nurbs_curve_non_rational() {
        let mut ctx = ExportContext::new();
        let nurbs = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 1.0, 1.0],
        };
        let id = ctx.export_nurbs_curve(&nurbs);
        let entity = ctx.entities.iter().find(|e| e.id == id);
        assert!(entity.is_some());
        assert_eq!(entity.unwrap().type_name, "B_SPLINE_CURVE_WITH_KNOTS");
        assert!(entity.unwrap().complex_parts.is_empty());
    }

    #[test]
    fn export_nurbs_curve_rational() {
        let mut ctx = ExportContext::new();
        let nurbs = NurbsCurve {
            degree: 2,
            knots: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points: vec![
                Point3::new(0.0, 0.0, 0.0),
                Point3::new(0.5, 1.0, 0.0),
                Point3::new(1.0, 0.0, 0.0),
            ],
            weights: vec![1.0, 0.707, 1.0],
        };
        let curve = CurveOnEdge {
            geometry: CurveGeometry::Nurbs(nurbs),
            t_start: 0.0,
            t_end: 1.0,
        };
        let id = ctx.export_curve(&curve);
        let entity = ctx.entities.iter().find(|e| e.id == id);
        assert!(entity.is_some());
        assert_eq!(entity.unwrap().complex_parts.len(), 1);
        assert_eq!(
            entity.unwrap().complex_parts[0].0,
            "RATIONAL_B_SPLINE_CURVE"
        );
    }

    #[test]
    fn compress_knots_basic() {
        let knots = vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0];
        let (vals, mults) = compress_knots(&knots);
        assert_eq!(vals, vec![0.0, 0.5, 1.0]);
        assert_eq!(mults, vec![3, 1, 3]);
    }

    #[test]
    fn compute_perpendicular_z_axis() {
        let n = Vec3::new(0.0, 0.0, 1.0);
        let p = compute_perpendicular(&n);
        assert!((p.length() - 1.0).abs() < 1e-10);
        assert!(n.dot(&p).abs() < 1e-10);
    }

    #[test]
    fn compute_perpendicular_x_axis() {
        let n = Vec3::new(1.0, 0.0, 0.0);
        let p = compute_perpendicular(&n);
        assert!((p.length() - 1.0).abs() < 1e-10);
        assert!(n.dot(&p).abs() < 1e-10);
    }

    #[test]
    fn compute_perpendicular_arbitrary() {
        let n = Vec3::new(1.0, 1.0, 1.0);
        let p = compute_perpendicular(&n);
        assert!((p.length() - 1.0).abs() < 1e-10);
        assert!(n.dot(&p).abs() < 1e-10);
    }
}
