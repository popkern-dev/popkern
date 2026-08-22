use crate::geometry::nurbs_surface::NurbsSurface;
use crate::geometry::{Point3, Vec3};
use serde::{Deserialize, Serialize};

/// Supported ACM 1032 polyhedral spline layout families.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum PolyhedralSplineFamily {
    RegularQuad,
    ExtraordinaryStar { valence: usize },
    Ngon { sides: usize },
    Polar { sectors: usize },
    T0,
    T1,
    T2,
    TwoTrianglesTwoQuads,
}

/// One tensor-product Bernstein-Bezier patch.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolyhedralPatch {
    pub family: PolyhedralSplineFamily,
    pub degree_u: usize,
    pub degree_v: usize,
    /// Row-major grid with `u` as the outer index and `v` as the inner index.
    pub control_points: Vec<Point3>,
}

impl PolyhedralPatch {
    pub fn new(
        family: PolyhedralSplineFamily,
        degree_u: usize,
        degree_v: usize,
        control_points: Vec<Point3>,
    ) -> Result<Self, PolyhedralSplineError> {
        let expected = (degree_u + 1) * (degree_v + 1);
        if control_points.len() != expected {
            return Err(PolyhedralSplineError::PatchControlPointCount {
                family,
                degree_u,
                degree_v,
                expected,
                got: control_points.len(),
            });
        }
        Ok(Self {
            family,
            degree_u,
            degree_v,
            control_points,
        })
    }

    pub fn cp(&self, u: usize, v: usize) -> Point3 {
        self.control_points[u * (self.degree_v + 1) + v]
    }

    pub fn set_cp(&mut self, u: usize, v: usize, p: Point3) {
        self.control_points[u * (self.degree_v + 1) + v] = p;
    }

    pub fn evaluate(&self, u: f64, v: f64) -> Point3 {
        let bu = bernstein_basis(self.degree_u, u);
        let bv = bernstein_basis(self.degree_v, v);
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for (i, &bi) in bu.iter().enumerate() {
            for (j, &bj) in bv.iter().enumerate() {
                let p = self.cp(i, j);
                out = out + (p - Point3::new(0.0, 0.0, 0.0)) * (bi * bj);
            }
        }
        Point3::new(out.x, out.y, out.z)
    }

    pub fn partial_u(&self, u: f64, v: f64) -> Vec3 {
        if self.degree_u == 0 {
            return Vec3::new(0.0, 0.0, 0.0);
        }
        let bu = bernstein_basis(self.degree_u - 1, u);
        let bv = bernstein_basis(self.degree_v, v);
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for (i, &bi) in bu.iter().enumerate() {
            for (j, &bj) in bv.iter().enumerate() {
                let delta = self.cp(i + 1, j) - self.cp(i, j);
                out = out + delta * (self.degree_u as f64 * bi * bj);
            }
        }
        out
    }

    pub fn partial_v(&self, u: f64, v: f64) -> Vec3 {
        if self.degree_v == 0 {
            return Vec3::new(0.0, 0.0, 0.0);
        }
        let bu = bernstein_basis(self.degree_u, u);
        let bv = bernstein_basis(self.degree_v - 1, v);
        let mut out = Vec3::new(0.0, 0.0, 0.0);
        for (i, &bi) in bu.iter().enumerate() {
            for (j, &bj) in bv.iter().enumerate() {
                let delta = self.cp(i, j + 1) - self.cp(i, j);
                out = out + delta * (self.degree_v as f64 * bi * bj);
            }
        }
        out
    }

    pub fn degree_raised_to_bicubic(&self) -> Self {
        let mut patch = self.clone();
        while patch.degree_u < 3 {
            patch = patch.raise_u_once();
        }
        while patch.degree_v < 3 {
            patch = patch.raise_v_once();
        }
        patch
    }

    pub fn to_bezier_nurbs(&self) -> NurbsSurface {
        NurbsSurface {
            degree_u: self.degree_u,
            degree_v: self.degree_v,
            knots_u: bezier_knots(self.degree_u),
            knots_v: bezier_knots(self.degree_v),
            control_points: self.control_points.clone(),
            num_u: self.degree_u + 1,
            num_v: self.degree_v + 1,
            weights: vec![1.0; self.control_points.len()],
        }
    }

    fn raise_u_once(&self) -> Self {
        let old_degree = self.degree_u;
        let new_degree = old_degree + 1;
        let mut out = Vec::with_capacity((new_degree + 1) * (self.degree_v + 1));
        for i in 0..=new_degree {
            for j in 0..=self.degree_v {
                out.push(raise_point_1d(
                    i,
                    new_degree,
                    old_degree,
                    self.cp(i.saturating_sub(1), j),
                    self.cp(i.min(old_degree), j),
                ));
            }
        }
        Self {
            family: self.family,
            degree_u: new_degree,
            degree_v: self.degree_v,
            control_points: out,
        }
    }

    fn raise_v_once(&self) -> Self {
        let old_degree = self.degree_v;
        let new_degree = old_degree + 1;
        let mut out = Vec::with_capacity((self.degree_u + 1) * (new_degree + 1));
        for i in 0..=self.degree_u {
            for j in 0..=new_degree {
                out.push(raise_point_1d(
                    j,
                    new_degree,
                    old_degree,
                    self.cp(i, j.saturating_sub(1)),
                    self.cp(i, j.min(old_degree)),
                ));
            }
        }
        Self {
            family: self.family,
            degree_u: self.degree_u,
            degree_v: new_degree,
            control_points: out,
        }
    }
}

/// A constructed local polyhedral spline complex.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PolyhedralSpline {
    pub family: PolyhedralSplineFamily,
    pub patches: Vec<PolyhedralPatch>,
}

/// Tensor-product patch edge selector.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PatchEdge {
    UMin,
    UMax,
    VMin,
    VMax,
}

/// Construction and validation failures.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
pub enum PolyhedralSplineError {
    #[error("{family:?} stencil has {got} points; expected {expected}")]
    StencilPointCount {
        family: PolyhedralSplineFamily,
        expected: usize,
        got: usize,
    },
    #[error(
        "{family:?} patch degree ({degree_u},{degree_v}) expects {expected} control points; got {got}"
    )]
    PatchControlPointCount {
        family: PolyhedralSplineFamily,
        degree_u: usize,
        degree_v: usize,
        expected: usize,
        got: usize,
    },
    #[error("unsupported {family} valence/sides/sectors {value}; supported: {supported}")]
    UnsupportedValence {
        family: &'static str,
        value: usize,
        supported: &'static str,
    },
    #[error(
        "mask table {name} has {got_rows}x{got_cols}; expected {expected_rows}x{expected_cols}"
    )]
    MaskShape {
        name: &'static str,
        expected_rows: usize,
        expected_cols: usize,
        got_rows: usize,
        got_cols: usize,
    },
    #[error("mask table {name} contains an invalid number at row {row}, column {col}")]
    MaskParse {
        name: &'static str,
        row: usize,
        col: usize,
    },
    #[error(
        "PnS {family} mask tables unavailable: the ACM ToMS 1032 package \
         (docs/code/ACM-ToMS-PnS-1032, deliberately untracked — license) \
         was absent when this build compiled, so the tables were not \
         embedded"
    )]
    ToMsTablesUnavailable { family: &'static str },
    #[error("{family:?} contains a non-finite control point at index {index}")]
    NonFinitePoint {
        family: PolyhedralSplineFamily,
        index: usize,
    },
    #[error("control net face {face} references missing vertex {vertex}")]
    MissingVertex { face: usize, vertex: usize },
    #[error("control net contains duplicate directed edge {from}->{to}")]
    DuplicateDirectedEdge { from: usize, to: usize },
    #[error(
        "control net vertex {vertex} is a boundary pinch: {boundary_edges} boundary edges meet there (a manifold rim has exactly 2)"
    )]
    BoundaryPinch {
        vertex: usize,
        boundary_edges: usize,
    },
    #[error("control net layout cannot provide stencil for {family:?}: {reason}")]
    StencilTopology {
        family: PolyhedralSplineFamily,
        reason: &'static str,
    },
    #[error("OBJ parse error at line {line}: {reason}")]
    ObjParse { line: usize, reason: String },
    #[error("{family:?} seam validation failed: {reason}")]
    SeamValidation {
        family: PolyhedralSplineFamily,
        reason: String,
    },
    #[error(
        "{family:?} patch {patch_index} interior Jacobian is degenerate: min normal length {min_normal_length}"
    )]
    DegeneratePatchInterior {
        family: PolyhedralSplineFamily,
        patch_index: usize,
        min_normal_length: f64,
    },
    #[error("{family:?} patch atlas is not BRep-ready: {reason}")]
    AtlasTopology {
        family: PolyhedralSplineFamily,
        reason: String,
    },
    #[error("{family:?} patch {patch_index} has invalid parameter domain: {reason}")]
    InvalidPatchDomain {
        family: PolyhedralSplineFamily,
        patch_index: usize,
        reason: String,
    },
    #[error("{family:?} BRep export failed: {reason}")]
    BrepExport {
        family: PolyhedralSplineFamily,
        reason: String,
    },
}

fn raise_point_1d(
    i: usize,
    new_degree: usize,
    _old_degree: usize,
    prev: Point3,
    cur: Point3,
) -> Point3 {
    if i == 0 {
        return cur;
    }
    if i == new_degree {
        return prev;
    }
    let a = i as f64 / new_degree as f64;
    let b = 1.0 - a;
    Point3::new(
        a * prev.x + b * cur.x,
        a * prev.y + b * cur.y,
        a * prev.z + b * cur.z,
    )
}

fn bernstein_basis(degree: usize, t: f64) -> Vec<f64> {
    let mut out = Vec::with_capacity(degree + 1);
    for i in 0..=degree {
        out.push(binom(degree, i) * t.powi(i as i32) * (1.0 - t).powi((degree - i) as i32));
    }
    out
}

fn binom(n: usize, k: usize) -> f64 {
    if k > n {
        return 0.0;
    }
    let mut num = 1.0;
    let mut den = 1.0;
    for i in 0..k {
        num *= (n - i) as f64;
        den *= (i + 1) as f64;
    }
    num / den
}

fn bezier_knots(degree: usize) -> Vec<f64> {
    let mut knots = vec![0.0; degree + 1];
    knots.extend(std::iter::repeat_n(1.0, degree + 1));
    knots
}
