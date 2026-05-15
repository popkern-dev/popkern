use super::types::{PolyhedralPatch, PolyhedralSplineError, PolyhedralSplineFamily};
use crate::geometry::Point3;

#[derive(Debug, Clone, Copy)]
pub struct MaskSpec {
    pub name: &'static str,
    pub rows: usize,
    pub cols: usize,
    pub text: &'static str,
}

pub fn apply_mask_to_patches(
    spec: MaskSpec,
    family: PolyhedralSplineFamily,
    degree_u: usize,
    degree_v: usize,
    points: &[Point3],
) -> Result<Vec<PolyhedralPatch>, PolyhedralSplineError> {
    if points.len() != spec.cols {
        return Err(PolyhedralSplineError::StencilPointCount {
            family,
            expected: spec.cols,
            got: points.len(),
        });
    }
    let values = parse_mask(spec)?;
    let out = multiply(&values, spec.rows, spec.cols, points);
    rows_to_patches(family, degree_u, degree_v, out)
}

pub fn apply_inline_mask_to_patch(
    mask: &[[f64; 9]; 9],
    family: PolyhedralSplineFamily,
    points: &[Point3; 9],
) -> Result<PolyhedralPatch, PolyhedralSplineError> {
    let mut out = Vec::with_capacity(9);
    for row in mask {
        let mut p = Point3::new(0.0, 0.0, 0.0);
        for (w, q) in row.iter().zip(points.iter()) {
            p.x += w * q.x;
            p.y += w * q.y;
            p.z += w * q.z;
        }
        out.push(p);
    }
    PolyhedralPatch::new(family, 2, 2, out)
}

fn parse_mask(spec: MaskSpec) -> Result<Vec<f64>, PolyhedralSplineError> {
    let mut values = Vec::with_capacity(spec.rows * spec.cols);
    let mut got_rows = 0;
    let mut got_cols = None;
    for (row_idx, line) in spec
        .text
        .lines()
        .filter(|line| {
            let trimmed = line.trim();
            !trimmed.is_empty() && !trimmed.starts_with('#')
        })
        .enumerate()
    {
        got_rows += 1;
        let cols: Vec<&str> = line.split(',').map(str::trim).collect();
        got_cols.get_or_insert(cols.len());
        if cols.len() != spec.cols {
            return Err(PolyhedralSplineError::MaskShape {
                name: spec.name,
                expected_rows: spec.rows,
                expected_cols: spec.cols,
                got_rows,
                got_cols: cols.len(),
            });
        }
        for (col_idx, col) in cols.iter().enumerate() {
            values.push(
                col.parse::<f64>()
                    .map_err(|_| PolyhedralSplineError::MaskParse {
                        name: spec.name,
                        row: row_idx,
                        col: col_idx,
                    })?,
            );
        }
    }
    let got_cols = got_cols.unwrap_or(0);
    if got_rows != spec.rows || got_cols != spec.cols {
        return Err(PolyhedralSplineError::MaskShape {
            name: spec.name,
            expected_rows: spec.rows,
            expected_cols: spec.cols,
            got_rows,
            got_cols,
        });
    }
    Ok(values)
}

fn multiply(values: &[f64], rows: usize, cols: usize, points: &[Point3]) -> Vec<Point3> {
    let mut out = Vec::with_capacity(rows);
    for r in 0..rows {
        let mut p = Point3::new(0.0, 0.0, 0.0);
        for c in 0..cols {
            let w = values[r * cols + c];
            let q = points[c];
            p.x += w * q.x;
            p.y += w * q.y;
            p.z += w * q.z;
        }
        out.push(p);
    }
    out
}

fn rows_to_patches(
    family: PolyhedralSplineFamily,
    degree_u: usize,
    degree_v: usize,
    rows: Vec<Point3>,
) -> Result<Vec<PolyhedralPatch>, PolyhedralSplineError> {
    let per_patch = (degree_u + 1) * (degree_v + 1);
    let mut patches = Vec::with_capacity(rows.len() / per_patch);
    for chunk in rows.chunks(per_patch) {
        patches.push(PolyhedralPatch::new(
            family,
            degree_u,
            degree_v,
            chunk.to_vec(),
        )?);
    }
    Ok(patches)
}
