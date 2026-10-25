//! VTK XML unstructured grid (`.vtu`, ASCII) of the analysis mesh with
//! per-vertex fields, for ParaView until popkern live has the overlay.
//!
//! Every solid's display mesh is written as triangles, solids
//! concatenated in order. Point data: `displacement` (3 components, mm),
//! `von_mises` and `max_principal` (MPa), `zz_error` (relative; −1 where
//! no estimate is honest), `temperature` (°C), `heat_flux` (3
//! components, W/m²). Cell data: `solid` and
//! `face` ids. Numbers are written in shortest round-trip form.

use std::fmt::Write;

use popkern_analysis::model::AnalysisModel;
use popkern_analysis::recover::SolidStresses;

/// `zz_error` value where a vertex has no estimate.
pub const NO_ESTIMATE: f64 = -1.0;

/// The fields to write; each is per solid, per vertex.
pub struct VtuFields<'a> {
    /// Displacement, mm.
    pub displacement: Option<&'a [Vec<[f64; 3]>]>,
    /// Recovered stresses.
    pub stresses: Option<&'a [SolidStresses]>,
    /// Temperature, °C.
    pub temperature: Option<&'a [Vec<f64>]>,
    /// Heat flux, W/m².
    pub heat_flux: Option<&'a [Vec<[f64; 3]>]>,
}

fn array(out: &mut String, name: &str, kind: &str, components: usize, values: &[f64]) {
    let _ = writeln!(
        out,
        "        <DataArray type=\"{kind}\" Name=\"{name}\" NumberOfComponents=\"{components}\" format=\"ascii\">"
    );
    out.push_str("          ");
    for (i, v) in values.iter().enumerate() {
        if i > 0 {
            out.push(' ');
        }
        if kind == "Float64" {
            let _ = write!(out, "{v:e}");
        } else {
            let _ = write!(out, "{}", *v as i64);
        }
    }
    out.push_str("\n        </DataArray>\n");
}

/// The `.vtu` text of `model` with `fields`.
pub fn write_vtu(model: &AnalysisModel, fields: &VtuFields<'_>) -> String {
    let points: usize = model.solids.iter().map(|s| s.vertices.len()).sum();
    let cells: usize = model.solids.iter().map(|s| s.triangles.len()).sum();
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\"?>\n");
    out.push_str(
        "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n",
    );
    out.push_str("  <UnstructuredGrid>\n");
    let _ = writeln!(
        out,
        "    <Piece NumberOfPoints=\"{points}\" NumberOfCells=\"{cells}\">"
    );
    out.push_str("      <Points>\n");
    let xyz: Vec<f64> = model
        .solids
        .iter()
        .flat_map(|s| s.vertices.iter().flatten().copied())
        .collect();
    array(&mut out, "points", "Float64", 3, &xyz);
    out.push_str("      </Points>\n      <Cells>\n");
    let mut connectivity = Vec::with_capacity(3 * cells);
    let mut base = 0usize;
    for s in &model.solids {
        for t in &s.triangles {
            connectivity.extend(t.iter().map(|&v| (base + v as usize) as f64));
        }
        base += s.vertices.len();
    }
    array(&mut out, "connectivity", "Int64", 1, &connectivity);
    let offsets: Vec<f64> = (1..=cells).map(|c| (3 * c) as f64).collect();
    array(&mut out, "offsets", "Int64", 1, &offsets);
    array(&mut out, "types", "UInt8", 1, &vec![5.0; cells]);
    out.push_str("      </Cells>\n      <PointData>\n");
    if let Some(d) = fields.displacement {
        let v: Vec<f64> = d.iter().flatten().flatten().copied().collect();
        array(&mut out, "displacement", "Float64", 3, &v);
    }
    if let Some(st) = fields.stresses {
        let all = || st.iter().flat_map(|s| s.vertices.iter());
        let vm: Vec<f64> = all().map(|v| v.von_mises).collect();
        let s1: Vec<f64> = all().map(|v| v.max_principal).collect();
        let zz: Vec<f64> = all().map(|v| v.error.unwrap_or(NO_ESTIMATE)).collect();
        array(&mut out, "von_mises", "Float64", 1, &vm);
        array(&mut out, "max_principal", "Float64", 1, &s1);
        array(&mut out, "zz_error", "Float64", 1, &zz);
    }
    if let Some(t) = fields.temperature {
        let v: Vec<f64> = t.iter().flatten().copied().collect();
        array(&mut out, "temperature", "Float64", 1, &v);
    }
    if let Some(q) = fields.heat_flux {
        let v: Vec<f64> = q.iter().flatten().flatten().copied().collect();
        array(&mut out, "heat_flux", "Float64", 3, &v);
    }
    out.push_str("      </PointData>\n      <CellData>\n");
    let solid: Vec<f64> = model
        .solids
        .iter()
        .flat_map(|s| std::iter::repeat_n(f64::from(s.id.0), s.triangles.len()))
        .collect();
    let face: Vec<f64> = model
        .solids
        .iter()
        .flat_map(|s| s.triangle_face.iter().map(|f| f64::from(f.0)))
        .collect();
    array(&mut out, "solid", "Int64", 1, &solid);
    array(&mut out, "face", "Int64", 1, &face);
    out.push_str("      </CellData>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n");
    out
}
