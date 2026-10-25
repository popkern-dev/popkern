//! A generative design as a VTK XML unstructured grid (`.vtu`, ASCII):
//! one hexahedron per design entry (its grid cell; the fragments of a
//! split cell overlap), with cell data `density` (the physical density
//! x̃), `fill` (the entry's volume over its cell's) and `floor` (the
//! kept material). In ParaView, a threshold on `density` shows the
//! design.

use std::fmt::Write;

use popkern_analysis::structural::generative::Design;

/// The `.vtu` text of `design`.
pub fn write_design_vtu(design: &Design) -> String {
    let h = design.cell_size;
    let n = design.entries.len();
    let mut out = String::new();
    out.push_str("<?xml version=\"1.0\"?>\n");
    out.push_str(
        "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\">\n",
    );
    out.push_str("  <UnstructuredGrid>\n");
    let _ = writeln!(
        out,
        "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{n}\">",
        8 * n
    );
    out.push_str("      <Points>\n");
    out.push_str(
        "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
    );
    // VTK_HEXAHEDRON corner order: the bottom face counter-clockwise,
    // then the top.
    const CORNERS: [[f64; 3]; 8] = [
        [0.0, 0.0, 0.0],
        [1.0, 0.0, 0.0],
        [1.0, 1.0, 0.0],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
        [1.0, 0.0, 1.0],
        [1.0, 1.0, 1.0],
        [0.0, 1.0, 1.0],
    ];
    for e in &design.entries {
        for c in CORNERS {
            let p = [0, 1, 2].map(|a| e.corner[a] + h * c[a]);
            let _ = writeln!(out, "          {:e} {:e} {:e}", p[0], p[1], p[2]);
        }
    }
    out.push_str("        </DataArray>\n      </Points>\n      <Cells>\n");
    out.push_str("        <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
    for i in 0..n {
        let b = 8 * i;
        let _ = writeln!(
            out,
            "          {} {} {} {} {} {} {} {}",
            b,
            b + 1,
            b + 2,
            b + 3,
            b + 4,
            b + 5,
            b + 6,
            b + 7
        );
    }
    out.push_str("        </DataArray>\n");
    out.push_str(
        "        <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n          ",
    );
    for i in 0..n {
        let _ = write!(out, "{}{}", if i > 0 { " " } else { "" }, 8 * (i + 1));
    }
    out.push_str("\n        </DataArray>\n");
    out.push_str("        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n          ");
    for i in 0..n {
        out.push_str(if i > 0 { " 12" } else { "12" });
    }
    out.push_str("\n        </DataArray>\n      </Cells>\n      <CellData Scalars=\"density\">\n");
    let cell = h * h * h;
    let fill: Vec<f64> = design
        .entries
        .iter()
        .map(|e| if cell > 0.0 { e.volume / cell } else { 0.0 })
        .collect();
    for (name, values) in [
        ("density", design.densities.as_slice()),
        ("fill", fill.as_slice()),
        ("floor", design.floors.as_slice()),
    ] {
        let _ = writeln!(
            out,
            "        <DataArray type=\"Float64\" Name=\"{name}\" format=\"ascii\">"
        );
        out.push_str("          ");
        for (i, v) in values.iter().enumerate() {
            let _ = write!(out, "{}{v:e}", if i > 0 { " " } else { "" });
        }
        out.push_str("\n        </DataArray>\n");
    }
    out.push_str("      </CellData>\n    </Piece>\n  </UnstructuredGrid>\n</VTKFile>\n");
    out
}
