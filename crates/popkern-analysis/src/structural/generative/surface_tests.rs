//! The design surface (PH-117): closed, oriented, and at the design's
//! outline.

use super::*;
use crate::structural::generative::{Design, DesignEntry};

/// A design of `n` cells of size `h` filling a box, at density `rho`.
fn box_design(n: [usize; 3], h: f64, rho: impl Fn([usize; 3]) -> f64) -> Design {
    let mut entries = Vec::new();
    let mut densities = Vec::new();
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                let corner = [i, j, k].map(|q| q as f64 * h);
                entries.push(DesignEntry {
                    solid: 0,
                    entry: entries.len(),
                    centre: corner.map(|c| c + 0.5 * h),
                    volume: h * h * h,
                    corner,
                });
                densities.push(rho([i, j, k]));
            }
        }
    }
    Design {
        floors: vec![0.0; entries.len()],
        entries,
        densities,
        cell_size: h,
        compliance: 0.0,
        volume_fraction: 1.0,
        history: Vec::new(),
        converged: true,
        unsupported_fraction: None,
        grey_fraction: 0.0,
        eroded_compliance: None,
        member_thickness: 0.0,
        report: crate::report::SolveReport::new(crate::options::AnalysisOptions::default()),
    }
}

/// A solid 8 × 4 × 2 box: the surface is closed and oriented, its faces
/// on the box's, and its volume converges to the box's as the cells
/// shrink, at second order (its edges are rounded by half a cell).
/// Measured: −5.1%, −1.3%, −0.34% at 1, 0.5 and 0.25 mm cells (−30%,
/// −8.3%, −2.2% when the field was averaged to the cell corners).
#[test]
fn a_solid_box_closes_and_converges_to_its_volume() {
    let mut errors = Vec::new();
    for h in [1.0, 0.5, 0.25] {
        let n = [8.0, 4.0, 2.0].map(|l: f64| (l / h).round() as usize);
        let m = design_surface(&box_design(n, h, |_| 1.0), 0.5);
        assert!(m.is_closed(), "h {h}");
        errors.push(m.volume() / 64.0 - 1.0);
    }
    assert!(errors.iter().all(|e| *e < 0.0), "{errors:?}");
    assert!(
        errors[1] / errors[0] < 0.33 && errors[2] / errors[1] < 0.33,
        "{errors:?}"
    );
    assert!(
        errors[0].abs() < 0.06 && errors[2].abs() < 0.004,
        "{errors:?}"
    );
}

/// A design with a void through it (density 0 in a slot) keeps the slot:
/// the surface is closed and its volume drops by the slot's. The
/// negative control: one triangle removed, the check sees the hole.
#[test]
fn a_slotted_design_keeps_its_slot_and_a_hole_is_seen() {
    let h = 0.25;
    let n = [32, 16, 8];
    // A slot 2 mm wide through the box at x in [3, 5).
    let slot = |q: [usize; 3]| {
        let x = q[0] as f64 * h;
        if (3.0..5.0).contains(&x) { 0.0 } else { 1.0 }
    };
    let solid = design_surface(&box_design(n, h, |_| 1.0), 0.5);
    let slotted = design_surface(&box_design(n, h, slot), 0.5);
    assert!(slotted.is_closed());
    // Two pieces of 3 × 4 × 2 and 3 × 4 × 2 against 8 × 4 × 2.
    let ratio = slotted.volume() / solid.volume();
    assert!((ratio - 0.75).abs() < 0.03, "{ratio}");
    let mut open = slotted.clone();
    open.triangles.pop();
    assert!(!open.is_closed());
}

/// The design as a part (`body`): on a solid box design each of the
/// box's six faces keeps its triangles (its area, less the rounded
/// edges), and only the rounding is new surface; a slot's two walls are
/// new surface. Measured at 0.25 mm cells: the faces keep 91–95% of
/// their area, the rounding is 4.8 mm² of 112, the slot walls 15.8 of
/// 16 mm² (their own edges rounded); averaged to the cell corners,
/// 78–89%, 9.2 mm² and 14.6 mm².
#[test]
fn the_design_body_inherits_the_part_faces() {
    use crate::mesh_test_helpers::{box_mesh, model_of};
    use crate::structural::generative::body::design_model;
    let model = model_of(box_mesh([0.0; 3], [8.0, 4.0, 2.0]), 6);
    let h = 0.25;
    let n = [32, 16, 8];
    let true_area = [8.0, 8.0, 16.0, 16.0, 32.0, 32.0];
    let perimeter = [12.0, 12.0, 20.0, 20.0, 24.0, 24.0];
    let solid = design_model(&model, &design_surface(&box_design(n, h, |_| 1.0), 0.5), h);
    assert_eq!(solid.faces.len(), 7);
    solid.validate().unwrap_or_else(|e| panic!("{e}"));
    for (f, area) in true_area.iter().enumerate() {
        let got = solid.faces[f].area;
        // The rounded band along the face's edges is at most a cell wide.
        let floor = area - perimeter[f] * h;
        assert!(
            got >= floor && got <= area * (1.0 + 1e-9),
            "face {f}: {got} of {area}"
        );
    }
    let rounding = solid.faces[6].area;
    assert!(rounding < 0.1 * 112.0, "{rounding}");
    let slot = |q: [usize; 3]| {
        let x = q[0] as f64 * h;
        if (3.0..5.0).contains(&x) { 0.0 } else { 1.0 }
    };
    let slotted = design_model(&model, &design_surface(&box_design(n, h, slot), 0.5), h);
    let walls = slotted.faces[6].area - rounding;
    assert!((walls - 16.0).abs() < 0.15 * 16.0, "{walls}");
}

/// A 2 × 2-cell column through a 6 mm cube, standing on its base and
/// reaching its top. Every layer of cell centres, the first and last
/// (against the cube's faces) included, carries the column's full
/// section, x and y from 2 to 4: the column keeps its width where it
/// meets the part's boundary. (Averaged to corner nodes, the section
/// was a diamond of half the area, shrinking to a point at the faces.)
/// Its volume is the column's less the half-cell rounding of its
/// edges.
#[test]
fn a_column_keeps_its_section_where_it_meets_the_part_boundary() {
    let column = |q: [usize; 3]| {
        if (2..4).contains(&q[0]) && (2..4).contains(&q[1]) {
            1.0
        } else {
            0.0
        }
    };
    let m = design_surface(&box_design([6, 6, 6], 1.0, column), 0.5);
    assert!(m.is_closed());
    for layer in 0..6 {
        let z = layer as f64 + 0.5;
        let on: Vec<[f64; 3]> = m
            .vertices
            .iter()
            .filter(|v| (v[2] - z).abs() < 1e-9)
            .copied()
            .collect();
        for a in 0..2 {
            let lo = on.iter().map(|v| v[a]).fold(f64::INFINITY, f64::min);
            let hi = on.iter().map(|v| v[a]).fold(f64::NEG_INFINITY, f64::max);
            assert!(
                (lo - 2.0).abs() < 1e-9 && (hi - 4.0).abs() < 1e-9,
                "z {z} axis {a}: {lo}..{hi}"
            );
        }
    }
    let v = m.volume();
    assert!(v > 0.8 * 24.0 && v <= 24.0, "{v}");
}

/// The material grid marches to the design's own surface (a host
/// marching a recorded design sees what the solver's surface shows),
/// and a grid whose material does not fill its dimensions marches to
/// nothing rather than reading out of bounds.
#[test]
fn a_material_grid_marches_to_the_design_surface() {
    let design = box_design(
        [6, 3, 3],
        0.5,
        |[i, j, _]| if i + j < 5 { 0.9 } else { 0.2 },
    );
    let grid = MaterialGrid::new(&design.entries, &design.densities, design.cell_size)
        .unwrap_or_else(|| panic!("no grid"));
    assert_eq!(grid.dims, [8, 5, 5]);
    assert_eq!(grid_surface(&grid, 0.5), design_surface(&design, 0.5));
    let mut short = grid.clone();
    short.material.pop();
    assert!(grid_surface(&short, 0.5).triangles.is_empty());
    let flat = MaterialGrid {
        dims: [1, 5, 5],
        material: vec![1.0; 25],
        ..grid
    };
    assert!(grid_surface(&flat, 0.5).triangles.is_empty());
}
