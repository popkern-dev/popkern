//! CC1 (PH-031): a slot narrower than a cell does not weld the material
//! on its two sides.

use super::assemble::FaceMechanics;
use super::*;
use crate::element::vector_boundary::PieceSupport;
use crate::mesh_test_helpers::{extrude_faces, model_of};
use crate::moments::Splitting;

/// A 6 × 2 × `t` plate clamped at x = 0 with a slot of width `w`
/// centred at y = 1 from x = 1 to the free end: two prongs on a root.
/// Faces: 0/1 z sides, 2 other sides, 3 lower prong end, 4 upper prong
/// end, 5 the clamped end.
pub(super) fn slotted(w: f64, t: f64) -> crate::model::AnalysisModel {
    let (a, b) = (1.0 - 0.5 * w, 1.0 + 0.5 * w);
    let outline = [
        [0.0, 0.0],
        [6.0, 0.0],
        [6.0, a],
        [1.0, a],
        [1.0, b],
        [6.0, b],
        [6.0, 2.0],
        [0.0, 2.0],
    ];
    let tris = [
        [0, 1, 2],
        [0, 2, 3],
        [0, 3, 4],
        [0, 4, 7],
        [4, 5, 6],
        [4, 6, 7],
    ];
    let mut mesh = extrude_faces(&outline, &tris, 0.0, t, &[2, 3, 2, 2, 2, 4, 2, 5]);
    mesh.exact_volume = Some(t * (12.0 - 5.0 * w));
    model_of(mesh, 6)
}

/// Plate thickness: 4 cells at h = 0.25, so the base discretisation
/// error stays small and CC1 measures the slit handling.
const T: f64 = 1.0;

/// The slot's opening (upper tip minus lower tip, along y) when the
/// prong ends are pulled apart by ±1e-3, at cell size `h`.
fn opening(w: f64, h: f64, shift: [f64; 3], split: Splitting) -> (f64, StructuralResult) {
    let model = slotted(w, T);
    let mut faces = vec![FaceMechanics::Free; 6];
    faces[5] = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    });
    faces[3] = FaceMechanics::Traction([0.0, -1e-3, 0.0]);
    faces[4] = FaceMechanics::Traction([0.0, 1e-3, 0.0]);
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces,
    };
    let options = AnalysisOptions {
        target_cell_size: Some(h),
        grid_shift: shift,
        max_active_cells: 200_000,
        coarsest_dofs: 600,
        cg_max_iters: 400,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic_with(
        &model,
        &problem,
        &options,
        crate::recover::Fallback::InCell,
        None,
        split,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let (a, b) = (1.0 - 0.5 * w, 1.0 + 0.5 * w);
    let at = |y: f64| {
        r.field
            .displacement_at([5.9, y, 0.5 * T])
            .unwrap_or_else(|e| panic!("{e}"))
            .unwrap_or_else(|| panic!("tip outside"))[1]
    };
    let open = at(0.5 * (b + 2.0)) - at(0.5 * a);
    (open, r)
}

const SHIFTS: [[f64; 3]; 5] = [
    [0.0; 3],
    [0.0, 0.125, 0.0],
    [0.07, 0.11, 0.03],
    [0.0, 0.05, 0.0],
    [0.13, 0.2, 0.09],
];

#[test]
fn cc1_a_slot_of_0_3h_matches_the_resolved_slot_and_welding_is_the_control() {
    // Each opening is taken against the converged solution of its own
    // plate (gap > one cell, so no cell spans it). The two-bar
    // reference on the same grid is the plate with a resolved 2-cell
    // slot; the 0.3h slot must read like it (C2 tolerance, 5%).
    // Measured at h = 0.25 over the shifts: resolved 0.981, 0.974,
    // 0.935, 0.972, 0.920; slit (exact fragment moments) 0.984, 0.986,
    // 0.907, 0.984, 0.910; welded 0.0018.
    let h = 0.25;
    let w = 0.3 * h;
    let reference = opening(w, w / 1.1, [0.0; 3], Splitting::Exact).0;
    let resolved_reference = opening(0.5, 0.0625, [0.0; 3], Splitting::Exact).0;
    for shift in SHIFTS {
        let (slit, _) = opening(w, h, shift, Splitting::Exact);
        let resolved = opening(0.5, h, shift, Splitting::Exact).0 / resolved_reference;
        let slit = slit / reference;
        assert!(
            (slit - resolved).abs() <= 0.05,
            "{shift:?}: {slit} vs {resolved}"
        );
        // Negative control: without splitting the prongs weld.
        let welded = opening(w, h, shift, Splitting::Welded).0 / reference;
        assert!(welded.abs() < 0.01, "{shift:?}: welded {welded}");
    }
}

#[test]
fn fragment_moments_sum_to_the_cell_and_stay_near_their_lattice_volume() {
    use crate::clip::clip_solid;
    use crate::moments::fragments::SPLIT_DEPTH;
    use crate::moments::rasterise_with;
    let model = slotted(0.075, T);
    let mesh = &model.solids[0];
    let h = 0.25;
    let mut split_cells = 0;
    for shift in SHIFTS {
        let options = AnalysisOptions {
            target_cell_size: Some(h),
            grid_shift: shift,
            ..AnalysisOptions::default()
        };
        let (mut lo, mut hi) = ([f64::INFINITY; 3], [f64::NEG_INFINITY; 3]);
        for v in &mesh.vertices {
            for a in 0..3 {
                lo[a] = lo[a].min(v[a]);
                hi[a] = hi[a].max(v[a]);
            }
        }
        let grid = Grid::for_model(lo, hi, &options).unwrap_or_else(|e| panic!("{e}"));
        let clipped = clip_solid(&grid, mesh).unwrap_or_else(|e| panic!("{e}"));
        let split = rasterise_with(
            &grid,
            mesh.id,
            &clipped,
            2,
            &|_| false,
            crate::moments::Splitting::Exact,
        )
        .unwrap_or_else(|e| panic!("{e}"));
        let whole = rasterise_with(
            &grid,
            mesh.id,
            &clipped,
            2,
            &|_| false,
            crate::moments::Splitting::Welded,
        )
        .unwrap_or_else(|e| panic!("{e}"));
        let s = h / f64::from(1u32 << SPLIT_DEPTH);
        for sc in &split.fragments.split {
            split_cells += 1;
            let entries = split.entries(sc.cell);
            assert!(entries.len() >= 2);
            let exact = &whole.cells[whole.entries(sc.cell).start].m;
            for (q, x) in exact.iter().enumerate() {
                let sum: f64 = entries.clone().map(|e| split.cells[e].m[q]).sum();
                assert!((sum - x).abs() <= 1e-12 * h * h * h, "{q}: {sum} vs {x}");
            }
            for e in entries {
                let tag = split.cells[e].fragment + 1;
                let lattice = sc.labels.iter().filter(|&&l| l == tag).count() as f64 * s * s * s;
                // Within one lattice layer over the cell's cross-section.
                assert!((split.cells[e].m[0] - lattice).abs() <= s * h * h, "{e}");
            }
        }
    }
    assert!(split_cells > 0, "no shift split a cell");
}

#[test]
fn the_resolved_slot_converges_at_its_re_entrant_root() {
    // PH-043: at h = 0.25 (3 cells across a prong) the resolved slot
    // reads 0.920-0.981 of converged. Isolated by moving the root wall
    // to x = 1.07 and probing shifts one axis at a time: tip sampling
    // (point vs face-averaged opening agree to 0.1%) and the clamp (cut
    // clamp, aligned root: 0.987) are not the cause. The loss comes from
    // the root's re-entrant corner cell cut in both x and y, which
    // stiffens the prong root's rotation (root cut + y cut 0.940 vs
    // aligned root 0.974), plus up to 3 points from aggregation at
    // gamma 0.25 at one shift (0.952 at gamma 0.1). It is ordinary
    // discretisation error: at h = 0.125 every shift is within 2.2%.
    // PH-044: roots chosen by material centroid keep the root corner's
    // prong sliver on its prong: 0.981, 0.974, 0.935, 0.972, 0.945
    // (the last was 0.920; 0.951 is the gamma -> 0 limit).
    let reference = opening(0.5, 0.0625, [0.0; 3], Splitting::Exact).0;
    for shift in SHIFTS {
        let coarse = opening(0.5, 0.25, shift, Splitting::Exact).0 / reference;
        let fine = opening(0.5, 0.125, shift, Splitting::Exact).0 / reference;
        assert!(coarse > 0.93, "{shift:?}: {coarse}");
        assert!((fine - 1.0).abs() <= 0.03, "{shift:?}: {fine}");
        assert!((fine - 1.0).abs() < (coarse - 1.0).abs(), "{shift:?}");
    }
}

/// Moments of the box `[lo, hi]` clipped to the cell at `cell_lo`, in
/// cell-local coordinates: the closed form of a split cell's fragment.
pub(super) fn box_moments(lo: [f64; 3], hi: [f64; 3], cell_lo: [f64; 3], h: f64) -> [f64; 27] {
    let one_d = |a: usize| -> [f64; 3] {
        let l = |x: f64| (2.0 * (x - cell_lo[a]) / h - 1.0).clamp(-1.0, 1.0);
        let (p, q) = (l(lo[a]), l(hi[a]));
        // ∫ ξ^k dx with dx = (h/2) dξ.
        [1i32, 2, 3].map(|k| 0.5 * h * (q.powi(k) - p.powi(k)) / f64::from(k))
    };
    let (x, y, z) = (one_d(0), one_d(1), one_d(2));
    let mut m = [0.0; 27];
    for c in 0..3 {
        for b in 0..3 {
            for a in 0..3 {
                m[crate::moments::moment_index(a, b, c)] = x[a] * y[b] * z[c];
            }
        }
    }
    m
}

#[test]
fn split_cells_get_the_exact_moments_of_each_prong() {
    use crate::clip::clip_solid;
    use crate::moments::rasterise_with;
    let h = 0.25;
    let w = 0.3 * h;
    let model = slotted(w, T);
    let mesh = &model.solids[0];
    let (a, b) = (1.0 - 0.5 * w, 1.0 + 0.5 * w);
    let lower = ([1.0, 0.0, 0.0], [6.0, a, T]);
    let upper = ([1.0, b, 0.0], [6.0, 2.0, T]);
    let mut checked = 0;
    for shift in SHIFTS {
        let options = AnalysisOptions {
            target_cell_size: Some(h),
            grid_shift: shift,
            ..AnalysisOptions::default()
        };
        let grid =
            Grid::for_model([0.0; 3], [6.0, 2.0, T], &options).unwrap_or_else(|e| panic!("{e}"));
        let clipped = clip_solid(&grid, mesh).unwrap_or_else(|e| panic!("{e}"));
        let m = rasterise_with(&grid, mesh.id, &clipped, 2, &|_| false, Splitting::Exact)
            .unwrap_or_else(|e| panic!("{e}"));
        assert_eq!(m.fragments.lattice_moment_cells, 0, "{shift:?}");
        for sc in &m.fragments.split {
            let ijk = grid.cell_ijk(sc.cell).unwrap_or_else(|| panic!("cell"));
            let lo = grid.cell_min(ijk);
            let forms = [lower, upper].map(|(p, q)| box_moments(p, q, lo, h));
            for e in m.entries(sc.cell) {
                let got = &m.cells[e].m;
                let err = |f: &[f64; 27]| {
                    got.iter()
                        .zip(f)
                        .map(|(x, y)| (x - y).abs())
                        .fold(0.0, f64::max)
                };
                let best = err(&forms[0]).min(err(&forms[1]));
                assert!(
                    best <= 1e-12 * h * h * h,
                    "{shift:?} cell {}: {best}",
                    sc.cell
                );
                checked += 1;
            }
        }
    }
    assert!(checked >= 2, "no split cell");
}

/// The worst nodal and interior error of the slotted plate with the
/// linear field u = G x prescribed on every face (relative to |G x|),
/// its split cells, and those left on lattice moments.
fn linear_field_error(shift: [f64; 3], split: Splitting) -> (f64, usize, usize) {
    let w = 0.075;
    let model = slotted(w, T);
    let g = [[1e-3, 2e-3, 0.0], [-1e-3, 5e-4, 1e-3], [0.0, 1e-3, -2e-3]];
    let fixed = FaceMechanics::Support(PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: g,
    });
    let problem = ElasticProblem {
        youngs: 1.0,
        poisson: 0.3,
        faces: vec![fixed; 6],
    };
    let options = AnalysisOptions {
        target_cell_size: Some(0.25),
        grid_shift: shift,
        cg_rel_tol: 1e-12,
        cg_max_iters: 400,
        coarsest_dofs: 600,
        ..AnalysisOptions::default()
    };
    let r = solve_elastic_with(
        &model,
        &problem,
        &options,
        crate::recover::Fallback::InCell,
        None,
        split,
    )
    .unwrap_or_else(|e| panic!("{e}"));
    let mut worst = 0.0f64;
    for (v, d) in model.solids[0]
        .vertices
        .iter()
        .zip(&r.vertex_displacement[0])
    {
        for i in 0..3 {
            let exact: f64 = (0..3).map(|j| g[i][j] * v[j]).sum();
            worst = worst.max((d[i] - exact).abs());
        }
    }
    // Points inside the prongs, by the slot, where cells split.
    for x in [2.03, 3.51, 4.77, 5.6] {
        for y in [1.0 - 0.5 * w - 0.02, 1.0 + 0.5 * w + 0.02] {
            let p = [x, y, 0.37];
            let d = r
                .field
                .displacement_at(p)
                .unwrap_or_else(|e| panic!("{e}"))
                .unwrap_or_else(|| panic!("outside"));
            for i in 0..3 {
                let exact: f64 = (0..3).map(|j| g[i][j] * p[j]).sum();
                worst = worst.max((d[i] - exact).abs());
            }
        }
    }
    (
        worst / 1.2e-2,
        r.report.split_cells,
        r.report.split_cells_lattice_moments,
    )
}

#[test]
fn lf1_is_exact_on_the_slit_plate_and_lattice_moments_are_the_control() {
    // |G x| reaches about 1.2e-2 on the plate. At zero shift the slot
    // is centred on a cell face (no split cell) and the top face lies on
    // a grid plane: its pieces sit in empty cells above the plate and
    // must bond to the plate's dofs (before that fix: 1.5e-3).
    // Measured: exact 3.5e-13 to 4.8e-12 over the shifts with 100-104
    // split cells; lattice moments 7.4e-5 to 2.1e-4.
    let (mut control_worst, mut split) = (0.0f64, 0);
    for shift in SHIFTS {
        let (exact, cells, lattice) = linear_field_error(shift, Splitting::Exact);
        assert!(exact <= 1e-8, "{shift:?}: {exact}");
        assert_eq!(lattice, 0, "{shift:?}");
        split += cells;
        control_worst = control_worst.max(linear_field_error(shift, Splitting::LatticeMoments).0);
    }
    assert!(split > 0, "no split cell");
    assert!(control_worst > 1e-6, "control {control_worst}");
}

#[test]
fn a_sliver_at_the_slot_root_corner_aggregates_to_its_own_prong() {
    // At this shift the cell x 1.12..1.37, y 1.05..1.30 holds a 0.05
    // strip of the upper prong (eta 0.2 < gamma 0.25). It touches the
    // prong cell above through a full face and the root-side cell
    // across the corner only through the strip's end; first-come BFS
    // took the latter (opening 0.920 of converged instead of 0.945).
    use crate::aggregate::aggregate;
    use crate::clip::clip_solid;
    let model = slotted(0.5, T);
    let options = AnalysisOptions {
        target_cell_size: Some(0.25),
        grid_shift: [0.13, 0.2, 0.09],
        ..AnalysisOptions::default()
    };
    let grid = Grid::for_model([0.0; 3], [6.0, 2.0, T], &options).unwrap_or_else(|e| panic!("{e}"));
    let mesh = &model.solids[0];
    let clipped = clip_solid(&grid, mesh).unwrap_or_else(|e| panic!("{e}"));
    let m =
        crate::moments::rasterise_with(&grid, mesh.id, &clipped, 2, &|_| false, Splitting::Exact)
            .unwrap_or_else(|e| panic!("{e}"));
    let area = clipped.pieces.iter().map(|p| p.area).sum();
    let a = aggregate(&grid, &m, 0.25, area).unwrap_or_else(|e| panic!("{e}"));
    let near = |x: f64, v: f64| (x - v).abs() < 0.01;
    let mut slivers = 0;
    for (i, c) in m.cells.iter().enumerate() {
        let lo = grid.cell_min(grid.cell_ijk(c.cell).unwrap_or_else(|| panic!("cell")));
        if !(near(lo[0], 1.12) && near(lo[1], 1.05)) || c.volume_fraction(0.25) >= 0.25 {
            continue;
        }
        let root = m.cells[a.root[i]].cell;
        let r = grid.cell_min(grid.cell_ijk(root).unwrap_or_else(|| panic!("root")));
        // The prong cell above, not the root-side cell at x 0.87.
        assert!(near(r[0], 1.12) && near(r[1], 1.30), "{lo:?} -> {r:?}");
        slivers += 1;
    }
    assert!(slivers >= 4, "{slivers} sliver cells");
}

#[test]
fn a_slot_face_just_inside_a_cell_wall_does_not_weld() {
    // PH-060: at a y shift of 0.2075 the upper prong's face (y = 1.0375)
    // lies 0.005 (0.02 h) below a cell wall, under the lattice spacing:
    // before, that sliver of prong took the lower prong's label and the
    // plate welded (opening 0.0018 of the reference, no split cell).
    // Measured now: 0.986 and 0.909, against the resolved slot's 0.940
    // and 0.945 at the same shifts (CC1's 5%).
    let h = 0.25;
    let w = 0.3 * h;
    let reference = opening(w, w / 1.1, [0.0; 3], Splitting::Exact).0;
    let resolved_reference = opening(0.5, 0.0625, [0.0; 3], Splitting::Exact).0;
    for shift in [[0.0, 0.2075, 0.0], [0.07, 0.2075, 0.03]] {
        let (slit, r) = opening(w, h, shift, Splitting::Exact);
        let slit = slit / reference;
        let resolved = opening(0.5, h, shift, Splitting::Exact).0 / resolved_reference;
        assert!(r.report.split_cells > 0, "{shift:?}: no split cell");
        assert!(
            (slit - resolved).abs() <= 0.05,
            "{shift:?}: {slit} vs {resolved}"
        );
    }
}
