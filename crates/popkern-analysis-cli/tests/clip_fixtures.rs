//! The facet clipper on the frozen STEP fixtures: per-face piece areas
//! equal the faces' tessellated areas, at several cell sizes and grid
//! offsets.

use std::path::Path;

use popkern_analysis::clip::clip_solid;
use popkern_analysis::grid::Grid;
use popkern_analysis_cli::step_model::model_from_step_file;

#[test]
fn piece_areas_match_face_areas_on_both_fixtures() {
    for name in ["block_two_bores", "block_dimple"] {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../tests/fixtures/analysis/{name}.step"));
        let m = model_from_step_file(&path).unwrap_or_else(|e| panic!("{name}: {e}"));
        let s = &m.solids[0];
        let mut lo = [f64::INFINITY; 3];
        let mut hi = [f64::NEG_INFINITY; 3];
        for v in &s.vertices {
            for a in 0..3 {
                lo[a] = lo[a].min(v[a]);
                hi[a] = hi[a].max(v[a]);
            }
        }
        for (h, shift) in [(0.5, 0.0), (0.77, 0.13), (1.3, 0.61)] {
            let g = Grid::covering(lo.map(|c| c - shift), hi, h).unwrap_or_else(|e| panic!("{e}"));
            let c = clip_solid(&g, s).unwrap_or_else(|e| panic!("{e}"));
            let mut areas = vec![0.0; m.faces.len()];
            for p in &c.pieces {
                areas[p.face.0 as usize] += p.area;
            }
            for (f, face) in m.faces.iter().enumerate() {
                let rel = (areas[f] - face.area).abs() / face.area;
                assert!(
                    rel <= 1e-12,
                    "{name} h {h} face {f}: {} vs {} ({rel:e})",
                    areas[f],
                    face.area
                );
            }
        }
    }
}

/// Gate X5 on the STEP fixtures: rasterised volume equals the mesh's
/// enclosed volume to 1e-12 (exact half) and the analytic volume within
/// chord x curved area (chord half), at several grids.
#[test]
fn x5_volume_parity_on_both_fixtures() {
    use popkern_analysis::model::SolidId;
    use popkern_analysis::moments::{rasterise_solid, volume_parity};
    use std::f64::consts::PI;
    let cases = [
        (
            "block_two_bores",
            1000.0 - 2.0 * PI * 1.5 * 1.5 * 10.0,
            2.0 * 2.0 * PI * 1.5 * 10.0,
        ),
        (
            "block_dimple",
            1000.0 - PI * 4.0 * 7.0 / 3.0,
            2.0 * PI * 3.0 * 2.0,
        ),
    ];
    for (name, exact, curved_area) in cases {
        let path = Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../tests/fixtures/analysis/{name}.step"));
        let m = model_from_step_file(&path).unwrap_or_else(|e| panic!("{name}: {e}"));
        let s = &m.solids[0];
        let mesh_volume = s.enclosed_volume();
        for (h, shift) in [(0.5, 0.0), (0.77, 0.13), (1.3, 0.61)] {
            let g = Grid::covering([-shift; 3], [10.0; 3], h).unwrap_or_else(|e| panic!("{e}"));
            let c = clip_solid(&g, s).unwrap_or_else(|e| panic!("{e}"));
            let r = rasterise_solid(&g, SolidId(0), &c, 2).unwrap_or_else(|e| panic!("{e}"));
            assert_eq!(
                volume_parity(&r, mesh_volume, 1e-12),
                Ok(()),
                "{name} h {h}"
            );
            let bound = m.chord_tolerance * curved_area;
            assert!(
                (r.volume - exact).abs() <= bound,
                "{name} h {h}: {} vs {exact}",
                r.volume
            );
        }
    }
}
