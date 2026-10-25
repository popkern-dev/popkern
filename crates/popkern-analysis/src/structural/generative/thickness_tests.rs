//! Local thickness: slabs read their thickness, a block's corners are a
//! small share, a plate thinner than the length is all thin.

use super::*;
use crate::structural::generative::DesignEntry;

/// The cells of a solid box `size` (mm) at cell size `h`.
fn block(size: [f64; 3], h: f64) -> CellField {
    let n = size.map(|l| (l / h).round() as usize);
    let mut entries = Vec::new();
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
            }
        }
    }
    let ones = vec![1.0; entries.len()];
    CellField::new(&entries, &ones, h).unwrap_or_else(|| panic!("no cells"))
}

/// A plate 1, 2 or 3 cells thick, the odd and the even, reads its
/// thickness everywhere under a length above it, and is all thin.
#[test]
fn a_slab_reads_its_thickness() {
    let h = 0.5;
    for n in [1.0, 2.0, 3.0] {
        let t = n * h;
        let thin = thin_region(&block([8.0, 8.0, t], h), t + 0.1).unwrap_or_else(|| panic!());
        assert!((thin.fraction - 1.0).abs() < 1e-12, "{n}: {thin:?}");
        assert!((thin.thickness - t).abs() < 1e-12, "{n}: {thin:?}");
    }
}

/// The 12 × 2 × 1 plate under a 1.5 mm length (a filter radius of
/// 0.75 mm) is all thin, at 1 mm; the 12 × 4 × 2 block (the negative
/// control) has no thin volume, its corners and edges included, at
/// coarse and fine cells.
#[test]
fn a_plate_is_thin_and_a_block_is_not() {
    for h in [1.0, 0.5, 0.25] {
        let plate = thin_region(&block([12.0, 2.0, 1.0], h), 1.5).unwrap_or_else(|| panic!());
        assert!((plate.fraction - 1.0).abs() < 1e-12, "{h}: {plate:?}");
        assert!((plate.thickness - 1.0).abs() < 1e-12, "{h}: {plate:?}");
        let thick = thin_region(&block([12.0, 4.0, 2.0], h), 1.5).unwrap_or_else(|| panic!());
        assert_eq!(thick.fraction, 0.0, "{h}: {thick:?}");
    }
}

/// A design on 8 × 8 × 10 cells of 0.5 mm: a slab `n` cells thick at
/// density 0.9, grey cells (0.4) on both its faces, void beyond.
fn slab_design(n: usize) -> CellField {
    let h = 0.5;
    let mut entries = Vec::new();
    let mut rho = Vec::new();
    for k in 0..10 {
        for j in 0..8 {
            for i in 0..8 {
                let corner = [i, j, k].map(|q| q as f64 * h);
                entries.push(DesignEntry {
                    solid: 0,
                    entry: entries.len(),
                    centre: corner.map(|c| c + 0.5 * h),
                    volume: h * h * h,
                    corner,
                });
                rho.push(if (3..3 + n).contains(&k) {
                    0.9
                } else if k == 2 || k == 3 + n {
                    0.4
                } else {
                    0.0
                });
            }
        }
    }
    CellField::new(&entries, &rho, h).unwrap_or_else(|| panic!("no cells"))
}

/// PH-128: a member's thickness inside a crisp design with grey cells at
/// its faces: a slab 2 cells thick reads 2 cells (1 mm) at every low
/// quantile, edges and corners included, and 4 cells reads 4.
#[test]
fn a_slab_in_a_design_reads_its_thickness() {
    for n in [2, 4] {
        let cells = slab_design(n);
        for share in [0.02, 0.1, 0.5] {
            let t = thickness_quantile(&cells, 3.0, share).unwrap_or(0.0);
            assert!(
                (t - 0.5 * n as f64).abs() < 1e-12,
                "{n} cells, {share}: {t}"
            );
        }
    }
}
