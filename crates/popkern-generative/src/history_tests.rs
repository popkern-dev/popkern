//! The history's record schedule, its blended frames and its bytes.

use super::*;
use crate::history_store::{decode, encode};

fn volume(mesh: &DisplayMesh) -> f64 {
    let p = |i: u32| {
        let i = 3 * i as usize;
        [0, 1, 2].map(|d| f64::from(mesh.positions[i + d]))
    };
    mesh.indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(p);
            (a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

/// Every update of the first dozen, then about 15% apart: 33 records
/// before the final design of a 300-update solve, none more than 15%
/// (plus one) past the last.
#[test]
fn records_are_dense_early_and_geometric_after() {
    let mut recorded = vec![1];
    while let Some(&last) = recorded.last() {
        let next = next_recorded(last);
        if next >= 300 {
            break;
        }
        recorded.push(next);
    }
    assert_eq!(recorded[..12], (1..=12).collect::<Vec<_>>()[..]);
    assert_eq!(recorded.len(), 33, "{recorded:?}");
    for w in recorded.windows(2).skip(11) {
        assert!(w[1] > w[0] && w[1] <= w[0] + w[0] * 3 / 20 + 1, "{w:?}");
    }
}

/// A 2 × 2 × 2 block of cells whose material goes from `from` to `to`.
fn cube_history(from: u8, to: u8) -> DesignHistory {
    let dims = [4, 4, 4];
    let mut active = Vec::new();
    for z in 1..3u32 {
        for y in 1..3u32 {
            for x in 1..3u32 {
                active.push(x + 4 * (y + 4 * z));
            }
        }
    }
    let frame = |update, m| HistoryFrame {
        update,
        compliance_nmm: 1.0 / update as f64,
        volume_fraction: f64::from(m) / 255.0,
        material: vec![m; 8],
    };
    DesignHistory {
        cell: 0.5,
        corner: [0.0; 3],
        dims,
        active,
        updates: 9,
        frames: vec![frame(1, from), frame(9, to)],
        space: display_mesh(
            &[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            &[[0, 1, 2]],
        ),
    }
}

/// A blend of two records grows with its weight: empty below ½, the
/// full block's surface at the second record, and never past either
/// end (positions are clamped; a non-finite one is the first record).
#[test]
fn frames_blend_between_records_and_clamp() {
    let history = cube_history(0, 255);
    assert!(history.frame(0.0).indices.is_empty());
    assert!(history.frame(0.4).indices.is_empty());
    let volumes: Vec<f64> = [0.6, 0.8, 1.0]
        .map(|at| volume(&history.frame(at)))
        .to_vec();
    assert!(volumes.windows(2).all(|w| w[1] > w[0]), "{volumes:?}");
    assert_eq!(history.frame(7.0), history.frame(1.0));
    assert_eq!(history.frame(f64::NAN), history.frame(0.0));
    assert!(volumes[2] > 0.0);
    // Normals are unit and point out of the material.
    let full = history.frame(1.0);
    for n in full.normals.as_chunks::<3>().0 {
        let len = n.iter().map(|x| x * x).sum::<f32>().sqrt();
        assert!((len - 1.0).abs() < 1e-5, "{n:?}");
    }
    assert_eq!(
        DesignHistory {
            frames: Vec::new(),
            ..history
        }
        .frame(0.0),
        DisplayMesh::default()
    );
}

/// Histories survive their bytes, and bytes that are not one, cut
/// short or naming cells outside the grid are refused.
#[test]
fn histories_round_trip_and_bad_bytes_are_refused() {
    let a = cube_history(0, 255);
    let b = cube_history(40, 200);
    let bytes = encode(&[("a".into(), &a), ("b".into(), &b)]);
    assert_eq!(bytes.len() % 4, 0);
    let back = decode(&bytes).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(
        back,
        vec![("a".to_string(), a.clone()), ("b".to_string(), b)]
    );
    assert!(decode(&bytes[..bytes.len() - 5]).is_err());
    assert!(decode(b"PKBAKE01\0\0\0\0").is_err());
    let mut outside = a.clone();
    outside.active[0] = 64;
    assert!(decode(&encode(&[("a".into(), &outside)])).is_err());
    let mut dangling = a;
    dangling.space.indices[0] = 3;
    assert!(decode(&encode(&[("a".into(), &dangling)])).is_err());
}
