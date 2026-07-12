//! Tests for the KP15 bi-6 deep-interior population.
//!
//! The primary gate is [`fixture_reproduces_all_deep_block_cps`]:
//! it consumes the committed Maple-script oracle fixture
//! (`tools/reference/kp15-bi6-cap/n6_cap_fixture.json`) — a
//! fixed-seed rational `n = 6` cap with all 16 deep-block CPs per
//! sector and zero-residual provenance — and asserts the Rust port
//! reproduces every deep-block CP to ≤ 1e-12.

use super::*;
use crate::geometry::Point3;

/// Parse a `"num/den"` rational string into `f64`.
fn parse_rat(s: &str) -> f64 {
    let (num, den) = s.split_once('/').expect("rational as num/den");
    num.trim().parse::<f64>().expect("numerator") / den.trim().parse::<f64>().expect("denominator")
}

/// Wrap a scalar oracle value as a `Point3` (x-channel only). The
/// construction is linear, so the x-channel reproduces the scalar
/// oracle exactly while y/z stay zero.
fn scalar(v: f64) -> Point3 {
    Point3::new(v, 0.0, 0.0)
}

fn load_fixture() -> serde_json::Value {
    let path = format!(
        "{}/../../tools/reference/kp15-bi6-cap/n6_cap_fixture.json",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("read {path}: {e}"));
    serde_json::from_str(&text).expect("parse fixture JSON")
}

#[test]
fn fixture_reproduces_all_deep_block_cps() {
    let fx = load_fixture();
    let n = fx["n"].as_u64().unwrap() as usize;
    assert_eq!(n, 6);

    let params = Bi6Reparams::for_valence(n);
    // The fixture is generated at c = 1/2, τ = 17/20, d2 = 103/1156.
    // Bi6Reparams::for_valence(6) reproduces these from the table.
    assert!((params.c - 0.5).abs() < 1e-15);
    assert!((params.tau - parse_rat(fx["tau"].as_str().unwrap())).abs() < 1e-15);
    assert!((params.d2 - parse_rat(fx["d2"].as_str().unwrap())).abs() < 1e-12);

    // ---- free set ----
    let fs = &fx["free_set"];
    let vec_of = |key: &str| -> Vec<Point3> {
        fs[key]
            .as_array()
            .unwrap()
            .iter()
            .map(|v| scalar(parse_rat(v.as_str().unwrap())))
            .collect()
    };
    let free = Bi6DeepInteriorFree {
        p00: scalar(parse_rat(fs["P00"].as_str().unwrap())),
        tan: [
            scalar(parse_rat(fs["TAN0"].as_str().unwrap())),
            scalar(parse_rat(fs["TAN1"].as_str().unwrap())),
        ],
        two: [
            scalar(parse_rat(fs["TWO0"].as_str().unwrap())),
            scalar(parse_rat(fs["TWO1"].as_str().unwrap())),
            scalar(parse_rat(fs["TWO2"].as_str().unwrap())),
        ],
        p220: scalar(parse_rat(fs["P220"].as_str().unwrap())),
        p210: scalar(parse_rat(fs["P210"].as_str().unwrap())),
        p23: vec_of("P23"),
        p33: vec_of("P33"),
    };

    // ---- input net ----
    let net_json = fx["input_net_4x4_per_sector"].as_array().unwrap();
    let net: Vec<[[Point3; 4]; 4]> = net_json
        .iter()
        .map(|blk| {
            let rows = blk.as_array().unwrap();
            let mut m = [[Point3::new(0.0, 0.0, 0.0); 4]; 4];
            for (r, row) in rows.iter().enumerate() {
                for (col, v) in row.as_array().unwrap().iter().enumerate() {
                    m[r][col] = scalar(parse_rat(v.as_str().unwrap()));
                }
            }
            m
        })
        .collect();
    assert_eq!(net.len(), n);

    // ---- build ----
    let grids = build_bi6_cap_from_net(&params, &net, &free).expect("build cap");
    assert_eq!(grids.sectors.len(), n);

    // The deferred lyc[0] residual must vanish (net is C²-consistent).
    assert!(
        grids.sal.x.abs() < 1e-10,
        "sal (deferred eq(8) residual) = {} should vanish",
        grids.sal.x
    );

    // ---- compare every deep-block CP (script i,j in [3,6]) ----
    let deep = fx["deep_block_cp_per_sector"].as_array().unwrap();
    let mut worst = 0.0_f64;
    let mut worst_at = (0usize, 0usize, 0usize);
    for (s, block) in deep.iter().enumerate() {
        for i in 3..7 {
            for j in 3..7 {
                let key = format!("{i},{j}");
                let expected = parse_rat(block[&key].as_str().unwrap());
                let got = grids.sectors[s][i][j].x;
                let err = (got - expected).abs();
                if err > worst {
                    worst = err;
                    worst_at = (s, i, j);
                }
                assert!(
                    err <= 1e-12,
                    "deep CP sector {s} [{i},{j}]: got {got}, expected {expected}, err {err:.3e}"
                );
            }
        }
    }
    // Sanity: y/z channels stayed zero (linearity of the port).
    for grid in grids.sectors.iter() {
        for i in 3..7 {
            for j in 3..7 {
                assert_eq!(grid[i][j].y, 0.0);
                assert_eq!(grid[i][j].z, 0.0);
            }
        }
    }
    println!(
        "KP15 deep-interior fixture: worst deep-block error {worst:.3e} at sector/i/j {worst_at:?}"
    );
}

#[test]
fn declines_low_valence_with_diagnostic() {
    // n < 5 is out of scope for this valence-general path.
    for n in [3usize, 4] {
        let params = Bi6Reparams::for_valence(n);
        let net = vec![[[Point3::new(0.0, 0.0, 0.0); 4]; 4]; n];
        let free = Bi6DeepInteriorFree {
            p00: Point3::new(0.0, 0.0, 0.0),
            tan: [Point3::new(0.0, 0.0, 0.0); 2],
            two: [Point3::new(0.0, 0.0, 0.0); 3],
            p220: Point3::new(0.0, 0.0, 0.0),
            p210: Point3::new(0.0, 0.0, 0.0),
            p23: vec![Point3::new(0.0, 0.0, 0.0); n],
            p33: vec![Point3::new(0.0, 0.0, 0.0); n],
        };
        let err = build_bi6_cap_from_net(&params, &net, &free);
        assert!(err.is_err(), "n = {n} should decline");
    }
}

#[test]
fn length_mismatch_rejected() {
    let params = Bi6Reparams::for_valence(6);
    let net = vec![[[Point3::new(0.0, 0.0, 0.0); 4]; 4]; 5]; // wrong length
    let free = Bi6DeepInteriorFree {
        p00: Point3::new(0.0, 0.0, 0.0),
        tan: [Point3::new(0.0, 0.0, 0.0); 2],
        two: [Point3::new(0.0, 0.0, 0.0); 3],
        p220: Point3::new(0.0, 0.0, 0.0),
        p210: Point3::new(0.0, 0.0, 0.0),
        p23: vec![Point3::new(0.0, 0.0, 0.0); 6],
        p33: vec![Point3::new(0.0, 0.0, 0.0); 6],
    };
    assert!(build_bi6_cap_from_net(&params, &net, &free).is_err());
}
