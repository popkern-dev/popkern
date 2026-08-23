use super::mask::MaskSpec;
use super::types::PolyhedralSplineError;

// The PnS mask tables are popkern's own, regenerated clean-room from the
// published papers (PK-317; generators and per-table provenance live in
// tools/reference/pns-tables/ and crates/popkern-core/data/pns-masks/
// README.md). The ACM ToMS 1032 reference package remains deliberately
// untracked (license — PK-315) and is only consulted, when present, by
// fixture-driven tests gated on cfg(toms_package).
macro_rules! table {
    ($name:literal) => {
        include_str!(concat!("../../../data/pns-masks/", $name))
    };
}

pub const REGULAR_MASK: [[f64; 9]; 9] = [
    [0.25, 0.25, 0.0, 0.25, 0.25, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.25, 0.25, 0.0, 0.25, 0.25, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.25, 0.25, 0.0, 0.25, 0.25, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.5, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.25, 0.25, 0.0, 0.25, 0.25],
];

pub fn extraordinary(valence: usize) -> Result<MaskSpec, PolyhedralSplineError> {
    match valence {
        3 => Ok(spec("eopSct3.csv", 48, 7, table!("eopSct3.csv"))),
        5 => Ok(spec("eopSct5.csv", 80, 11, table!("eopSct5.csv"))),
        6 => Ok(spec("eopSct6.csv", 384, 13, table!("eopSct6.csv"))),
        7 => Ok(spec("eopSct7.csv", 448, 15, table!("eopSct7.csv"))),
        8 => Ok(spec("eopSct8.csv", 512, 17, table!("eopSct8.csv"))),
        _ => Err(PolyhedralSplineError::UnsupportedValence {
            family: "extraordinary star",
            value: valence,
            supported: "3, 5, 6, 7, 8",
        }),
    }
}

pub fn ngon(sides: usize) -> Result<MaskSpec, PolyhedralSplineError> {
    match sides {
        3 => Ok(spec("ngonSct3.csv", 48, 12, table!("ngonSct3.csv"))),
        5 => Ok(spec("ngonSct5.csv", 80, 20, table!("ngonSct5.csv"))),
        6 => Ok(spec("ngonSct6.csv", 384, 24, table!("ngonSct6.csv"))),
        7 => Ok(spec("ngonSct7.csv", 448, 28, table!("ngonSct7.csv"))),
        8 => Ok(spec("ngonSct8.csv", 512, 32, table!("ngonSct8.csv"))),
        _ => Err(PolyhedralSplineError::UnsupportedValence {
            family: "n-gon",
            value: sides,
            supported: "3, 5, 6, 7, 8",
        }),
    }
}

pub fn polar(sectors: usize) -> Result<MaskSpec, PolyhedralSplineError> {
    match sectors {
        3 => Ok(spec("polarSct3.csv", 36, 4, table!("polarSct3.csv"))),
        4 => Ok(spec("polarSct4.csv", 48, 5, table!("polarSct4.csv"))),
        5 => Ok(spec("polarSct5.csv", 60, 6, table!("polarSct5.csv"))),
        6 => Ok(spec("polarSct6.csv", 72, 7, table!("polarSct6.csv"))),
        7 => Ok(spec("polarSct7.csv", 84, 8, table!("polarSct7.csv"))),
        8 => Ok(spec("polarSct8.csv", 96, 9, table!("polarSct8.csv"))),
        _ => Err(PolyhedralSplineError::UnsupportedValence {
            family: "polar",
            value: sectors,
            supported: "3, 4, 5, 6, 7, 8",
        }),
    }
}

pub fn t0() -> Result<MaskSpec, PolyhedralSplineError> {
    Ok(spec("T0.csv", 64, 14, table!("T0.csv")))
}

pub fn t1() -> Result<MaskSpec, PolyhedralSplineError> {
    Ok(spec("T1.csv", 128, 18, table!("T1.csv")))
}

pub fn t2() -> Result<MaskSpec, PolyhedralSplineError> {
    Ok(spec("T2.csv", 256, 20, table!("T2.csv")))
}

fn spec(name: &'static str, rows: usize, cols: usize, text: &'static str) -> MaskSpec {
    MaskSpec {
        name,
        rows,
        cols,
        text,
    }
}
