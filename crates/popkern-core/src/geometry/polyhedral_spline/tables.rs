use super::mask::MaskSpec;
use super::types::PolyhedralSplineError;

// The PnS mask tables come from the ACM ToMS 1032 package, which is
// deliberately NOT in the repository (license — PK-315). When the
// package is present at build time (cfg emitted by build.rs) the
// tables embed exactly as before; when absent, every table lookup
// returns a structured `ToMsTablesUnavailable` error instead of the
// whole crate failing to compile (measured: every CI leg died at
// `include_str!` for over a week before this gate, PK-314).
#[cfg(toms_package)]
macro_rules! table {
    ($name:literal) => {
        include_str!(concat!(
            "../../../../../docs/code/ACM-ToMS-PnS-1032/1032/Src/Source/src/Patch/Table/",
            $name
        ))
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

#[cfg(not(toms_package))]
fn tables_unavailable(family: &'static str) -> PolyhedralSplineError {
    PolyhedralSplineError::ToMsTablesUnavailable { family }
}

#[cfg(not(toms_package))]
pub fn extraordinary(_valence: usize) -> Result<MaskSpec, PolyhedralSplineError> {
    Err(tables_unavailable("extraordinary star"))
}

#[cfg(not(toms_package))]
pub fn ngon(_sides: usize) -> Result<MaskSpec, PolyhedralSplineError> {
    Err(tables_unavailable("n-gon"))
}

#[cfg(not(toms_package))]
pub fn polar(_sectors: usize) -> Result<MaskSpec, PolyhedralSplineError> {
    Err(tables_unavailable("polar"))
}

#[cfg(not(toms_package))]
pub fn t0() -> Result<MaskSpec, PolyhedralSplineError> {
    Err(tables_unavailable("T0"))
}

#[cfg(not(toms_package))]
pub fn t1() -> Result<MaskSpec, PolyhedralSplineError> {
    Err(tables_unavailable("T1"))
}

#[cfg(not(toms_package))]
pub fn t2() -> Result<MaskSpec, PolyhedralSplineError> {
    Err(tables_unavailable("T2"))
}

#[cfg(toms_package)]
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

#[cfg(toms_package)]
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

#[cfg(toms_package)]
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

#[cfg(toms_package)]
pub fn t0() -> Result<MaskSpec, PolyhedralSplineError> {
    Ok(spec("T0.csv", 64, 14, table!("T0.csv")))
}

#[cfg(toms_package)]
pub fn t1() -> Result<MaskSpec, PolyhedralSplineError> {
    Ok(spec("T1.csv", 128, 18, table!("T1.csv")))
}

#[cfg(toms_package)]
pub fn t2() -> Result<MaskSpec, PolyhedralSplineError> {
    Ok(spec("T2.csv", 256, 20, table!("T2.csv")))
}

#[cfg(toms_package)]
fn spec(name: &'static str, rows: usize, cols: usize, text: &'static str) -> MaskSpec {
    MaskSpec {
        name,
        rows,
        cols,
        text,
    }
}
