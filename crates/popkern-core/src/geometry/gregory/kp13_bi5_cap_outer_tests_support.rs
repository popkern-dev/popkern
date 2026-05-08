use super::*;

pub(super) const EPS: f64 = 1e-10;

pub(super) fn approx_p(a: Point3, b: Point3, eps: f64) -> bool {
    (a.x - b.x).abs() < eps && (a.y - b.y).abs() < eps && (a.z - b.z).abs() < eps
}

// ---- Round-trip: set_corner_jet_xy → extract jet → original jet ----

pub(super) fn extract_jet_10(piece: &Bi5CapPiece) -> TaylorJet2 {
    let p50 = piece.cps[5][0];
    let p40 = piece.cps[4][0];
    let p51 = piece.cps[5][1];
    let p30 = piece.cps[3][0];
    let p52 = piece.cps[5][2];
    let p41 = piece.cps[4][1];
    TaylorJet2 {
        value: p50,
        d_u: Point3::new(
            5.0 * (p50.x - p40.x),
            5.0 * (p50.y - p40.y),
            5.0 * (p50.z - p40.z),
        ),
        d_v: Point3::new(
            5.0 * (p51.x - p50.x),
            5.0 * (p51.y - p50.y),
            5.0 * (p51.z - p50.z),
        ),
        d_uu: Point3::new(
            20.0 * (p50.x - 2.0 * p40.x + p30.x),
            20.0 * (p50.y - 2.0 * p40.y + p30.y),
            20.0 * (p50.z - 2.0 * p40.z + p30.z),
        ),
        d_vv: Point3::new(
            20.0 * (p52.x - 2.0 * p51.x + p50.x),
            20.0 * (p52.y - 2.0 * p51.y + p50.y),
            20.0 * (p52.z - 2.0 * p51.z + p50.z),
        ),
        d_uv: Point3::new(
            25.0 * (p51.x - p50.x - p41.x + p40.x),
            25.0 * (p51.y - p50.y - p41.y + p40.y),
            25.0 * (p51.z - p50.z - p41.z + p40.z),
        ),
    }
}

pub(super) fn extract_jet_01(piece: &Bi5CapPiece) -> TaylorJet2 {
    let p05 = piece.cps[0][5];
    let p15 = piece.cps[1][5];
    let p04 = piece.cps[0][4];
    let p25 = piece.cps[2][5];
    let p03 = piece.cps[0][3];
    let p14 = piece.cps[1][4];
    TaylorJet2 {
        value: p05,
        d_u: Point3::new(
            5.0 * (p15.x - p05.x),
            5.0 * (p15.y - p05.y),
            5.0 * (p15.z - p05.z),
        ),
        d_v: Point3::new(
            5.0 * (p05.x - p04.x),
            5.0 * (p05.y - p04.y),
            5.0 * (p05.z - p04.z),
        ),
        d_uu: Point3::new(
            20.0 * (p25.x - 2.0 * p15.x + p05.x),
            20.0 * (p25.y - 2.0 * p15.y + p05.y),
            20.0 * (p25.z - 2.0 * p15.z + p05.z),
        ),
        d_vv: Point3::new(
            20.0 * (p05.x - 2.0 * p04.x + p03.x),
            20.0 * (p05.y - 2.0 * p04.y + p03.y),
            20.0 * (p05.z - 2.0 * p04.z + p03.z),
        ),
        d_uv: Point3::new(
            25.0 * (p15.x - p05.x - p14.x + p04.x),
            25.0 * (p15.y - p05.y - p14.y + p04.y),
            25.0 * (p15.z - p05.z - p14.z + p04.z),
        ),
    }
}

pub(super) fn extract_jet_11(piece: &Bi5CapPiece) -> TaylorJet2 {
    let p55 = piece.cps[5][5];
    let p45 = piece.cps[4][5];
    let p54 = piece.cps[5][4];
    let p35 = piece.cps[3][5];
    let p53 = piece.cps[5][3];
    let p44 = piece.cps[4][4];
    TaylorJet2 {
        value: p55,
        d_u: Point3::new(
            5.0 * (p55.x - p45.x),
            5.0 * (p55.y - p45.y),
            5.0 * (p55.z - p45.z),
        ),
        d_v: Point3::new(
            5.0 * (p55.x - p54.x),
            5.0 * (p55.y - p54.y),
            5.0 * (p55.z - p54.z),
        ),
        d_uu: Point3::new(
            20.0 * (p55.x - 2.0 * p45.x + p35.x),
            20.0 * (p55.y - 2.0 * p45.y + p35.y),
            20.0 * (p55.z - 2.0 * p45.z + p35.z),
        ),
        d_vv: Point3::new(
            20.0 * (p55.x - 2.0 * p54.x + p53.x),
            20.0 * (p55.y - 2.0 * p54.y + p53.y),
            20.0 * (p55.z - 2.0 * p54.z + p53.z),
        ),
        d_uv: Point3::new(
            25.0 * (p55.x - p45.x - p54.x + p44.x),
            25.0 * (p55.y - p45.y - p54.y + p44.y),
            25.0 * (p55.z - p45.z - p54.z + p44.z),
        ),
    }
}

pub(super) fn sample_jet(seed: f64) -> TaylorJet2 {
    TaylorJet2 {
        value: Point3::new(seed, 2.0 * seed, -seed),
        d_u: Point3::new(seed * 0.3, 0.5, 1.1),
        d_v: Point3::new(0.2, seed * 0.6, -0.3),
        d_uu: Point3::new(0.1, -0.2, seed * 0.05),
        d_vv: Point3::new(seed * 0.07, 0.3, -0.1),
        d_uv: Point3::new(0.04, 0.08, seed * 0.02),
    }
}
