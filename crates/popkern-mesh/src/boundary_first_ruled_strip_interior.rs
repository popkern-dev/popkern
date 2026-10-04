//! Interior-resolution check for doubly curved four-sided charts on
//! the structured tensor path (split from `boundary_first_ruled_strip`
//! under the 500-line budget; PK-650).

use popkern_core::geometry::Point3;

/// Whether consecutive cross-chain stations (in a) and rail stations
/// (in b) keep the surface's iso-line chords within `chord_height` at
/// the chart's middle — where a doubly curved surface may bend more
/// than along its boundary.
pub(super) fn resolves_interior(
    chart: &[(f64, f64)],
    crosses: [&[usize]; 2],
    rail: &[usize],
    eval: &dyn Fn(f64, f64) -> Point3,
    chord_height: f64,
) -> bool {
    // The tensor lattice's columns are the UNION of both cross chains'
    // stations (zipped when misaligned), so that is what must resolve
    // the interior (a torus band's smaller rim alone under-samples the
    // band's middle iso-line).
    let mut stations: Vec<f64> = crosses
        .iter()
        .flat_map(|c| c.iter().map(|&i| chart[i].0))
        .collect();
    stations.sort_by(f64::total_cmp);
    stations.dedup_by(|a, b| (*a - *b).abs() <= 1e-12 * (1.0 + b.abs()));
    let cross = crosses[0];
    let a_mid = 0.5 * (chart[cross[0]].0 + chart[cross[cross.len() - 1]].0);
    let b_mid = 0.5 * (chart[rail[0]].1 + chart[rail[rail.len() - 1]].1);
    let sag = |p: Point3, m: Point3, q: Point3| {
        m.distance_to(&Point3::new(
            0.5 * (p.x + q.x),
            0.5 * (p.y + q.y),
            0.5 * (p.z + q.z),
        ))
    };
    let along_a = stations.windows(2).all(|w| {
        let (a0, a1) = (w[0], w[1]);
        sag(
            eval(a0, b_mid),
            eval(0.5 * (a0 + a1), b_mid),
            eval(a1, b_mid),
        ) <= chord_height
    });
    let along_b = rail.windows(2).all(|w| {
        let (b0, b1) = (chart[w[0]].1, chart[w[1]].1);
        sag(
            eval(a_mid, b0),
            eval(a_mid, 0.5 * (b0 + b1)),
            eval(a_mid, b1),
        ) <= chord_height
    });
    along_a && along_b
}
