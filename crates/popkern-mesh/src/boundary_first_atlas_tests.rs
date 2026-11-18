//! Atlas assembly contracts.

use super::*;
use popkern_core::tolerance::Tolerance;

/// A ring's source names the segment leaving its sample: in an assembled
/// piece, a crossing sample's next segment runs along the row and has no
/// edge. Charged to the edge the crossing leaves along (its continuation
/// into the other piece), a split the band requested on a row segment
/// landed on a ball's seam inside the cap, and the slit meridian read a
/// row segment as a real-edge crossing.
#[test]
fn a_segment_along_a_row_carries_no_edge() {
    let brep = BRep::sphere(1.0, Tolerance::default()).unwrap();
    let (_, face) = brep.faces.iter().next().unwrap();
    let on_face = face.surface.as_ref().unwrap();
    let (surface, orientation) = (&on_face.geometry, on_face.orientation);
    let (seam, _) = brep.edges.iter().next().unwrap();
    let descriptor = ChartDescriptor::of(surface, &brep.tolerance);
    let cuts: Vec<Cut> = descriptor
        .singular_rows
        .iter()
        .filter(|row| row.kind == SingularKind::Chart)
        .map(|row| {
            cut_at(
                &descriptor,
                cap_row_at(&descriptor, surface, row, CAP_TURNS[0]).unwrap(),
            )
        })
        .collect();
    assert_eq!(cuts.len(), 2, "a ball has two chart singular rows");
    // The ball's seam pair: down one side from pole to pole across both
    // rows, up the other a period on, each sample on the seam edge.
    let (north, south) = (cuts[0].cap.singular.at, cuts[1].cap.singular.at);
    let rows = [cuts[0].cap.at, 0.5 * (north + south), cuts[1].cap.at];
    let p = cuts[0].period.unwrap();
    let mut vs = vec![north];
    vs.extend(rows);
    vs.push(south);
    let down: Vec<(f64, f64)> = vs.iter().map(|&v| (0.0, v)).collect();
    let up: Vec<(f64, f64)> = vs.iter().rev().map(|&v| (p, v)).collect();
    let uv: Vec<(f64, f64)> = down.into_iter().chain(up).collect();
    let n = uv.len();
    let ring = BoundaryRing {
        pos: uv.iter().map(|&(u, v)| surface.point_at(u, v)).collect(),
        sources: (0..n)
            .map(|i| (seam, i as f64 / n as f64, (i + 1) as f64 / n as f64))
            .collect(),
        uv,
    };
    let sides: Vec<Side> = ring.uv.iter().map(|&uv| side(&cuts, uv)).collect();
    assert!(sides.iter().filter(|s| matches!(s, Side::Row(_))).count() == 4);
    let pieces = assemble(
        &cuts,
        (std::slice::from_ref(&ring), &[sides]),
        (orientation, 0.01, None, &[0.5, 0.5]),
        (surface, &descriptor),
        &|_| false,
    )
    .or_else(|| {
        // The face lies left of its loop: the other winding when this
        // one is the complement.
        let mut ring = ring.clone();
        ring.uv.reverse();
        ring.pos.reverse();
        let sides: Vec<Side> = ring.uv.iter().map(|&uv| side(&cuts, uv)).collect();
        assemble(
            &cuts,
            (std::slice::from_ref(&ring), &[sides]),
            (orientation, 0.01, None, &[0.5, 0.5]),
            (surface, &descriptor),
            &|_| false,
        )
    })
    .expect("the ball's seam ring assembles into pieces");
    assert_eq!(pieces.len(), 3, "two caps and the band");
    let on_row = |uv: (f64, f64)| cuts.iter().any(|c| c.on_row(uv));
    for (piece, rings) in &pieces {
        for r in rings {
            let m = r.uv.len();
            for k in 0..m {
                if on_row(r.uv[k]) && on_row(r.uv[(k + 1) % m]) {
                    assert_eq!(
                        r.sources[k].0,
                        EdgeHandle::default(),
                        "{piece:?}: segment {k} runs along a row but names an edge"
                    );
                }
            }
        }
    }
}
