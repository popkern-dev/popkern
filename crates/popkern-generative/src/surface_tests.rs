use super::*;
use popkern_analysis::options::AnalysisOptions;
use popkern_analysis::report::SolveReport;
use popkern_analysis::structural::generative::DesignEntry;

/// A design on unit cells at `cells` (all inside the part), the cells in
/// `solid` at density 0.9, the rest 0.1 except `dense` at 0.4.
fn design(cells: &[[i64; 3]], solid: &[[i64; 3]], dense: &[[i64; 3]]) -> Design {
    let entries: Vec<DesignEntry> = cells
        .iter()
        .enumerate()
        .map(|(i, c)| {
            let corner = c.map(|x| x as f64);
            DesignEntry {
                solid: 0,
                entry: i,
                centre: corner.map(|x| x + 0.5),
                volume: 1.0,
                corner,
            }
        })
        .collect();
    let densities = cells
        .iter()
        .map(|c| {
            if solid.contains(c) {
                0.9
            } else if dense.contains(c) {
                0.4
            } else {
                0.1
            }
        })
        .collect();
    Design {
        floors: vec![0.0; entries.len()],
        entries,
        densities,
        cell_size: 1.0,
        compliance: 0.0,
        volume_fraction: 0.0,
        history: Vec::new(),
        converged: true,
        unsupported_fraction: None,
        grey_fraction: 0.0,
        member_thickness: 0.0,
        eroded_compliance: None,
        report: SolveReport::new(AnalysisOptions::default()),
    }
}

/// Every cell of an n × n × n block.
fn block(n: i64) -> Vec<[i64; 3]> {
    (0..n)
        .flat_map(|z| (0..n).flat_map(move |y| (0..n).map(move |x| [x, y, z])))
        .collect()
}

fn shells(mesh: &SurfaceMesh) -> usize {
    // Edge-connected triangle sets.
    let n = mesh.triangles.len();
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], i: usize) -> usize {
        let mut i = i;
        while p[i] != i {
            p[i] = p[p[i]];
            i = p[i];
        }
        i
    }
    let mut by_edge: std::collections::HashMap<(u32, u32), usize> = Default::default();
    for (t, tri) in mesh.triangles.iter().enumerate() {
        for k in 0..3 {
            let (a, b) = (tri[k], tri[(k + 1) % 3]);
            if let Some(&o) = by_edge.get(&(a.min(b), a.max(b))) {
                let (x, y) = (find(&mut parent, t), find(&mut parent, o));
                parent[x] = y;
            } else {
                by_edge.insert((a.min(b), a.max(b)), t);
            }
        }
    }
    (0..n).filter(|&i| find(&mut parent, i) == i).count()
}

#[test]
fn an_edge_contact_closes_through_its_densest_cell() {
    // (1,1,1) and (2,2,1) touch along an edge; (2,1,1) is the denser of
    // the two cells between them.
    let cells = block(4);
    let d = design(&cells, &[[1, 1, 1], [2, 2, 1]], &[[2, 1, 1]]);
    let body = body_surface(&d);
    assert_eq!((body.closed_cells, body.open_contacts), (1, 0));
    assert!(body.mesh.is_closed());
    assert_eq!(shells(&body.mesh), 1);
    // The same body as the three cells made solid outright.
    let explicit = body_surface(&design(&cells, &[[1, 1, 1], [2, 1, 1], [2, 2, 1]], &[]));
    assert_eq!(body.mesh, explicit.mesh);
}

#[test]
fn a_corner_contact_closes_through_a_face_path() {
    let cells = block(4);
    let d = design(&cells, &[[1, 1, 1], [2, 2, 2]], &[]);
    let body = body_surface(&d);
    assert_eq!((body.closed_cells, body.open_contacts), (2, 0));
    assert_eq!(shells(&body.mesh), 1);
    // Unclosed, the pair meets in a neck along the diagonal: the closed
    // body holds two cells' worth more.
    let crisp = Design {
        densities: d
            .densities
            .iter()
            .map(|&x| if x >= 0.5 { 1.0 } else { 0.0 })
            .collect(),
        ..d.clone()
    };
    let unclosed = design_surface(&crisp, 0.5).volume();
    assert!(
        body.mesh.volume() > unclosed + 1.0,
        "{} vs {unclosed}",
        body.mesh.volume()
    );
}

#[test]
fn face_connected_designs_are_left_alone() {
    let cells = block(4);
    let d = design(&cells, &[[1, 1, 1], [2, 1, 1], [2, 2, 1]], &[]);
    let body = body_surface(&d);
    assert_eq!((body.closed_cells, body.open_contacts), (0, 0));
    assert_eq!(shells(&body.mesh), 1);
}

#[test]
fn a_contact_whose_cells_are_outside_the_part_stays_open() {
    // Only the two touching cells are in the part.
    let cells = [[1, 1, 1], [2, 2, 1]];
    let d = design(&cells, &cells, &[]);
    let body = body_surface(&d);
    assert_eq!((body.closed_cells, body.open_contacts), (0, 1));
}
