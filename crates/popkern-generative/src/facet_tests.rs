use super::*;

/// A box [0, a] × [0, b] × [0, c] split into an n × n grid of quads per
/// face, two triangles each, counter-clockwise from outside; vertices
/// shared.
fn grid_box(size: V3, n: usize) -> (Vec<V3>, Vec<[u32; 3]>) {
    let mut vertices: Vec<V3> = Vec::new();
    let mut index: HashMap<[i64; 3], u32> = HashMap::new();
    let mut id = |q: [usize; 3], vertices: &mut Vec<V3>| -> u32 {
        let k = q.map(|x| x as i64);
        *index.entry(k).or_insert_with(|| {
            vertices.push([0, 1, 2].map(|d| size[d] * q[d] as f64 / n as f64));
            (vertices.len() - 1) as u32
        })
    };
    let mut triangles = Vec::new();
    // Each face: the fixed axis, its side, and the two running axes in
    // the order that makes (u × v) point outward.
    for axis in 0..3 {
        for side in [0, n] {
            let (u, v) = if side == n {
                ((axis + 1) % 3, (axis + 2) % 3)
            } else {
                ((axis + 2) % 3, (axis + 1) % 3)
            };
            for i in 0..n {
                for j in 0..n {
                    let at = |di: usize, dj: usize| {
                        let mut q = [0; 3];
                        q[axis] = side;
                        q[u] = i + di;
                        q[v] = j + dj;
                        q
                    };
                    let [a, b, c, d] =
                        [at(0, 0), at(1, 0), at(1, 1), at(0, 1)].map(|q| id(q, &mut vertices));
                    triangles.push([a, b, c]);
                    triangles.push([a, c, d]);
                }
            }
        }
    }
    (vertices, triangles)
}

fn volume(vertices: &[V3], triangles: &[[u32; 3]]) -> f64 {
    triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| vertices[i as usize]);
            dot(a, cross(b, c)) / 6.0
        })
        .sum()
}

fn options() -> FacetOptions {
    FacetOptions::for_tolerance(&Tolerance::default())
}

#[test]
fn a_subdivided_box_merges_to_six_faces_and_four_edge_vertices() {
    let (v, t) = grid_box([2.0, 3.0, 5.0], 4);
    assert!(volume(&v, &t) > 29.99);
    let owners = vec![0; t.len()];
    let body =
        facet(&v, &t, &owners, &options(), Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(body.stats.faces, 6);
    assert_eq!(body.brep.faces.len(), 6);
    // Only the 8 corners stay: the edge-midpoint vertices sit on
    // straight edges between two faces; the face-interior ones are gone
    // with the merge.
    assert_eq!(body.brep.vertices.len(), 8);
    assert_eq!(body.brep.edges.len(), 12);
}

#[test]
fn owners_split_faces() {
    let (v, t) = grid_box([1.0, 1.0, 1.0], 2);
    // The +z face's triangles (each face has 8): half owner 1.
    let owners: Vec<u32> = (0..t.len())
        .map(|i| {
            let tri = t[i];
            let top = tri.iter().all(|&k| v[k as usize][2] == 1.0);
            let left = tri.iter().all(|&k| v[k as usize][0] <= 0.5);
            u32::from(top && left)
        })
        .collect();
    let body =
        facet(&v, &t, &owners, &options(), Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(body.stats.faces, 7);
    assert_eq!(body.owners.iter().filter(|&&o| o == 1).count(), 1);
}

#[test]
fn an_open_mesh_declines_where_it_is_open() {
    let (v, mut t) = grid_box([1.0, 1.0, 1.0], 1);
    t.pop();
    let owners = vec![0; t.len()];
    match facet(&v, &t, &owners, &options(), Tolerance::default()) {
        Err(FacetError::NotClosed { uses: 1, .. }) => {}
        other => panic!("{other:?}"),
    }
}

#[test]
fn two_boxes_touching_at_a_corner_decline_as_pinched() {
    let (mut v, mut t) = grid_box([1.0, 1.0, 1.0], 1);
    let (v2, t2) = grid_box([1.0, 1.0, 1.0], 1);
    let base = v.len() as u32;
    // The second box's origin corner sits on the first's (1, 1, 1).
    for p in v2 {
        v.push([p[0] + 1.0, p[1] + 1.0, p[2] + 1.0]);
    }
    t.extend(t2.iter().map(|tri| tri.map(|i| i + base)));
    let owners = vec![0; t.len()];
    match facet(&v, &t, &owners, &options(), Tolerance::default()) {
        Err(FacetError::PinchedVertex { at, .. }) => assert_eq!(at, [1.0, 1.0, 1.0]),
        other => panic!("{other:?}"),
    }
}

#[test]
fn two_cubes_glued_on_a_wall_merge_across_the_seam() {
    let cube = |o: V3| {
        let (v, t) = grid_box([1.0, 1.0, 1.0], 1);
        (
            v.iter()
                .map(|p| [p[0] + o[0], p[1] + o[1], p[2] + o[2]])
                .collect::<Vec<_>>(),
            t,
        )
    };
    let (a, ta) = cube([0.0, 0.0, 0.0]);
    let (b, tb) = cube([1.0, 0.0, 0.0]);
    let mut vertices = a.clone();
    let mut triangles: Vec<[u32; 3]> = Vec::new();
    let mut map_b = Vec::new();
    for p in &b {
        match vertices.iter().position(|q| q == p) {
            Some(i) => map_b.push(i as u32),
            None => {
                vertices.push(*p);
                map_b.push((vertices.len() - 1) as u32);
            }
        }
    }
    let tb: Vec<[u32; 3]> = tb.iter().map(|t| t.map(|i| map_b[i as usize])).collect();
    // Drop the shared x = 1 wall from both.
    let on_wall = |t: &[u32; 3]| t.iter().all(|&i| vertices[i as usize][0] == 1.0);
    triangles.extend(ta.iter().filter(|t| !on_wall(t)));
    triangles.extend(tb.iter().filter(|t| !on_wall(t)));
    let owners = vec![0; triangles.len()];
    let body = facet(
        &vertices,
        &triangles,
        &owners,
        &options(),
        Tolerance::default(),
    )
    .unwrap_or_else(|e| panic!("{e}"));
    // The 2 × 1 × 1 box: six faces, the four x = 1 vertices dropped.
    assert_eq!(body.stats.faces, 6);
    assert_eq!(body.brep.vertices.len(), 8);
}

/// The boundary of unit voxels as a closed mesh, vertices shared.
fn voxels(cells: &[[i64; 3]]) -> (Vec<V3>, Vec<[u32; 3]>) {
    let set: std::collections::HashSet<[i64; 3]> = cells.iter().copied().collect();
    let mut index: HashMap<[i64; 3], u32> = HashMap::new();
    let mut v: Vec<V3> = Vec::new();
    let mut id = |q: [i64; 3], v: &mut Vec<V3>| {
        *index.entry(q).or_insert_with(|| {
            v.push(q.map(|x| x as f64));
            (v.len() - 1) as u32
        })
    };
    let mut t = Vec::new();
    for &c in cells {
        for axis in 0..3 {
            for side in [0i64, 1] {
                let mut n = c;
                n[axis] += if side == 1 { 1 } else { -1 };
                if set.contains(&n) {
                    continue;
                }
                let (u, w) = if side == 1 {
                    ((axis + 1) % 3, (axis + 2) % 3)
                } else {
                    ((axis + 2) % 3, (axis + 1) % 3)
                };
                let at = |du: i64, dw: i64| {
                    let mut q = c;
                    q[axis] += side;
                    q[u] += du;
                    q[w] += dw;
                    q
                };
                let [a, b, cc, d] = [at(0, 0), at(1, 0), at(1, 1), at(0, 1)].map(|q| id(q, &mut v));
                t.push([a, b, cc]);
                t.push([a, cc, d]);
            }
        }
    }
    (v, t)
}

#[test]
fn a_ring_keeps_its_top_as_one_face_with_a_hole() {
    // A 3 × 3 × 1 ring of voxels: top and bottom are square annuli.
    let cells: Vec<[i64; 3]> = (0..3)
        .flat_map(|y| (0..3).map(move |x| [x, y, 0]))
        .filter(|c| *c != [1, 1, 0])
        .collect();
    let (v, t) = voxels(&cells);
    let owners = vec![0; t.len()];
    let body =
        facet(&v, &t, &owners, &options(), Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    // Top, bottom, four outer walls, four inner walls.
    assert_eq!(body.brep.faces.len(), 10);
    let holed = body
        .brep
        .faces
        .values()
        .filter(|f| f.inner_loops.len() == 1)
        .count();
    assert_eq!(holed, 2, "top and bottom keep their hole");
    assert!((body.stats.mesh_volume - 8.0).abs() < 1e-12);
}

#[test]
fn separate_pieces_become_separate_solids() {
    let (v, t) = voxels(&[[0, 0, 0], [3, 0, 0]]);
    let owners = vec![0; t.len()];
    let body =
        facet(&v, &t, &owners, &options(), Tolerance::default()).unwrap_or_else(|e| panic!("{e}"));
    assert_eq!(body.brep.solids.len(), 2);
    assert_eq!(body.brep.faces.len(), 12);
}
