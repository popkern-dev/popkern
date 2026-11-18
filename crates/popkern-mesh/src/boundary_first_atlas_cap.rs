//! Meshing the atlas's pieces (chart-service.md §6a, steps 5 and 6): the
//! band through the face-chart pipeline, each cap in its pole projection
//! (an isotropic lattice there, the shared constrained triangulation and
//! deviation refinement through `ChartedSurface`), merged into one face
//! mesh with the row samples the pieces share welded.

use popkern_math::PortableMath;
use popkern_math::hash::FastHashMap as HashMap;

use popkern_core::geometry::chart::atlas::cap_chart;
use popkern_core::geometry::chart::ops::{SINGULAR_EXACT, chart_normal};
use popkern_core::geometry::chart::{ChartDescriptor, ChartedSurface};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, FaceHandle};

use super::super::{BoundaryRing, SplitRequest};
use super::{Cut, Piece, Pieces};

#[path = "boundary_first_atlas_cap_rows.rs"]
mod rows;
use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};
use rows::{exposed_row_edge, join_seam_rings, merge, with_band_row_samples};

#[path = "boundary_first_atlas_cap_lattice.rs"]
mod lattice;
use lattice::cap_lattice;
pub(super) use lattice::{band_row_cells, estimated_cap_cell};

/// Mesh every piece and merge them. `None` when a piece declines.
pub(super) fn mesh_pieces(
    brep: &BRep,
    face: FaceHandle,
    (surface, orientation): (&SurfaceGeometry, bool),
    (descriptor, options): (&ChartDescriptor, &TessellationOptions),
    (cuts, pieces, meridian): (&[Cut], Pieces, Option<f64>),
    split_requests: &mut Vec<SplitRequest>,
) -> Option<FaceMesh> {
    let _ = brep;
    let mut meshes = Vec::new();
    // The band first: its slit can add samples where it crosses a row
    // circle, and each cap must carry every sample of the rows it shares
    // with the band (else a T-junction on the row).
    let mut pieces = pieces;
    // The caps the face reaches: the band is sized over its box extended
    // to their singular rows, the face's box (`mesh_rings`'s `sizing`).
    let caps: Vec<usize> = pieces
        .iter()
        .filter_map(|(p, _)| match p {
            Piece::Cap(k) => Some(*k),
            Piece::Band => None,
        })
        .collect();
    pieces.sort_by_key(|(piece, _)| !matches!(piece, Piece::Band));
    let mut band_row_samples: Vec<((f64, f64), Point3)> = Vec::new();
    for (piece, rings) in pieces {
        let rings = match piece {
            Piece::Cap(k) => with_band_row_samples(rings, &cuts[k], &band_row_samples),
            Piece::Band => join_seam_rings(rings),
        };
        let piece_name = piece;
        let before = split_requests.len();
        let mesh = match piece {
            Piece::Band => {
                if popkern_core::trace::trace_enabled("atlas-dump") {
                    for (r, ring) in rings.iter().enumerate() {
                        for (i, (uv, src)) in ring.uv.iter().zip(&ring.sources).enumerate() {
                            eprintln!(
                                "    [atlas-dump] band ring {r} {i} uv=({:.5},{:.5}) src={:?}",
                                uv.0, uv.1, src
                            );
                        }
                    }
                }
                let periods = ChartDescriptor::evaluation_periods_of(surface);
                match super::super::chart_rings(
                    (face, surface, orientation),
                    rings,
                    (periods, meridian, true),
                    options,
                )? {
                    super::super::Charted::Rings(mut rings) => {
                        // The pipeline reads ring 0 as the outer boundary (hole
                        // placement aligns the others to its branch): in a band
                        // piece that is the ring enclosing the most chart area,
                        // whatever order assembly gave them (the turned Viviani
                        // ball's figure-eight lobes came first, landed a period
                        // off the band and its hole meshed as material).
                        let area = |r: &BoundaryRing| {
                            let n = r.uv.len();
                            (0..n)
                                .map(|i| {
                                    let (a, b) = (r.uv[i], r.uv[(i + 1) % n]);
                                    a.0 * b.1 - b.0 * a.1
                                })
                                .sum::<f64>()
                                .abs()
                        };
                        if let Some(outer) = (0..rings.len())
                            .max_by(|&a, &b| area(&rings[a]).total_cmp(&area(&rings[b])))
                        {
                            rings.swap(0, outer);
                        }
                        let sizing = face_box(&rings, cuts, &caps);
                        super::super::mesh_rings(
                            (face, surface, orientation),
                            rings,
                            options,
                            split_requests,
                            descriptor,
                            sizing,
                        )
                    }
                    super::super::Charted::Mesh(fm) => Some(fm),
                }
            }
            Piece::Cap(k) => mesh_cap(
                face,
                (surface, orientation),
                (descriptor, options),
                &cuts[k],
                rings,
                split_requests,
            ),
        };
        // A row has no edge: a split request on one cannot be served.
        let null = popkern_core::topology::EdgeHandle::default();
        let mut i = before;
        while i < split_requests.len() {
            if split_requests[i].0 == null {
                split_requests.remove(i);
            } else {
                i += 1;
            }
        }
        if mesh.is_none() && popkern_core::trace::trace_enabled("atlas") {
            eprintln!("    [atlas] face {face:?}: piece {piece_name:?} declined");
        }
        let mesh = mesh?;
        if matches!(piece, Piece::Band) {
            band_row_samples.extend(
                mesh.vertices
                    .iter()
                    .filter(|v| cuts.iter().any(|c| c.on_row(v.uv)))
                    .map(|v| (v.uv, v.position)),
            );
        }
        meshes.push(mesh);
    }
    let merged = merge(face, meshes);
    // Postcondition: the pieces close up along every row (a row edge is
    // used by the two pieces it separates). A boundary edge lying along a
    // row means a piece was assembled wrong or is missing: decline.
    if let Some(edge) = exposed_row_edge(&merged, cuts) {
        if popkern_core::trace::trace_enabled("atlas") {
            let (a, b) = (&merged.vertices[edge.0], &merged.vertices[edge.1]);
            eprintln!(
                "    [atlas] face {face:?}: row edge {edge:?} exposed after merge: uv {:?} {:?} pos {:?} {:?}",
                a.uv, b.uv, a.position, b.position
            );
        }
        return None;
    }
    Some(merged)
}

/// The band's outer ring box extended across to the singular rows of the
/// caps `caps`: the face's box in its chart.
fn face_box(
    rings: &[BoundaryRing],
    cuts: &[Cut],
    caps: &[usize],
) -> Option<((f64, f64), (f64, f64))> {
    let outer = rings.first()?;
    let (mut lo, mut hi) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for &(u, v) in &outer.uv {
        lo = (lo.0.min(u), lo.1.min(v));
        hi = (hi.0.max(u), hi.1.max(v));
    }
    for &k in caps {
        let at = cuts.get(k)?.cap.singular.at;
        match cuts[k].cap.singular.axis {
            popkern_core::geometry::chart::Axis::V => {
                (lo.1, hi.1) = (lo.1.min(at), hi.1.max(at));
            }
            popkern_core::geometry::chart::Axis::U => {
                (lo.0, hi.0) = (lo.0.min(at), hi.0.max(at));
            }
        }
    }
    (lo.0 <= hi.0 && lo.1 <= hi.1).then_some(((lo.0, hi.0), (lo.1, hi.1)))
}

/// Mesh one cap piece in the pole projection about its singular row.
fn mesh_cap(
    face: FaceHandle,
    (surface, orientation): (&SurfaceGeometry, bool),
    (descriptor, options): (&ChartDescriptor, &TessellationOptions),
    cut: &Cut,
    rings: Vec<BoundaryRing>,
    split_requests: &mut Vec<SplitRequest>,
) -> Option<FaceMesh> {
    let why = |reason: &str| {
        if popkern_core::trace::trace_enabled("atlas") {
            eprintln!("    [atlas] face {face:?}: cap declined: {reason}");
        }
    };
    let Some(chart) = cap_chart(descriptor, surface, &cut.cap) else {
        why("no cap chart");
        return None;
    };
    let Some(centre) = chart.origin() else {
        why("the cap chart has no origin");
        return None;
    };
    let project = |p: &Point3| chart.local_at((0.0, 0.0), p);
    // The rings in the projection: positions are canonical, so their
    // projection is exact; face parameters are kept for the vertices.
    // Each ring sample's face parameters, by its canonical position (the
    // slit dissolution below reorders and drops samples).
    let bits = |p: &Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    let mut uv_at: HashMap<[u64; 3], (f64, f64)> = HashMap::default();
    for r in &rings {
        for (p, &uv) in r.pos.iter().zip(&r.uv) {
            uv_at.entry(bits(p)).or_insert(uv);
        }
    }
    // Coincident double traversals (a seam's two uses from the row to the
    // pole) are slits in the projection, which has no period: dissolved
    // into the genuine boundary components, as the band's rings are.
    let rings_xy: Vec<BoundaryRing> = rings
        .into_iter()
        .flat_map(|r| {
            let ring = BoundaryRing {
                uv: r.pos.iter().map(project).collect(),
                pos: r.pos,
                sources: r.sources,
            };
            super::super::boundary_first_slit::dissolve_slits(ring, (None, None))
        })
        .filter(|r| r.uv.len() >= 3)
        .collect();
    let (mut lo, mut hi) = (
        (f64::INFINITY, f64::INFINITY),
        (f64::NEG_INFINITY, f64::NEG_INFINITY),
    );
    for p in rings_xy.iter().flat_map(|r| r.uv.iter()) {
        lo = (lo.0.min(p.0), lo.1.min(p.1));
        hi = (hi.0.max(p.0), hi.1.max(p.1));
    }
    let extent = (hi.0 - lo.0).max(hi.1 - lo.1);
    if !(extent.is_finite() && extent > 0.0) {
        why("empty extent");
        return None;
    }
    let charted = ChartedSurface::new(surface, chart, ((lo.0, hi.0), (lo.1, hi.1)));
    let reach = rings_xy
        .iter()
        .flat_map(|r| r.uv.iter())
        .map(|p| p.0.phypot(p.1))
        .fold(0.0_f64, f64::max);
    let lattice = cap_lattice(&charted, (lo, hi), reach, options);
    let (n, cell) = (lattice.n, lattice.cell);

    let mut uv_pts: Vec<(f64, f64)> = Vec::new();
    let mut pos_pts: Vec<Option<Point3>> = Vec::new();
    let mut ring_ranges = Vec::new();
    for ring in &rings_xy {
        let start = uv_pts.len();
        uv_pts.extend(ring.uv.iter().copied());
        pos_pts.extend(ring.pos.iter().copied().map(Some));
        ring_ranges.push((start, uv_pts.len()));
    }
    // Ring segments longer than one and a half cells request a midpoint
    // split, as the band's do (a row segment is sampled to the cap's cell
    // by assembly and has no edge to split); lattice points keep the
    // band's clearance, a third of a cell.
    // All of a segment's splits in one round (a midpoint a round converged
    // a cap's rim over many encroachment rounds, each re-meshing the face),
    // and none the driver's dedup would refuse (requested again every round).
    let dedup = crate::encroachment::SAMPLE_DEDUP;
    for ring in &rings_xy {
        let m = ring.uv.len();
        for k in 0..m {
            let (edge, f0, f1) = ring.sources[k];
            // In space: the cell is arc length (the lattice's warp).
            let length = ring.pos[k].distance_to(&ring.pos[(k + 1) % m]);
            if length <= 1.5 * cell || edge == popkern_core::topology::EdgeHandle::default() {
                continue;
            }
            let pieces = (length / cell).ceil().max(2.0) as usize;
            if (f1 - f0).abs() / pieces as f64 <= 2.0 * dedup {
                continue;
            }
            for q in 1..pieces {
                split_requests.push((edge, f0 + (f1 - f0) * q as f64 / pieces as f64));
            }
        }
    }
    let ring_index = super::super::boundary_first_ring_index::RingIndex::new(
        &rings_xy,
        1.0,
        1.0,
        (lattice.rim_cell, lattice.rim_cell),
    );
    // The singular point itself, when the piece holds it: every triangle
    // round it then has it as a corner, whose face parameters stand for
    // the whole pole row (the measuring quadrature's pole corner), and no
    // triangle wraps round it in the face chart.
    let holds_origin = ring_index.inside((0.0, 0.0)) && !ring_index.within_clearance((0.0, 0.0));
    if holds_origin {
        uv_pts.push((0.0, 0.0));
        pos_pts.push(Some(centre));
    }
    // Stations offset half a cell from the box's corners, which are ring
    // samples: a corner-aligned lattice put points exactly on a diagonal
    // ring segment (a corner ball's row chord), T-points the constrained
    // recovery cannot resolve.
    for i in 0..n {
        for j in 0..n {
            let p = lattice.station(i, j);
            let clear_of_origin = !holds_origin || p.0.phypot(p.1) > cell / 3.0;
            if clear_of_origin && ring_index.inside(p) && !ring_index.within_clearance(p) {
                uv_pts.push(p);
                pos_pts.push(None);
            }
        }
    }
    let Some((mut triangles, alias)) = super::super::boundary_first_triangulate::triangulate_region(
        face,
        &uv_pts,
        (1.0, 0.0, 1.0),
        &ring_index,
        &ring_ranges,
    ) else {
        if popkern_core::trace::trace_enabled("atlas-dump") {
            for (r, ring) in rings_xy.iter().enumerate() {
                for (i, (xy, src)) in ring.uv.iter().zip(&ring.sources).enumerate() {
                    eprintln!(
                        "    [atlas-dump] tri ring {r} {i} xy=({:.9},{:.9}) src={:?}",
                        xy.0, xy.1, src
                    );
                }
            }
        }
        why("triangulation");
        return None;
    };
    let deviation = super::super::boundary_first_deviation::Chart {
        face,
        surface: &charted,
        planar: false,
        scale: (1.0, 1.0),
        tri_scale: (1.0, 0.0, 1.0),
        ring_index: &ring_index,
        ring_ranges: &ring_ranges,
    };
    super::super::boundary_first_deviation::refine_by_deviation(
        &deviation,
        &mut uv_pts,
        &mut pos_pts,
        &mut triangles,
        options.chord_height,
    );
    // Coverage: every ring segment an edge, else request the encroached
    // span's midpoint (an edge's; a row segment declines).
    let mut edges = popkern_math::hash::FastHashSet::default();
    for t in &triangles {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            edges.insert((x.min(y), x.max(y)));
        }
    }
    let sources: Vec<_> = rings_xy
        .iter()
        .flat_map(|r| r.sources.iter().copied())
        .collect();
    let mut uncovered = false;
    for &(start, end) in &ring_ranges {
        let m = end - start;
        for k in 0..m {
            let (a, b) = (start + k, start + (k + 1) % m);
            let (ca, cb) = (alias[a], alias[b]);
            if ca != cb && !edges.contains(&(ca.min(cb), ca.max(cb))) {
                uncovered = true;
                let (edge, f0, f1) = sources[a];
                split_requests.push((edge, 0.5 * (f0 + f1)));
            }
        }
    }
    if uncovered {
        if popkern_core::trace::trace_enabled("atlas-dump") {
            for (r, ring) in rings_xy.iter().enumerate() {
                for (i, (xy, src)) in ring.uv.iter().zip(&ring.sources).enumerate() {
                    eprintln!(
                        "    [atlas-dump] ring {r} {i} xy=({:.6},{:.6}) src={:?}",
                        xy.0, xy.1, src
                    );
                }
            }
            eprintln!("    [atlas-dump] interior {} cell {cell}", uv_pts.len());
        }
        why("a ring segment uncovered");
        return None;
    }
    // Vertices: ring samples keep their face parameters and positions;
    // interior points take theirs through the chart.
    let exact = |p: &Point3| SINGULAR_EXACT * p.x.abs().max(p.y.abs()).max(p.z.abs()).max(1.0);
    let n_ring = ring_ranges.last().map_or(0, |r| r.1);
    let mut vertices = Vec::with_capacity(uv_pts.len());
    for (i, &xy) in uv_pts.iter().enumerate() {
        let known = (i < n_ring)
            .then(|| pos_pts[i])
            .flatten()
            .and_then(|p| uv_at.get(&bits(&p)).map(|&uv| (uv, p)));
        let (uv, position) = match known {
            Some(k) => k,
            None => {
                let Some(uv) = charted.face_uv(xy.0, xy.1) else {
                    why(&format!("chart inverse at {xy:?} did not converge"));
                    return None;
                };
                (
                    uv,
                    pos_pts[i].unwrap_or_else(|| charted.point_at(xy.0, xy.1)),
                )
            }
        };
        let n = chart_normal(descriptor, surface, uv, position, exact(&position))
            .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
        let normal = if orientation { n } else { -n };
        vertices.push(MeshVertex {
            position,
            normal,
            uv,
        });
    }
    if !orientation {
        for t in &mut triangles {
            t.swap(1, 2);
        }
    }
    Some(FaceMesh {
        face,
        vertices,
        triangles,
    })
}
