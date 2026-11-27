//! Meshing a cone face in its development (split from
//! `boundary_first_cone_development.rs` under the file-size rule): the
//! rings developed from their canonical positions, the shared
//! constrained triangulation and deviation refinement, and the face's
//! parameters back on every vertex.

use popkern_math::PortableMath;
use popkern_math::hash::FastHashMap as HashMap;

use popkern_core::geometry::chart::ChartDescriptor;
use popkern_core::geometry::chart::ops::{SINGULAR_EXACT, chart_normal};
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{EdgeHandle, FaceHandle};

use super::super::{BoundaryRing, SplitRequest};
use super::Development;
use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};

/// Mesh the face's rings in the development. `None` when the
/// triangulation declines or a ring segment stays uncovered (its span's
/// midpoint requested, as the face chart does).
pub(in super::super) fn mesh_developed(
    face: FaceHandle,
    (surface, orientation): (&SurfaceGeometry, bool),
    (descriptor, options): (&ChartDescriptor, &TessellationOptions),
    (dev, rings, azimuths): (&Development<'_>, &[BoundaryRing], &[Vec<f64>]),
    split_requests: &mut Vec<SplitRequest>,
) -> Option<FaceMesh> {
    let why = |reason: &str| {
        if popkern_core::trace::trace_enabled("boundary-first") {
            eprintln!("    [bf] face {face:?}: cone development declined: {reason}");
        }
    };
    let bits = |p: &Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    let mut uv_at: HashMap<[u64; 3], (f64, f64)> = HashMap::default();
    for r in rings {
        for (p, &uv) in r.pos.iter().zip(&r.uv) {
            uv_at.entry(bits(p)).or_insert(uv);
        }
    }
    // Ring samples in the development, read from their canonical
    // positions (the unwrapped azimuth picks the representative): samples
    // the face chart spread along the apex row at one position are one
    // point here, merged.
    let rings_xy: Vec<BoundaryRing> = rings
        .iter()
        .zip(azimuths)
        .flat_map(|(r, az)| {
            let ring = merge_coincident(BoundaryRing {
                uv: r
                    .pos
                    .iter()
                    .zip(az)
                    .map(|(p, &u)| dev.develop_point(p, u))
                    .collect(),
                pos: r.pos.clone(),
                sources: r.sources.clone(),
            });
            // Every apex sample lands on the origin: a ring visiting the
            // apex twice pinches there, its lobes meshed apart.
            super::super::boundary_first_slit::dissolve_slits(ring, (None, None))
                .into_iter()
                .flat_map(super::super::boundary_first_slit::dissolve_pinches)
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
    // The lattice cell: the chord a facet may span at the face's far
    // side, where the cone's normal curvature radius is `s tan α`
    // (smaller toward the apex, where the deviation refinement adds
    // points), within the face's subdivision bounds.
    let reach = rings_xy
        .iter()
        .flat_map(|r| r.uv.iter())
        .map(|p| p.0.phypot(p.1))
        .fold(0.0_f64, f64::max);
    let radius = reach * dev.sin_a / dev.cos_a;
    let chord_cell = (8.0 * options.chord_height.max(0.0) * radius).sqrt();
    let subdivisions = (extent / chord_cell.max(f64::MIN_POSITIVE)).clamp(
        options.min_subdivisions.max(1) as f64,
        options.max_subdivisions.max(1) as f64,
    );
    let cell = extent / subdivisions;
    let ring_index =
        super::super::boundary_first_ring_index::RingIndex::new(&rings_xy, 1.0, 1.0, (cell, cell));

    let mut uv_pts: Vec<(f64, f64)> = Vec::new();
    let mut pos_pts: Vec<Option<Point3>> = Vec::new();
    let mut ring_ranges = Vec::new();
    for ring in &rings_xy {
        let start = uv_pts.len();
        uv_pts.extend(ring.uv.iter().copied());
        pos_pts.extend(ring.pos.iter().copied().map(Some));
        ring_ranges.push((start, uv_pts.len()));
    }
    request_long_segment_splits(&rings_xy, cell, split_requests);
    // Stations offset half a cell from the box's corners (ring samples).
    let n = subdivisions.ceil() as usize + 1;
    for i in 0..n {
        for j in 0..n {
            let p = (
                lo.0 + (i as f64 + 0.5) * cell,
                lo.1 + (j as f64 + 0.5) * cell,
            );
            if ring_index.inside(p) && !ring_index.within_clearance(p) {
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
        why("triangulation");
        return None;
    };
    let deviation = super::super::boundary_first_deviation::Chart {
        face,
        surface: dev,
        planar: false,
        scale: (1.0, 1.0),
        tri_scale: (1.0, 0.0, 1.0),
        ring_index: &ring_index,
        ring_ranges: &ring_ranges,
        alias: &alias,
    };
    super::super::boundary_first_deviation::refine_by_deviation(
        &deviation,
        &mut uv_pts,
        &mut pos_pts,
        &mut triangles,
        options.chord_height,
    );
    if !covers_every_segment(&rings_xy, &ring_ranges, &alias, &triangles, split_requests) {
        why("a ring segment uncovered");
        return None;
    }
    // Vertices: ring samples keep their face parameters and positions;
    // interior points take theirs through the development.
    let exact = |p: &Point3| SINGULAR_EXACT * p.x.abs().max(p.y.abs()).max(p.z.abs()).max(1.0);
    let n_ring = ring_ranges.last().map_or(0, |r| r.1);
    let mut vertices = Vec::with_capacity(uv_pts.len());
    for (i, &xy) in uv_pts.iter().enumerate() {
        let known = (i < n_ring)
            .then(|| pos_pts[i])
            .flatten()
            .and_then(|p| uv_at.get(&bits(&p)).map(|&uv| (uv, p)));
        let (uv, position) = known.unwrap_or_else(|| {
            let uv = dev.face_uv(xy.0, xy.1);
            (
                uv,
                pos_pts[i].unwrap_or_else(|| surface.point_at(uv.0, uv.1)),
            )
        });
        let n = chart_normal(descriptor, surface, uv, position, exact(&position))
            .unwrap_or(Vec3::new(0.0, 0.0, 0.0));
        vertices.push(MeshVertex {
            position,
            normal: if orientation { n } else { -n },
            uv,
        });
    }
    // The development keeps the face chart's handedness (φ runs against
    // u), so the face chart's winding rule applies.
    if !orientation {
        for t in &mut triangles {
            t.swap(1, 2);
        }
    }
    if popkern_core::trace::trace_enabled("boundary-first") {
        eprintln!(
            "    [bf] face {face:?} served in its cone development: {} rings, {} triangles",
            rings_xy.len(),
            triangles.len()
        );
    }
    Some(FaceMesh {
        face,
        vertices,
        triangles,
    })
}

/// A ring with each run of consecutive samples at one position kept as
/// its first: the merged segment runs from that sample to the next
/// distinct one, its source spanning both edge spans when they are one
/// edge's (else the later span's: a split there lands beyond the run).
pub(super) fn merge_coincident(ring: BoundaryRing) -> BoundaryRing {
    let n = ring.uv.len();
    let mut out = BoundaryRing {
        uv: Vec::with_capacity(n),
        pos: Vec::with_capacity(n),
        sources: Vec::with_capacity(n),
    };
    for k in 0..n {
        if let Some(last) = out.pos.last()
            && last.distance_to(&ring.pos[k]) <= 0.0
        {
            if let Some(src) = out.sources.last_mut() {
                let next = ring.sources[k];
                *src = if src.0 == next.0 {
                    (src.0, src.1, next.2)
                } else {
                    next
                };
            }
            continue;
        }
        out.uv.push(ring.uv[k]);
        out.pos.push(ring.pos[k]);
        out.sources.push(ring.sources[k]);
    }
    // The run closing the ring onto its first sample.
    while out.pos.len() > 1
        && out
            .pos
            .last()
            .zip(out.pos.first())
            .is_some_and(|(a, b)| a.distance_to(b) <= 0.0)
    {
        out.uv.pop();
        out.pos.pop();
        out.sources.pop();
    }
    out
}

/// Ring segments longer than one and a half cells request their
/// splits, all in one round (in space: the development is an isometry).
fn request_long_segment_splits(
    rings_xy: &[BoundaryRing],
    cell: f64,
    split_requests: &mut Vec<SplitRequest>,
) {
    let dedup = crate::encroachment::SAMPLE_DEDUP;
    for ring in rings_xy {
        let m = ring.uv.len();
        for k in 0..m {
            let (edge, f0, f1) = ring.sources[k];
            let length = ring.pos[k].distance_to(&ring.pos[(k + 1) % m]);
            if length <= 1.5 * cell || edge == EdgeHandle::default() {
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
}

/// Whether every ring segment is a triangle edge; an uncovered one
/// requests its span's midpoint.
fn covers_every_segment(
    rings_xy: &[BoundaryRing],
    ring_ranges: &[(usize, usize)],
    alias: &[usize],
    triangles: &[[usize; 3]],
    split_requests: &mut Vec<SplitRequest>,
) -> bool {
    let mut edges = popkern_math::hash::FastHashSet::default();
    for t in triangles {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            edges.insert((x.min(y), x.max(y)));
        }
    }
    let sources: Vec<_> = rings_xy
        .iter()
        .flat_map(|r| r.sources.iter().copied())
        .collect();
    let mut covered = true;
    for &(start, end) in ring_ranges {
        let m = end - start;
        for k in 0..m {
            let (a, b) = (start + k, start + (k + 1) % m);
            let (ca, cb) = (alias[a], alias[b]);
            if ca != cb && !edges.contains(&(ca.min(cb), ca.max(cb))) {
                covered = false;
                let (edge, f0, f1) = sources[a];
                split_requests.push((edge, 0.5 * (f0 + f1)));
            }
        }
    }
    covered
}
