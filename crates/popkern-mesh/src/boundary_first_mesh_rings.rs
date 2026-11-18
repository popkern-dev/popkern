//! The boundary-first pipeline from charted rings to a face mesh: the
//! interior lattice, the constrained triangulation, deviation refinement,
//! coverage and materialisation. A face's rings reach here charted (lifted
//! across cuts, wraps joined); an atlas band piece takes the same path.

use popkern_core::geometry::chart::ChartDescriptor;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{EdgeHandle, FaceHandle};

use super::*;
use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};

/// Mesh a face from its charted `rings` (ring 0 the outer boundary).
pub(super) fn mesh_rings(
    (face_handle, surface, orientation): (FaceHandle, &SurfaceGeometry, bool),
    mut rings: Vec<BoundaryRing>,
    options: &TessellationOptions,
    split_requests: &mut Vec<SplitRequest>,
    descriptor: &ChartDescriptor,
    sizing: Option<((f64, f64), (f64, f64))>,
) -> Option<FaceMesh> {
    if rings.is_empty() {
        return None;
    }
    let chord = options.chord_height;
    boundary_first_unwrap::place_rings_on_chart(&mut rings, surface, orientation, chord);
    if popkern_core::trace::trace_enabled("boundary-first-ring") {
        for (ri, ring) in rings.iter().enumerate() {
            eprintln!(
                "    [bf-ring] face {face_handle:?} ring {ri} ({} pts)",
                ring.uv.len()
            );
            for (i, ((u, v), p)) in ring.uv.iter().zip(ring.pos.iter()).enumerate() {
                eprintln!(
                    "      [bf-ring] {i} uv=({u:.4},{v:.4}) pos=({:.2},{:.2},{:.2})",
                    p.x, p.y, p.z
                );
            }
        }
    }
    // Four-sided ruled charts (fillet/chamfer strips) take the
    // structured tensor path; over-chord cells request boundary splits.
    if let Some(fm) = boundary_first_ruled_strip::try_tessellate(
        face_handle,
        surface,
        orientation,
        &rings,
        options,
    ) {
        let refine = boundary_first_strip_refine::request_deviation_splits;
        refine(
            &fm,
            (surface, descriptor),
            &rings[0],
            chord,
            options.angular_bound(),
            split_requests,
        );
        return Some(fm);
    }

    // UV conditioning for the incircle test (`boundary_first_ring_condition`).
    let exact_uv = boundary_first_ring_condition::condition_for_incircle(&mut rings);

    // Point set: boundary points (fixed) + interior grid points.
    let mut uv_pts: Vec<(f64, f64)> = Vec::new();
    let mut pos_pts: Vec<Option<Point3>> = Vec::new(); // None → eval surface
    let mut ring_ranges: Vec<(usize, usize)> = Vec::new();
    let mut sources: Vec<(EdgeHandle, f64, f64)> = Vec::new();
    for ring in &rings {
        let start = uv_pts.len();
        uv_pts.extend(ring.uv.iter().copied());
        pos_pts.extend(ring.pos.iter().copied().map(Some));
        sources.extend(ring.sources.iter().copied());
        ring_ranges.push((start, uv_pts.len()));
    }

    // Interior grid from the outer ring's bbox; keep points strictly
    // inside the trim region and away from the boundary (a third of
    // a cell) so Delaunay doesn't create slivers against the rings.
    let (mut u0, mut u1, mut v0, mut v1) = (
        f64::INFINITY,
        f64::NEG_INFINITY,
        f64::INFINITY,
        f64::NEG_INFINITY,
    );
    for &(u, v) in &rings[0].uv {
        u0 = u0.min(u);
        u1 = u1.max(u);
        v0 = v0.min(v);
        v1 = v1.max(v);
    }
    // Station counts (`boundary_first_lattice::face_stations`): a piece of
    // a face (an atlas band) sized over the face's box, laid over its own.
    let (n, n_u, n_v) = boundary_first_lattice::face_stations(
        (surface, descriptor),
        ((u0, u1), (v0, v1)),
        sizing,
        options,
    );
    // Chart metric scale (mm per parameter unit at the chart
    // centre). The chart mixes units (radians × millimetres on a
    // cylinder), and every proximity quantity below — cell size,
    // clearance band, protecting disks, segment refinement — must be
    // measured in the surface metric, not raw UV: the raw clearance
    // on a long thin strip (v in mm, u in radians) exceeded the
    // whole chart's u half-span and evacuated the ENTIRE interior
    // lattice (measured: plate-rim strips meshed boundary-only,
    // every triangle a rail-to-rail chord, volume 0.5% low — PK-218).
    // The same scale feeds the Delaunay incircle predicate below.
    let (su, sv, shear, sv_sheared) = {
        let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let p = surface.point_at(um, vm);
        let eps_u = ((u1 - u0).abs() * 1e-4).max(1e-9);
        let eps_v = ((v1 - v0).abs() * 1e-4).max(1e-9);
        let pu = (surface.point_at(um + eps_u, vm) - p) * (1.0 / eps_u);
        let pv = (surface.point_at(um, vm + eps_v) - p) * (1.0 / eps_v);
        let (du, dv) = (pu.length().max(1e-9), pv.length().max(1e-9));
        // The chart's shear (F of the first fundamental form) enters
        // the Delaunay metric through its Cholesky factor: a twisted
        // blade's chart maps uv-good triangles to near-collinear 3D
        // caps under a diagonal scale (E4: propeller facets 176-178 deg
        // obtuse, leaning 30-40 deg off their vertex normals).
        let f = pu.dot(&pv);
        let shear = f / du;
        let rest = (dv * dv - shear * shear).max((1e-3 * dv).powi(2)).sqrt();
        (du, dv, shear, rest)
    };
    // Per-direction lattice spacings (metric).
    let cells = ((su * (u1 - u0)) / n_u as f64, (sv * (v1 - v0)) / n_v as f64);
    // Lattice cell (diagnostics): the finer per-direction spacing.
    let cell = ((su * (u1 - u0)) / n_u as f64)
        .min((sv * (v1 - v0)) / n_v as f64)
        .max(1e-12);
    let clearance = cell / 3.0;
    // Long boundary segments split toward the lattice spacing
    // (`boundary_first_segment_split`). Planar charts skip it.
    if !matches!(
        surface,
        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
    ) {
        let lattice = boundary_first_segment_split::SegmentLattice {
            scale: (su, sv),
            cells,
            square_cell: ((su * (u1 - u0)) / n as f64)
                .max((sv * (v1 - v0)) / n as f64)
                .max(1e-12),
            directional: n_u != n || n_v != n,
            turn: options
                .angular_bound()
                .filter(|_| {
                    matches!(
                        surface,
                        popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)
                    )
                })
                .map(|bound| (surface, bound)),
        };
        boundary_first_segment_split::request_long_segment_splits(&rings, &lattice, split_requests);
    }
    // Lattice points are filtered by region and a per-segment
    // clearance band (a third of the lattice spacing across the
    // segment). Points inside a segment's diametral disk are KEPT: the
    // triangulation recovers every ring segment by constrained flips
    // (`boundary_first_cdt`), where the unconstrained Delaunay used to
    // need disk-free segments and so evacuated the interior beside
    // every long one (PK-651). The ring index answers by bucket,
    // identically to scanning every segment (PK-589).
    let ring_index = boundary_first_ring_index::RingIndex::new(&rings, su, sv, cells);
    // RULED charts (v-midline deviation below chord height: cylinder
    // and cone walls, extruded NURBS) skip the interior grid: full-
    // height boundary strips are EXACT along the ruling, and the
    // adaptive edge sampling makes them as dense in u as the rims —
    // relative volume deficit ~1e-4. Interior grid points are capped
    // at max_subdivisions (u-spacing 2π/64), so any that survive the
    // protection filters REPLACE those strips with triangles carrying
    // 16× the chord error (measured: the chamfered column's volume
    // gate at 4.6e-3 vs its 2e-3 budget once the boundary splits
    // shrank the seam's protecting disks enough for grid points to
    // survive). Doubly-curved charts keep the grid — strips would
    // leave their interiors flat.
    let ruled_v = !matches!(
        surface,
        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
    ) && boundary_first_ruled_lattice::v_line_straight(
        surface,
        0.5 * (u0 + u1),
        (v0, v1),
        &boundary_first_ruled_lattice::v_line_fractions(surface, v0, v1),
        options.chord_height,
    );
    // Delaunay metric: the surface metric, v compressed on ruled lattices.
    let mut tri_scale = (su, shear, sv_sheared);
    if ruled_v {
        let compress = boundary_first_ruled_lattice::ruled_lattice(
            &rings,
            &ring_index,
            (u0, u1, v0, v1),
            (su, sv),
            cells,
            options,
            &mut uv_pts,
            &mut pos_pts,
        );
        if boundary_first_ruled_lattice::ruled_everywhere(surface, (u0, u1, v0, v1), chord) {
            tri_scale.1 *= compress;
            tri_scale.2 *= compress;
        }
    } else {
        let mut census = [0usize; 3]; // region, clearance, kept
        for i in 1..n_u {
            for j in 1..n_v {
                let u = u0 + (u1 - u0) * i as f64 / n_u as f64;
                let v = v0 + (v1 - v0) * j as f64 / n_v as f64;
                if !ring_index.inside((u, v)) {
                    census[0] += 1;
                } else if ring_index.within_clearance((u, v)) {
                    census[1] += 1;
                } else {
                    census[2] += 1;
                    uv_pts.push((u, v));
                    pos_pts.push(None);
                }
            }
        }
        // Lattice kill census + boundary-segment metric stats: the
        // diagnosis surface for evacuated interiors (a depth-first
        // sampling cap once left ONE seam segment spanning half the
        // chart whose diametral disk killed 50% of the lattice —
        // "max" names that class instantly).
        if popkern_core::trace::trace_enabled("boundary-first-ring") {
            let mut seg_lens: Vec<f64> = Vec::new();
            let mut worst = ((0.0, 0.0), (0.0, 0.0), 0.0_f64);
            for ring in &rings {
                let m = ring.uv.len();
                for k in 0..m {
                    let a = ring.uv[k];
                    let b = ring.uv[(k + 1) % m];
                    let l = ((su * (b.0 - a.0)).powi(2) + (sv * (b.1 - a.1)).powi(2)).sqrt();
                    if l > worst.2 {
                        worst = (a, b, l);
                    }
                    seg_lens.push(l);
                }
            }
            seg_lens.sort_by(f64::total_cmp);
            let med = seg_lens.get(seg_lens.len() / 2).copied().unwrap_or(0.0);
            let long = seg_lens.iter().filter(|&&l| l > 3.0 * cell).count();
            eprintln!(
                "    [bf-census] face {face_handle:?} n_u={n_u} n_v={n_v}: region={} clearance={} kept={} | segs {} med {med:.4} max {:.4} at ({:.4},{:.4})->({:.4},{:.4}) long(>3cell) {long}",
                census[0],
                census[1],
                census[2],
                seg_lens.len(),
                worst.2,
                worst.0.0,
                worst.0.1,
                worst.1.0,
                worst.1.1
            );
        }
    }

    if popkern_core::trace::trace_enabled("boundary-first-ring") {
        let boundary_n: usize = ring_ranges.iter().map(|&(s, e)| e - s).sum();
        eprintln!(
            "    [bf-grid] face {face_handle:?}: boundary {boundary_n} interior {} (n={n} cell={cell:.4} clearance={clearance:.4})",
            uv_pts.len() - boundary_n
        );
    }

    // Metric-scaled Delaunay filtered to the region.
    let (mut triangles, alias) = boundary_first_triangulate::triangulate_region(
        face_handle,
        &uv_pts,
        tri_scale,
        &ring_index,
        &ring_ranges,
    )?;
    let chart = boundary_first_deviation::Chart {
        face: face_handle,
        surface,
        planar: matches!(surface, SurfaceGeometry::Plane { .. }),
        scale: (su, sv),
        tri_scale,
        ring_index: &ring_index,
        ring_ranges: &ring_ranges,
    };
    boundary_first_deviation::refine_by_deviation(
        &chart,
        &mut uv_pts,
        &mut pos_pts,
        &mut triangles,
        chord,
    );
    if popkern_core::trace::trace_enabled("boundary-first-ring") {
        let mut areas: Vec<f64> = triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
                ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs()
            })
            .collect();
        let mut sorted = areas.clone();
        sorted.sort_by(f64::total_cmp);
        let median = sorted.get(sorted.len() / 2).copied().unwrap_or(0.0);
        for (t, area) in triangles.iter().zip(areas.iter_mut()) {
            if *area > median * 50.0 && median > 0.0 {
                let (a, b, c) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
                eprintln!(
                    "    [bf-outlier] face {face_handle:?} area {area:.3} (median {median:.4}) uv ({:.3},{:.3})({:.3},{:.3})({:.3},{:.3})",
                    a.0, a.1, b.0, b.1, c.0, c.1
                );
            }
        }
    }
    // Boundary-edge coverage check: every ring segment must appear as
    // a triangle edge, else the Delaunay cut a concave corner — fall
    // back to the uniform-grid path rather than ship a gap.
    let mut edge_set = popkern_math::hash::FastHashSet::default();
    for t in &triangles {
        for (x, y) in [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])] {
            edge_set.insert((x.min(y), x.max(y)));
        }
    }
    let mut uncovered = false;
    for &(start, end) in &ring_ranges {
        let m = end - start;
        for k in 0..m {
            let a = start + k;
            let b = start + (k + 1) % m;
            let (ca, cb) = (alias[a], alias[b]);
            if ca != cb && !edge_set.contains(&(ca.min(cb), ca.max(cb))) {
                uncovered = true;
                // Nearest OTHER point to the segment: the encroacher
                // (or a T-point) that displaced it from the Delaunay.
                let (pa, pb) = (uv_pts[a], uv_pts[b]);
                let seg_proj = |p: (f64, f64)| -> (f64, f64) {
                    let (dx, dy) = (pb.0 - pa.0, pb.1 - pa.1);
                    let l2 = (dx * dx + dy * dy).max(1e-300);
                    let t = (((p.0 - pa.0) * dx + (p.1 - pa.1) * dy) / l2).clamp(0.0, 1.0);
                    let d = ((p.0 - pa.0 - dx * t).powi(2) + (p.1 - pa.1 - dy * t).powi(2)).sqrt();
                    (d, t)
                };
                let nearest = uv_pts
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| *i != a && *i != b)
                    .map(|(i, &p)| {
                        let (d, t) = seg_proj(p);
                        (i, d, t)
                    })
                    .min_by(|x, y| x.1.partial_cmp(&y.1).unwrap_or(std::cmp::Ordering::Equal));
                if popkern_core::trace::trace_enabled("boundary-first") {
                    eprintln!("    [bf] span: a={:?}", sources.get(a));
                    eprintln!(
                        "    [bf] ring segment {a}-{b} uncovered at uv ({:.4},{:.4})-({:.4},{:.4}) nearest pt {:?} d={:.4} uv=({:.4},{:.4})",
                        pa.0,
                        pa.1,
                        pb.0,
                        pb.1,
                        nearest.map(|(i, _, _)| i),
                        nearest.map(|(_, d, _)| d).unwrap_or(f64::NAN),
                        nearest.map(|(i, _, _)| uv_pts[i].0).unwrap_or(f64::NAN),
                        nearest.map(|(i, _, _)| uv_pts[i].1).unwrap_or(f64::NAN),
                    );
                }
                // Refinement. A T-POINT (a foreign ring point lying
                // ON the segment) can never be resolved by midpoint
                // bisection — the point stays interior to some
                // sub-segment forever (measured livelock: an arm
                // socket's ring point at v=12 on a full-height
                // straight run). CONFORM instead: split the
                // canonical span at the T-point's own projection
                // parameter, so both incident faces gain a sample AT
                // the point and the sub-segments become coverable.
                // Genuine encroachers (off-segment) keep the
                // Ruppert-style midpoint split.
                let (ea, f_lo, f_hi) = sources[a];
                let seg_len = ((pb.0 - pa.0).powi(2) + (pb.1 - pa.1).powi(2)).sqrt();
                let conforming = nearest
                    .filter(|(_, d, t)| {
                        *d < 1e-6 * seg_len.max(1.0) && *t > 1e-6 && *t < 1.0 - 1e-6
                    })
                    .map(|(_, _, t)| t);
                let f = match conforming {
                    Some(t) => f_lo + (f_hi - f_lo) * t,
                    None => 0.5 * (f_lo + f_hi),
                };
                split_requests.push((ea, f));
            }
        }
    }
    if uncovered {
        if popkern_core::trace::trace_enabled("boundary-first-dump") {
            eprintln!("    [bf-dump] face {face_handle:?} rings {ring_ranges:?}");
            for (i, &(u, v)) in uv_pts.iter().enumerate() {
                eprintln!("    [bf-dump] pt {i} {u:.17e} {v:.17e}");
            }
        }
        return None;
    }

    // Facing here follows the stored `orientation` flag; per-face
    // witnesses are NOT universally honest (extrude caps lie in the
    // flag, chamfer cones in the loop winding), so the final
    // authority is the global winding-coherence pass in `tessellate`
    // (`orient_face_meshes_coherently`), which propagates manifold
    // consistency across shared canonical boundary vertices.
    // Materialize vertices.
    boundary_first_ring_condition::restore_exact_uv(&mut uv_pts, &exact_uv);
    let mut vertices = Vec::with_capacity(uv_pts.len());
    for (i, &(u, v)) in uv_pts.iter().enumerate() {
        let position = pos_pts[i].unwrap_or_else(|| surface.point_at(u, v));
        let mut normal = surface.normal_at(u, v);
        let len = normal.length();
        if len > 1e-300 {
            normal = Vec3::new(normal.x / len, normal.y / len, normal.z / len);
        }
        if !orientation {
            normal = -normal;
        }
        vertices.push(MeshVertex {
            position,
            normal,
            uv: (u, v),
        });
    }
    if !orientation {
        for t in &mut triangles {
            t.swap(1, 2);
        }
    }
    boundary_first_fold::drop_degenerate_slivers(&mut triangles, &vertices);
    boundary_first_fidelity::trace_bf_volume(face_handle, &vertices, &triangles, &rings);
    if !closed_chart_fidelity_ok(
        face_handle,
        surface,
        descriptor,
        &vertices,
        &triangles,
        options,
    ) {
        return None;
    }
    Some(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}
