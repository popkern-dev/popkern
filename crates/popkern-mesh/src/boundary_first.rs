//! Boundary-first face tessellation (PK-095 closure).
//!
//! Face meshes sampled independently leave T-junction cracks along
//! every shared edge (the cap-rim "flap ring"). This module makes
//! adjacent faces watertight BY CONSTRUCTION: every edge is sampled
//! ONCE into a canonical polyline (adaptive to the 3D curve's
//! flatness), and each face's boundary vertices ARE those canonical
//! points — both sides of an edge use bit-identical positions. The
//! face interior is filled with grid points and triangulated against
//! the fixed boundary by a Bowyer–Watson Delaunay in UV, filtered to
//! the trim polygon.
//!
//! Faces without complete pcurve loops (or without edge curves) fall
//! back to the uniform grid+clip path in `tessellate_face`.
//!
//! Observability: `POPKERN_TRACE=boundary-first` logs, per face falling
//! back to the uniform grid, WHY (missing pcurve, no canonical samples,
//! off-surface residual, uncovered ring segment, unchartable wrap).

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};

#[path = "boundary_first_cdt.rs"]
mod boundary_first_cdt;
#[path = "boundary_first_delaunay.rs"]
mod boundary_first_delaunay;
#[path = "boundary_first_deviation.rs"]
mod boundary_first_deviation;
#[path = "boundary_first_fidelity.rs"]
mod boundary_first_fidelity;
#[path = "boundary_first_fold.rs"]
mod boundary_first_fold;
use boundary_first_fidelity::closed_chart_fidelity_ok;
#[path = "boundary_first_edge_samples.rs"]
mod boundary_first_edge_samples;
#[path = "boundary_first_lattice.rs"]
mod boundary_first_lattice;
#[path = "boundary_first_pole_cap.rs"]
mod boundary_first_pole_cap;
#[path = "boundary_first_revolution.rs"]
mod boundary_first_revolution;
#[path = "boundary_first_revolution_bipolar.rs"]
mod boundary_first_revolution_bipolar;
#[path = "boundary_first_ring_condition.rs"]
mod boundary_first_ring_condition;
#[path = "boundary_first_ring_index.rs"]
mod boundary_first_ring_index;
#[path = "boundary_first_ruled_lattice.rs"]
mod boundary_first_ruled_lattice;
#[path = "boundary_first_ruled_strip.rs"]
mod boundary_first_ruled_strip;
#[path = "boundary_first_segment_split.rs"]
mod boundary_first_segment_split;
#[path = "boundary_first_slit.rs"]
mod boundary_first_slit;
#[path = "boundary_first_strip_refine.rs"]
mod boundary_first_strip_refine;
#[path = "boundary_first_triangulate.rs"]
mod boundary_first_triangulate;
#[path = "boundary_first_unwrap.rs"]
mod boundary_first_unwrap;
#[path = "boundary_first_wrap_pair.rs"]
mod boundary_first_wrap_pair;
#[path = "boundary_first_wrap_slit.rs"]
mod boundary_first_wrap_slit;
pub(crate) use boundary_first_edge_samples::{EdgeSamples, sample_edges};
use boundary_first_wrap_pair::{join_wrapped_rings, stash_wrapped_ring};

/// One boundary ring: ordered UV polygon with matched 3D positions
/// and, per point, the canonical edge + curve fraction it came from
/// (the handle encroachment refinement needs to split the RIGHT
/// shared edge span).
struct BoundaryRing {
    uv: Vec<(f64, f64)>,
    pos: Vec<Point3>,
    /// Per point: the canonical edge span the ring segment STARTING
    /// at this point traverses — (edge, fraction at this point,
    /// fraction at the next walk sample). Encroachment refinement
    /// splits at the span midpoint; recording both ends makes the
    /// request correct regardless of the walk direction over the
    /// edge (reversed walks descend in fraction).
    sources: Vec<(EdgeHandle, f64, f64)>,
}

/// A request to refine the canonical samples of one edge at one
/// curve fraction, issued when Delaunay coverage fails because a
/// boundary segment is encroached (Ruppert-style split-and-retry;
/// both incident faces re-tessellate with the shared refined
/// samples, so watertightness is preserved by construction).
pub(crate) type SplitRequest = (EdgeHandle, f64);

/// Tessellate a face against canonical edge samples. Returns `None`
/// when the face lacks the data (pcurves / curves) — caller falls
/// back to the uniform-grid path.
pub(crate) fn tessellate_face_boundary_first(
    brep: &BRep,
    face_handle: FaceHandle,
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    options: &TessellationOptions,
    split_requests: &mut Vec<SplitRequest>,
) -> Option<FaceMesh> {
    let face = brep.faces.get(face_handle)?;
    let surface_on_face = face.surface.as_ref()?;
    let surface = &surface_on_face.geometry;
    let orientation = surface_on_face.orientation;

    // A face whose outer loop traverses EVERY edge twice is a pure
    // seam spine (a full sphere's meridian pair): its chart "polygon"
    // is two overlapping border lines that bound no interior, and
    // Delaunay on it emits only degenerate slivers (measured: the
    // primitive sphere meshed 0/46 non-degenerate once its seam
    // carried real arc geometry — the chord-era ring build failed and
    // silently routed these faces to the pcurve-grid path,
    // which handles the periodic rectangle directly). Keep that
    // routing explicit.
    {
        let mut edge_counts: HashMap<EdgeHandle, usize> = HashMap::new();
        let start = brep.loops.get(face.outer_loop)?.half_edge;
        let mut cur = start;
        let mut steps = 0;
        loop {
            let he = brep.half_edges.get(cur)?;
            *edge_counts.entry(he.edge).or_insert(0) += 1;
            cur = he.next;
            steps += 1;
            if cur == start || steps > brep.half_edges.len() {
                break;
            }
        }
        if !edge_counts.is_empty() && edge_counts.values().all(|&c| c >= 2) {
            // PK-339: the two-apex revolution arm serves the classic
            // member of this class (a full sphere: two poles joined by
            // one meridian seam) by construction — attempt it before
            // falling back. Anything it declines keeps the pcurve-grid
            // route (Delaunay stays skipped for seam spines).
            return boundary_first_revolution_bipolar::seam_spine_revolution(
                brep,
                face_handle,
                edge_samples,
                options,
            );
        }
    }

    // Build boundary rings from canonical samples.
    let mut rings: Vec<BoundaryRing> = Vec::new();
    let (mut wrapped_rings, mut wrapped_period) = (Vec::<BoundaryRing>::new(), None::<(f64, bool)>); // PK-468
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());
    for lh in loops {
        let built = build_ring(brep, lh, edge_samples)?;
        // Slit dissolution (PK-335): coincident double-traversals
        // (shell rim slit arms) are chart artifacts — cut the ring
        // into its genuine boundary components before charting, or
        // Delaunay coverage chases the zero-width arms forever.
        for mut ring in boundary_first_slit::dissolve_slits(built, surface.parameter_periodicity())
        {
            // Periodic surfaces (analytic cylinders: u = angle) can have
            // pcurves authored in different period branches; a ring that
            // crosses the seam jumps ~2π in the UV chart and Delaunay
            // cannot cover the jump segment. Unwrap u to the branch
            // nearest each previous sample (positions are canonical 3D
            // points and unaffected; cos/sin accept any branch).
            {
                let (period_u, period_v) = surface.parameter_periodicity();
                // Chart-degenerate samples (poles, apexes) and sphere
                // pole crossings (`chart_degenerate`, `resolve_pole_crossings`).
                let mut u_degenerate =
                    boundary_first_unwrap::chart_degenerate(&ring, surface, period_u);
                boundary_first_pole_cap::resolve_pole_crossings(
                    &mut ring,
                    &mut u_degenerate,
                    surface,
                    orientation,
                );
                if let Some(p) = period_u
                    && !boundary_first_unwrap::unwrap_dir(
                        &mut ring,
                        p,
                        |uv| &mut uv.0,
                        &u_degenerate,
                    )
                {
                    // A genuine full wrap through an apex point is the
                    // structured-revolution class (PK-336); sole-loop
                    // faces only — holes keep the pcurve-grid route.
                    if face.inner_loops.is_empty()
                        && rings.is_empty()
                        && let Some(fm) = boundary_first_revolution::try_tessellate(
                            face_handle,
                            surface,
                            orientation,
                            &ring,
                            &u_degenerate,
                            p,
                            options,
                        )
                    {
                        return Some(fm);
                    }
                    stash_wrapped_ring(ring, &mut wrapped_rings); // PK-468 band slit
                    boundary_first_wrap_pair::record_wrap(&mut wrapped_period, (p, false))?;
                    continue;
                }
                if let Some(p) = period_v
                    && !boundary_first_unwrap::unwrap_dir(&mut ring, p, |uv| &mut uv.1, &[])
                {
                    boundary_first_wrap_pair::record_wrap(&mut wrapped_period, (p, true))?;
                    boundary_first_wrap_pair::stash_transposed(ring, &mut wrapped_rings);
                    continue; // a v-wrapped band (an imported pipe elbow) slits transposed
                }
            }
            for lobe in boundary_first_slit::dissolve_pinches(ring) {
                if lobe.uv.len() >= 3 {
                    rings.push(lobe);
                }
            }
        }
    }
    if let Some((p, transposed)) = wrapped_period {
        // The joined band polygon is the face's OUTER boundary: ring 0
        // by convention (hole-sheet alignment and the lattice's bbox
        // read it). Appended after the holes, a cross-drilled wall's
        // holes stayed a period off the band and their region meshed
        // as material (bracket1: facets 2.11 off; PK-652).
        let joined = join_wrapped_rings(
            face_handle,
            &wrapped_rings,
            p,
            surface,
            orientation,
            options.chord_height,
            transposed,
        )?;
        rings.splice(0..0, joined);
    }
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
            surface,
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
    // Refine like the uniform-grid path would for curved surfaces: sample
    // the surface diagonal and bump density if it deviates from the
    // chord more than chord_height (then the angular facet bound).
    let mut n = options.min_subdivisions.max(1);
    {
        let p00 = surface.point_at(u0, v0);
        let p11 = surface.point_at(u1, v1);
        let pm = surface.point_at(0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let chord_mid = Point3::new(
            0.5 * (p00.x + p11.x),
            0.5 * (p00.y + p11.y),
            0.5 * (p00.z + p11.z),
        );
        let dev = pm.distance_to(&chord_mid);
        while (dev / (n as f64 * n as f64)) > options.chord_height && n < options.max_subdivisions {
            n *= 2;
        }
        n = crate::facet_angle::densify_lattice(surface, ((u0, u1), (v0, v1)), n, options);
    }
    // Per-direction station counts (see `boundary_first_lattice::lattice_counts`).
    let (n_u, n_v) =
        boundary_first_lattice::lattice_counts(surface, (u0, u1), (v0, v1), n, options);
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
    let mut edge_set = std::collections::HashSet::new();
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
    if !closed_chart_fidelity_ok(face_handle, surface, &vertices, &triangles, options) {
        return None;
    }
    Some(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
}

/// The surface-on-face of a half-edge's loop (ring-walk helper).
fn sof_ring_geometry<'a>(
    brep: &'a BRep,
    data: &popkern_core::topology::HalfEdge,
) -> Option<&'a popkern_core::geometry::surface::SurfaceOnFace> {
    brep.faces[brep.loops[data.loop_].face].surface.as_ref()
}

/// Walk one loop, concatenating canonical edge samples in loop order.
/// UV comes from each half-edge's pcurve at the matching fraction
/// (reversed when the half-edge runs against the edge curve); 3D is
/// the shared canonical point.
fn build_ring(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
) -> Option<BoundaryRing> {
    let start = brep.loops.get(loop_handle)?.half_edge;
    let start = boundary_first_edge_samples::rotate_ring_start(brep, start);
    // Continuity-based closed-edge entry applies only to MIXED loops
    // (a closed edge following an open predecessor — the PK-454
    // seam-split defect class). All-closed loops (a ring tube's
    // seam + rim pairs) carry pcurves authored as one consistent
    // family and keep the historical direction rules — re-deciding
    // them by raw uv distance broke a served pipe ring (127 open
    // mesh edges).
    let loop_has_open = boundary_first_edge_samples::edge_endpoint_positions(brep, start)
        .is_some_and(|(a, b)| a.distance_to(&b) >= 1.0e-9);
    // Chain-continuity period lift + folded-seam unfold (PK-457, same
    // scope as the PK-454 entry rule; see
    // `boundary_first_fold::fold_shift_for_leg`).
    let ring_periods = boundary_first_fold::ring_chain_periods(brep, loop_handle, loop_has_open);
    let mut uv = Vec::new();
    let mut pos = Vec::new();
    let mut sources = Vec::new();
    let mut first_leg_starts: HashMap<EdgeHandle, (f64, f64)> = HashMap::new();
    let mut he = start;
    for _ in 0..10_000 {
        let data = &brep.half_edges[he];
        let eh = data.edge;
        // ZERO-LENGTH edges (degenerate boolean residue at discarded
        // tool corners) contribute nothing to the ring: skip them
        // entirely rather than bailing the whole face to the uniform-grid
        // path. Endpoints resolve through twin-or-next so boundary
        // sides with a dead twin do not panic.
        let ends = boundary_first_edge_samples::edge_endpoint_positions(brep, he);
        if ends.is_some_and(|(va, vb)| va.distance_to(&vb) < 1.0e-9)
            && !edge_samples.contains_key(&eh)
        {
            he = data.next;
            if he == start {
                break;
            }
            continue;
        }
        let Some(pc) = data.pcurve.as_ref() else {
            if popkern_core::trace::trace_enabled("boundary-first") {
                eprintln!("    [bf] edge {eh:?}: missing pcurve");
            }
            return None;
        };
        let Some(samples) = edge_samples.get(&eh) else {
            if popkern_core::trace::trace_enabled("boundary-first") {
                eprintln!("    [bf] edge {eh:?}: no canonical samples");
            }
            return None;
        };
        let origin = brep.vertices[data.origin].position;

        let (closed_edge, closed_entry_forward, forward) =
            boundary_first_edge_samples::ring_walk_direction(
                pc,
                &sof_ring_geometry(brep, data)?.geometry,
                samples,
                if loop_has_open {
                    uv.last().copied()
                } else {
                    None
                },
                &origin,
            )?;
        let m = samples.points.len();
        let walk = (closed_edge, forward, closed_entry_forward);
        boundary_first_fold::trace_ring_leg(he, eh, m, walk, pc, uv.last().copied());
        // Emit all but the last point (the next half-edge provides
        // its own start, which is the shared vertex).
        let (mut leg_shift, mut walk_forward) = (None::<(f64, f64)>, true);
        for k in 0..(m - 1) {
            let idx = if forward { k } else { m - 1 - k };
            let f = samples.fractions[idx];
            let sof = brep.faces[brep.loops[data.loop_].face].surface.as_ref()?;
            let pcurve_forward_along_walk = match closed_entry_forward {
                Some(entry_forward) => entry_forward,
                None => boundary_first_edge_samples::pcurve_forward_along_walk(
                    pc,
                    &sof.geometry,
                    &origin,
                    closed_edge,
                ),
            };
            // Seed at the sample's CURVE-WINDOW fraction, not the
            // uniform walk fraction k/(m-1): canonical samples are
            // ADAPTIVE (encroachment splits refine specific spans),
            // so index and curve fraction diverge — on the coil seam
            // a uniform-index seed sat half a pitch from its target
            // and the Newton refinement converged onto the ADJACENT
            // TURN of the helix (residual 0.156 ≈ the pitch; the
            // whole wall fell back to the uniform-grid path at fine
            // tolerances). The pcurve direction relative to the
            // CURVE is the walk-relative direction composed with the
            // sample orientation.
            let curve_f = samples.fractions[idx];
            let pcurve_along_curve = pcurve_forward_along_walk == forward;
            let tp = if pcurve_along_curve {
                pc.t_start + (pc.t_end - pc.t_start) * curve_f
            } else {
                pc.t_end - (pc.t_end - pc.t_start) * curve_f
            };
            // The pcurve is only a SEED: refine each UV by point
            // inversion of the canonical 3D sample so imprecise
            // pcurves (straight chords for bulging curves) still
            // yield exact boundary UVs. If the refined image stays
            // far from the canonical point, the curve genuinely does
            // not lie on this surface (sphere seam pole chord) —
            // fall back to the uniform-grid path.
            let uvp = pc.point_at(tp);
            let target = samples.points[idx];
            let (uu, vv, residual, gate) = boundary_first_edge_samples::refine_boundary_uv(
                &sof.geometry,
                (uvp.x, uvp.y),
                &target,
                samples,
            )?;
            // NOTE: widening this gate by recorded EDGE tolerance (which
            // blends record) was tried twice and is FALSIFIED both times:
            // it repartitions the boundary-first/uniform-grid face split
            // and the mixed boundaries crack far worse (hexagon6_blunt
            // 0 -> 482; pentagon5_blunt/G2 106 -> 492). Do not reintroduce
            // it. A FACE tolerance is import provenance only: a surface
            // missing its boundary at the file's precision (PK-644).
            let face_tolerance = brep.faces[brep.loops[data.loop_].face].tolerance;
            if residual > gate.max(face_tolerance.unwrap_or(0.0)) {
                if popkern_core::trace::trace_enabled("boundary-first") {
                    eprintln!(
                        "    [bf] edge {eh:?}: off-surface residual {residual:.3e} at sample {idx}/{m} \
                         target ({:.4},{:.4},{:.4}) seed uv ({:.4},{:.4}) refined ({uu:.4},{vv:.4})",
                        target.x, target.y, target.z, uvp.x, uvp.y
                    );
                }
                return None;
            }
            // Fold placement is decided at the leg's FIRST sample and
            // applied to the whole leg (a folded seam pair's second
            // leg belongs one period over; the chain lift would fight
            // the deliberate jump, so a shifted leg skips it). The
            // jump crosses a degenerate (pole) vertex, and the ring
            // gets BOTH its chart images — the crossing then runs
            // along the zero-3D-length pole row (a convex-hull edge
            // Delaunay always keeps) instead of a full-period
            // diagonal no refinement can ever cover (measured: one
            // seam span split per round toward the pole, all eight
            // rounds, on every revolved-NURBS dimple cap).
            let (uu, vv) = match leg_shift {
                Some(shift) if shift != (0.0, 0.0) => (uu + shift.0, vv + shift.1),
                Some(_) => {
                    boundary_first_fold::lift_to_chain((uu, vv), uv.last().copied(), ring_periods)
                }
                None => {
                    let shift = boundary_first_fold::fold_shift_for_leg(
                        &sof.geometry,
                        first_leg_starts.get(&eh).copied(),
                        (uu, vv),
                        uv.last().copied(),
                        ring_periods,
                    );
                    leg_shift = Some(shift);
                    if shift == (0.0, 0.0) {
                        boundary_first_fold::lift_to_chain(
                            (uu, vv),
                            uv.last().copied(),
                            ring_periods,
                        )
                    } else {
                        uv.push((uu, vv));
                        pos.push(boundary_first_fold::ring_point(target));
                        sources.push((eh, f, f));
                        (uu + shift.0, vv + shift.1)
                    }
                }
            };
            first_leg_starts.entry(eh).or_insert((uu, vv));
            uv.push((uu, vv));
            pos.push(boundary_first_fold::ring_point(target)); // PK-344 signed zeros
            let next_idx = if forward { idx + 1 } else { idx - 1 };
            sources.push((eh, f, samples.fractions[next_idx]));
            walk_forward = pcurve_forward_along_walk;
        }
        let end = if forward { m - 1 } else { 0 };
        let geometry = &sof_ring_geometry(brep, data)?.geometry;
        let last = uv.last().copied();
        if let Some(image) = boundary_first_fold::collapsed_row_arrival(
            geometry,
            pc,
            walk_forward,
            last,
            leg_shift,
            ring_periods,
        ) {
            uv.push(image);
            pos.push(boundary_first_fold::ring_point(samples.points[end]));
            sources.push((eh, samples.fractions[end], samples.fractions[end]));
        }
        he = data.next;
        if he == start {
            break;
        }
    }
    boundary_first_fold::weld_consecutive_positions(&mut pos);
    let mut ring = BoundaryRing { uv, pos, sources };
    boundary_first_fold::drop_repeated_samples(&mut ring);
    boundary_first_fold::close_through_collapsed_runs(&mut ring, ring_periods);
    Some(ring)
}

/// Even-odd point-in-region across all rings (outer + holes).
fn point_in_region(p: (f64, f64), rings: &[BoundaryRing]) -> bool {
    let mut inside = false;
    for ring in rings {
        let n = ring.uv.len();
        for i in 0..n {
            let (x1, y1) = ring.uv[i];
            let (x2, y2) = ring.uv[(i + 1) % n];
            if (y1 > p.1) != (y2 > p.1) {
                let x_int = x1 + (p.1 - y1) / (y2 - y1) * (x2 - x1);
                if p.0 < x_int {
                    inside = !inside;
                }
            }
        }
    }
    inside
}

/// Metric distance from a raw-UV point to the nearest ring segment:
/// both are scaled by (su, sv) so the clearance comparison is in
/// millimetres regardless of the chart's parameter units.
fn dist_to_rings(p: (f64, f64), rings: &[BoundaryRing], su: f64, sv: f64) -> f64 {
    let p = (p.0 * su, p.1 * sv);
    let mut best = f64::INFINITY;
    for ring in rings {
        let n = ring.uv.len();
        for i in 0..n {
            let a = (ring.uv[i].0 * su, ring.uv[i].1 * sv);
            let b = (ring.uv[(i + 1) % n].0 * su, ring.uv[(i + 1) % n].1 * sv);
            let (dx, dy) = (b.0 - a.0, b.1 - a.1);
            let l2 = dx * dx + dy * dy;
            let t = if l2 < 1e-300 {
                0.0
            } else {
                (((p.0 - a.0) * dx + (p.1 - a.1) * dy) / l2).clamp(0.0, 1.0)
            };
            let (qx, qy) = (a.0 + dx * t, a.1 + dy * t);
            best = best.min(((p.0 - qx).powi(2) + (p.1 - qy).powi(2)).sqrt());
        }
    }
    best
}
