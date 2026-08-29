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
//! Observability: set `POPKERN_TRACE=boundary-first` to log, per face falling back
//! to the uniform-grid tessellator, WHY it fell back (missing pcurve, no
//! canonical samples, off-surface residual, uncovered ring segment,
//! periodic full-wrap ring).

use std::collections::HashMap;

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::TessellationOptions;
use crate::mesh::{FaceMesh, MeshVertex};

#[path = "boundary_first_delaunay.rs"]
mod boundary_first_delaunay;
#[path = "boundary_first_edge_samples.rs"]
mod boundary_first_edge_samples;
#[path = "boundary_first_revolution.rs"]
mod boundary_first_revolution;
#[path = "boundary_first_revolution_bipolar.rs"]
mod boundary_first_revolution_bipolar;
#[path = "boundary_first_ruled_strip.rs"]
mod boundary_first_ruled_strip;
#[path = "boundary_first_slit.rs"]
mod boundary_first_slit;
#[path = "boundary_first_unwrap.rs"]
mod boundary_first_unwrap;
use boundary_first_delaunay::bowyer_watson;
pub(crate) use boundary_first_edge_samples::{EdgeSamples, sample_all_edges};

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
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());
    for lh in loops {
        let built = build_ring(brep, lh, edge_samples)?;
        // Slit dissolution (PK-335): coincident double-traversals
        // (shell rim slit arms) are chart artifacts — cut the ring
        // into its genuine boundary components before charting, or
        // Delaunay coverage chases the zero-width arms forever.
        for mut ring in boundary_first_slit::dissolve_slits(built, {
            use popkern_core::geometry::surface::SurfaceEval as _;
            surface.parameter_periodicity()
        }) {
            // Periodic surfaces (analytic cylinders: u = angle) can have
            // pcurves authored in different period branches; a ring that
            // crosses the seam jumps ~2π in the UV chart and Delaunay
            // cannot cover the jump segment. Unwrap u to the branch
            // nearest each previous sample (positions are canonical 3D
            // points and unaffected; cos/sin accept any branch).
            {
                use popkern_core::geometry::surface::SurfaceEval as _;
                let (period_u, period_v) = surface.parameter_periodicity();
                // A sample AT a chart pole (sphere pole, cone apex) has
                // undefined azimuth: shifting u by half a period leaves
                // the 3D point fixed. Its pcurve u is arbitrary, and if
                // it participates in the consecutive unwrap it chains
                // the rest of the ring into a spurious ±period branch —
                // a phantom "full wrap" (flagship corner caps hold one
                // vertex exactly at the cap sphere's pole; the phantom
                // winding declined 2 of 4 caps). Degenerate samples are
                // skipped as unwrap anchors and their u is rewritten
                // from the flanking regular samples afterwards.
                let ring_scale = ring
                    .pos
                    .iter()
                    .zip(ring.pos.iter().cycle().skip(1))
                    .map(|(a, b)| a.distance_to(b))
                    .fold(0.0_f64, f64::max)
                    .max(1e-9);
                // Gate at 1e-4 of the ring scale: a sample whose
                // half-period u-shift moves it by less than a
                // ten-thousandth of a boundary chord has no meaningful
                // azimuth for charting, and rewriting its u moves
                // geometry sub-chord. A pole VERTEX projected at
                // v ≈ 1e-6 shifts by ~2·r·v ≈ 1e-8 — above the old
                // 1e-7·scale gate — and its arbitrary atan2 u chained
                // PK-325's corner-cap rings into a phantom full wrap
                // (all four wavy-plate top caps declined to uniform-grid
                // sampling and their chains tore off the strips).
                let u_degenerate: Vec<bool> = match period_u {
                    Some(p) => ring
                        .uv
                        .iter()
                        .zip(ring.pos.iter())
                        .map(|((u, v), pos)| {
                            surface.point_at(u + 0.5 * p, *v).distance_to(pos) < 1.0e-4 * ring_scale
                        })
                        .collect(),
                    None => Vec::new(),
                };
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
                    if popkern_core::trace::trace_enabled("boundary-first") {
                        eprintln!("    [bf] face {face_handle:?} u-wrap decline");
                    }
                    boundary_first_unwrap::dump_ring(&ring, "u");
                    return None;
                }
                if let Some(p) = period_v
                    && !boundary_first_unwrap::unwrap_dir(&mut ring, p, |uv| &mut uv.1, &[])
                {
                    boundary_first_unwrap::dump_ring(&ring, "v");
                    return None;
                }
            }
            if ring.uv.len() >= 3 {
                rings.push(ring);
            }
        }
    }
    if rings.is_empty() {
        return None;
    }
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
    // structured tensor path: exact boundary conformity, interior on
    // the cross-station lattice, chord error = the edge sampling's
    // own tolerance. Declines (the boundary-first trace topic reason) continue below.
    if let Some(fm) = boundary_first_ruled_strip::try_tessellate(
        face_handle,
        surface,
        orientation,
        &rings,
        options,
    ) {
        return Some(fm);
    }

    // Collinear boundary runs (straight edges sampled at several
    // stations) sit exactly on one line in UV; Bowyer–Watson then
    // connects the run's endpoints directly and the interior
    // stations never join the hull — the coverage check below would
    // reject the whole face. Nudge each strictly-collinear boundary
    // point INWARD by an invisible epsilon in UV (3D positions stay
    // canonical, so shared-edge closure is unaffected).
    for ring in &mut rings {
        let n = ring.uv.len();
        if n < 3 {
            continue;
        }
        let mut scale = 0.0_f64;
        for i in 0..n {
            let a = ring.uv[i];
            let b = ring.uv[(i + 1) % n];
            scale = scale.max((a.0 - b.0).abs().max((a.1 - b.1).abs()));
        }
        let eps = scale.max(1e-9) * 1.0e-7;
        let mut area2 = 0.0;
        for i in 0..n {
            let (x1, y1) = ring.uv[i];
            let (x2, y2) = ring.uv[(i + 1) % n];
            area2 += x1 * y2 - x2 * y1;
        }
        let inward = if area2 > 0.0 { 1.0 } else { -1.0 };
        let orig: Vec<(f64, f64)> = ring.uv.clone();
        for i in 0..n {
            let p = orig[i];
            let q = orig[(i + n - 1) % n];
            let r = orig[(i + 1) % n];
            let cross = (p.0 - q.0) * (r.1 - q.1) - (p.1 - q.1) * (r.0 - q.0);
            let span = ((r.0 - q.0).powi(2) + (r.1 - q.1).powi(2))
                .sqrt()
                .max(1e-300);
            if (cross / span).abs() < eps {
                // Left normal of q→r is inward for CCW rings.
                let (dx, dy) = ((r.0 - q.0) / span, (r.1 - q.1) / span);
                ring.uv[i] = (p.0 - -(dy * eps * inward), p.1 + -(dx * eps * inward));
            }
        }
    }

    // Deterministic micro-jitter on every boundary point's UV
    // (≈1e-7 of the ring scale, direction hashed from the index):
    // exact COCIRCULAR point groups — font-glyph arcs flatten into
    // them by construction — read as degenerate to the incircle
    // test the same way exactly-collinear runs do, and cascade
    // coverage holes through the Delaunay. The jitter is orders of
    // magnitude below any chord tolerance; 3D positions stay
    // canonical, so shared-edge closure is unaffected.
    {
        let mut scale = 0.0_f64;
        for ring in &rings {
            let n = ring.uv.len();
            for i in 0..n {
                let a = ring.uv[i];
                let b = ring.uv[(i + 1) % n];
                scale = scale.max((a.0 - b.0).abs().max((a.1 - b.1).abs()));
            }
        }
        let eps = scale.max(1e-9) * 1.0e-7;
        let mut k = 0u64;
        for ring in &mut rings {
            for uv in &mut ring.uv {
                // Cheap deterministic hash → angle.
                k = k
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                let angle = (k >> 11) as f64 / (1u64 << 53) as f64 * std::f64::consts::TAU;
                uv.0 += eps * angle.cos();
                uv.1 += eps * angle.sin();
            }
        }
    }

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
    let n = options.min_subdivisions.max(1);
    // Refine like the uniform-grid path would for curved surfaces: sample
    // the surface diagonal and bump density if it deviates from the
    // chord more than chord_height.
    let mut n = n;
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
    }
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
    let (su, sv) = {
        let (um, vm) = (0.5 * (u0 + u1), 0.5 * (v0 + v1));
        let p = surface.point_at(um, vm);
        let eps_u = ((u1 - u0).abs() * 1e-4).max(1e-9);
        let eps_v = ((v1 - v0).abs() * 1e-4).max(1e-9);
        let du = surface.point_at(um + eps_u, vm).distance_to(&p) / eps_u;
        let dv = surface.point_at(um, vm + eps_v).distance_to(&p) / eps_v;
        (du.max(1e-9), dv.max(1e-9))
    };
    let cell = ((su * (u1 - u0)) / n as f64)
        .max((sv * (v1 - v0)) / n as f64)
        .max(1e-12);
    let clearance = cell / 3.0;
    // LONG boundary segments on CURVED charts must be refined before
    // the grid: the diametral-disk exclusion below evacuates every
    // grid point within a segment's protecting disk, and a 2-sample
    // straight edge (a cylinder ruling or seam — zero 3D chord
    // height, so adaptive edge sampling never subdivides it) owns a
    // disk covering much of the chart. The boundary-only Delaunay
    // that results is legal in UV but its triangles are 3D CHORDS
    // through the solid (measured: the breached-bore wall's vertical
    // trim anchored a membrane capping the void, mesh volume +0.7%;
    // the column wall's seam evacuated its ENTIRE 2π×60 chart).
    // Ruppert's protection is split-the-segment, not evacuate-the-
    // interior: request splits until segments are ~cell-sized; the
    // rounds in `tessellate` re-run all faces with the shared
    // refined samples, so watertightness is preserved. Planar charts
    // skip this — their chords lie on the surface.
    if !matches!(
        surface,
        popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }
    ) {
        for ring in &rings {
            let m = ring.uv.len();
            for k in 0..m {
                let a = ring.uv[k];
                let b = ring.uv[(k + 1) % m];
                let len2 = (su * (b.0 - a.0)).powi(2) + (sv * (b.1 - a.1)).powi(2);
                if len2 > (1.5 * cell) * (1.5 * cell) {
                    let (eh, f_lo, f_hi) = ring.sources[k];
                    split_requests.push((eh, 0.5 * (f_lo + f_hi)));
                }
            }
        }
    }
    // A boundary segment survives unconstrained Delaunay if its
    // DIAMETRAL disk is point-free (the classic protecting-ball
    // condition); a plain distance clearance admits grid points that
    // sit outside the clearance band yet inside the diametral disk
    // of a diagonal segment, and the coverage check then rejects the
    // whole face (measured on a plate face with a circular hole
    // ring).
    let in_any_diametral_disk = |u: f64, v: f64| -> bool {
        let (x, y) = (u * su, v * sv);
        for ring in &rings {
            let m = ring.uv.len();
            for k in 0..m {
                let a = (ring.uv[k].0 * su, ring.uv[k].1 * sv);
                let b = (ring.uv[(k + 1) % m].0 * su, ring.uv[(k + 1) % m].1 * sv);
                let cx = 0.5 * (a.0 + b.0);
                let cy = 0.5 * (a.1 + b.1);
                let r2 = 0.25 * ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2));
                if (x - cx).powi(2) + (y - cy).powi(2) < r2 * 1.0201 {
                    return true;
                }
            }
        }
        false
    };
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
    ) && {
        let u_mid = 0.5 * (u0 + u1);
        let pa = surface.point_at(u_mid, v0);
        let pb = surface.point_at(u_mid, v1);
        let pm = surface.point_at(u_mid, 0.5 * (v0 + v1));
        let chord_mid = Point3::new(
            0.5 * (pa.x + pb.x),
            0.5 * (pa.y + pb.y),
            0.5 * (pa.z + pb.z),
        );
        pm.distance_to(&chord_mid) <= options.chord_height
    };
    if ruled_v {
        // RULED charts (cylinder/cone walls, extruded NURBS): the
        // interior lattice reuses the BOUNDARY's u-stations — the
        // adaptive edge sampling already made those as dense as the
        // curvature demands, while the bbox grid is capped at
        // max_subdivisions and its coarser u-spacing carries 16× the
        // chord error (measured: 4.6e-3 vs a 2e-3 volume gate on the
        // chamfered column). Along the ruling v is STRAIGHT, so a few
        // rows cost no accuracy and give the Delaunay interior
        // support against trim bites (measured: the breached-bore
        // wall membraned with boundary-only triangulation even after
        // the trim splits).
        let mut u_stations: Vec<f64> = rings
            .iter()
            .flat_map(|r| r.uv.iter().map(|&(u, _)| u))
            .collect();
        u_stations.sort_by(f64::total_cmp);
        // Dedup PROPORTIONALLY to the chart span: the two rims invert
        // to near-identical u (differing ~1e-8 from independent
        // atan2 of separately-sampled 3D points), and near-duplicate
        // interior points corrupt the Bowyer-Watson cavity
        // (measured: non-Delaunay 0.63-rad chords beside the seam
        // with a dense kept lattice — the PK-146 seam band).
        let station_merge = ((u1 - u0).abs() * 1.0e-7).max(1e-9);
        u_stations.dedup_by(|a, b| (*a - *b).abs() < station_merge);
        if u_stations.len() > 512 {
            let step = u_stations.len() as f64 / 512.0;
            u_stations = (0..512)
                .map(|i| u_stations[(i as f64 * step) as usize])
                .collect();
        }
        let n_v = options.min_subdivisions.max(4);
        let seam_trace =
            popkern_core::trace::trace_enabled("boundary-first-seam") && u_stations.len() > 100;
        let mut kill = [0usize; 4]; // region, clearance, disk, kept — seam band only
        // Interior v-stations CONFORM to the boundary's v-cadence:
        // with a fixed sparse row count, ruling-side boundary samples
        // (a refined seam at half the row spacing) sit exactly midway
        // between rows and the whole seam neighborhood becomes one
        // near-cocircular family — every triangulation is legal and
        // Bowyer-Watson picks long-u fans through the empty inter-row
        // band (measured: 0.63-rad chords beside the bore seam,
        // volume +1.3e-3 at any chord tolerance, ZERO circumcircle
        // violations — the fans are legal). Taking the rows from the
        // ring's own interior v-values tiles the lattice conformally
        // with the boundary; the fixed grid remains the fallback for
        // boundaries with no interior v-samples.
        let mut v_stations: Vec<f64> = rings
            .iter()
            .flat_map(|r| r.uv.iter().map(|&(_, v)| v))
            .filter(|&v| v > v0 + 1e-9 && v < v1 - 1e-9)
            .collect();
        v_stations.sort_by(f64::total_cmp);
        let v_merge = ((v1 - v0).abs() * 1.0e-7).max(1e-9);
        v_stations.dedup_by(|a, b| (*a - *b).abs() < v_merge);
        if v_stations.len() < n_v - 1 {
            v_stations = (1..n_v)
                .map(|j| v0 + (v1 - v0) * j as f64 / n_v as f64)
                .collect();
        } else if v_stations.len() > 64 {
            let step = v_stations.len() as f64 / 64.0;
            v_stations = (0..64)
                .map(|i| v_stations[(i as f64 * step) as usize])
                .collect();
        }
        for &u in &u_stations {
            if u <= u0 || u >= u1 {
                continue;
            }
            for &v in &v_stations {
                let in_band = u - u0 < 0.8 || u1 - u < 0.8;
                let detail = seam_trace && u1 - u < 0.9 && (v - 1.0).abs() < 0.01;
                if !point_in_region((u, v), &rings) {
                    if seam_trace && in_band {
                        kill[0] += 1;
                    }
                    if detail {
                        eprintln!("    [bf-seam] ({u:.3},{v:.3}) KILLED region");
                    }
                } else if dist_to_rings((u, v), &rings, su, sv) <= clearance {
                    if seam_trace && in_band {
                        kill[1] += 1;
                    }
                    if detail {
                        eprintln!(
                            "    [bf-seam] ({u:.3},{v:.3}) KILLED clearance d={:.4}",
                            dist_to_rings((u, v), &rings, su, sv)
                        );
                    }
                } else if in_any_diametral_disk(u, v) {
                    if seam_trace && in_band {
                        kill[2] += 1;
                    }
                    if detail {
                        eprintln!("    [bf-seam] ({u:.3},{v:.3}) KILLED disk");
                    }
                } else {
                    if seam_trace && in_band {
                        kill[3] += 1;
                    }
                    if detail {
                        eprintln!("    [bf-seam] ({u:.3},{v:.3}) KEPT");
                    }
                    uv_pts.push((u, v));
                    pos_pts.push(None);
                }
            }
        }
        if seam_trace {
            eprintln!(
                "  [bf-seam] stations={} n_v={n_v} cell={cell:.4} clearance={clearance:.4} band kills: region={} clearance={} disk={} kept={}",
                u_stations.len(),
                kill[0],
                kill[1],
                kill[2],
                kill[3]
            );
        }
    } else {
        for i in 1..n {
            for j in 1..n {
                let u = u0 + (u1 - u0) * i as f64 / n as f64;
                let v = v0 + (v1 - v0) * j as f64 / n as f64;
                if point_in_region((u, v), &rings)
                    && dist_to_rings((u, v), &rings, su, sv) > clearance
                    && !in_any_diametral_disk(u, v)
                {
                    uv_pts.push((u, v));
                    pos_pts.push(None);
                }
            }
        }
    }

    if popkern_core::trace::trace_enabled("boundary-first-ring") {
        let boundary_n: usize = ring_ranges.iter().map(|&(s, e)| e - s).sum();
        eprintln!(
            "    [bf-grid] face {face_handle:?}: boundary {boundary_n} interior {} (n={n} cell={cell:.4} clearance={clearance:.4})",
            uv_pts.len() - boundary_n
        );
    }

    // Delaunay on METRIC-SCALED coordinates: the chart mixes units
    // (radians × millimetres on a cylinder — 1 rad ≈ 35 mm on the
    // column), and the incircle predicate degrades continuously with
    // anisotropy (measured on the full-period wall: aspect ~78 gave
    // a sane mesh, ~312 lost 30% of the wedge volume, ~1250 went
    // NEGATIVE — inside-out junk). Scaling u,v by |∂S/∂u|,|∂S/∂v| at
    // the chart centre makes the triangulation topology follow the
    // true surface metric; (su, sv) is the same scale the protection
    // filters above use. Stored vertices stay in raw UV.
    let bw_pts: Vec<(f64, f64)> = uv_pts.iter().map(|&(u, v)| (u * su, v * sv)).collect();
    let tris = bowyer_watson(&bw_pts)?;
    boundary_first_edge_samples::delaunay_audit(face_handle, &bw_pts, &tris);

    // Filter to the trim region and enforce winding.
    let mut triangles: Vec<[usize; 3]> = Vec::new();
    for t in tris {
        let (a, b, c) = (uv_pts[t[0]], uv_pts[t[1]], uv_pts[t[2]]);
        let cx = (a.0 + b.0 + c.0) / 3.0;
        let cy = (a.1 + b.1 + c.1) / 3.0;
        if !point_in_region((cx, cy), &rings) {
            continue;
        }
        let signed = (b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1);
        if signed.abs() < 1e-16 {
            continue;
        }
        if signed > 0.0 {
            triangles.push(t);
        } else {
            triangles.push([t[0], t[2], t[1]]);
        }
    }
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
            if !edge_set.contains(&(a.min(b), a.max(b))) {
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
    if popkern_core::trace::trace_enabled("boundary-first-ring") {
        let mut vol6 = 0.0;
        for t in &triangles {
            let (a, b, c) = (
                vertices[t[0]].position,
                vertices[t[1]].position,
                vertices[t[2]].position,
            );
            vol6 += a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                + a.z * (b.x * c.y - b.y * c.x);
        }
        let mut uv_area2 = 0.0;
        let mut ring_area2 = 0.0;
        for t in &triangles {
            let (a, b, c) = (vertices[t[0]].uv, vertices[t[1]].uv, vertices[t[2]].uv);
            uv_area2 += ((b.0 - a.0) * (c.1 - a.1) - (c.0 - a.0) * (b.1 - a.1)).abs();
        }
        for ring in &rings {
            let m = ring.uv.len();
            let mut s = 0.0;
            for i in 0..m {
                let (x1, y1) = ring.uv[i];
                let (x2, y2) = ring.uv[(i + 1) % m];
                s += x1 * y2 - x2 * y1;
            }
            ring_area2 += s;
        }
        eprintln!(
            "    [bf-vol] face {face_handle:?} tris {} wedge6 {vol6:.1} uv_area2 {uv_area2:.2} ring_area2 {:.2}",
            triangles.len(),
            ring_area2.abs()
        );
    }
    Some(FaceMesh {
        face: face_handle,
        vertices,
        triangles,
    })
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
    let mut uv = Vec::new();
    let mut pos = Vec::new();
    let mut sources = Vec::new();
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

        // Direction: does the canonical polyline start at this
        // half-edge's origin? For CLOSED edges (both endpoints are
        // the same point) that test is degenerate — use the pcurve's
        // walk direction instead: its early image must match the
        // polyline's second sample rather than its second-to-last.
        let closed_edge = samples
            .points
            .first()
            .zip(samples.points.last())
            .map(|(a, b)| a.distance_to(b) < 1e-9)
            .unwrap_or(false);
        let forward = if closed_edge && samples.points.len() >= 3 {
            let sof_dir = brep.faces[brep.loops[data.loop_].face].surface.as_ref()?;
            let (pt0, pt1) = pc.parameter_range();
            let early_uv = pc.point_at(pt0 + (pt1 - pt0) * 0.05);
            let early = sof_dir.geometry.point_at(early_uv.x, early_uv.y);
            early.distance_to(&samples.points[1])
                <= early.distance_to(&samples.points[samples.points.len() - 2])
        } else {
            samples.points.first()?.distance_to(&origin)
                <= samples.points.last()?.distance_to(&origin)
        };
        let m = samples.points.len();
        if popkern_core::trace::trace_enabled("boundary-first-ring") {
            eprintln!("      [bf-ring] he={he:?} edge={eh:?} m={m} closed={closed_edge}");
        }
        // Emit all but the last point (the next half-edge provides
        // its own start, which is the shared vertex).
        for k in 0..(m - 1) {
            let idx = if forward { k } else { m - 1 - k };
            let f = samples.fractions[idx];
            // Pcurve param at the same fraction of ITS window,
            // mirrored when the pcurve runs against the walk
            // direction. The pcurve's own direction is detected from
            // its endpoint images.
            let p_start = pc.point_at(pc.t_start);
            let p_end = pc.point_at(pc.t_end);
            let sof = brep.faces[brep.loops[data.loop_].face].surface.as_ref()?;
            let s_start = sof.geometry.point_at(p_start.x, p_start.y);
            let s_end = sof.geometry.point_at(p_end.x, p_end.y);
            let pcurve_forward_along_walk =
                s_start.distance_to(&origin) <= s_end.distance_to(&origin);
            // Fraction along the WALK is k/(m-1); convert to the
            // curve fraction (f) orientation… the pcurve maps the
            // same geometric edge, so use the walk fraction with the
            // pcurve's own direction.
            let walk_f = k as f64 / (m - 1) as f64;
            let tp = if pcurve_forward_along_walk {
                pc.t_start + (pc.t_end - pc.t_start) * walk_f
            } else {
                pc.t_end - (pc.t_end - pc.t_start) * walk_f
            };
            // The pcurve is only a SEED: refine each UV by point
            // inversion of the canonical 3D sample so imprecise
            // pcurves (straight chords for bulging curves) still
            // yield exact boundary UVs. If the refined image stays
            // far from the canonical point, the curve genuinely does
            // not lie on this surface (sphere seam pole chord) —
            // fall back to the uniform-grid path.
            let uvp = pc.point_at(tp);
            let sof = brep.faces[brep.loops[data.loop_].face].surface.as_ref()?;
            let target = samples.points[idx];
            let (mut uu, mut vv) = (uvp.x, uvp.y);
            // Analytic cylinders: invert analytically (atan2 in the
            // surface's own basis) — the generic Newton projection
            // wanders across period branches and reversed-branch
            // pcurve seeds walked rings the long way around the seam
            // (span 2π). Canonical branch here; the ring-level
            // unwrap restores continuity.
            if let popkern_core::geometry::SurfaceGeometry::Cylinder { origin, axis, .. } =
                &sof.geometry
            {
                let (perp1, perp2) =
                    popkern_core::geometry::analytic_surfaces::orthonormal_basis(*axis);
                let d = target - *origin;
                (uu, vv) = (d.dot(&perp2).atan2(d.dot(&perp1)), d.dot(axis));
            } else if sof.geometry.point_at(uu, vv).distance_to(&target) > 1.0e-9
                && let Some((ru, rv)) =
                    popkern_core::geometry::uv_projection::project_point_to_uv_near(
                        &sof.geometry,
                        &target,
                        uu,
                        vv,
                    )
            {
                (uu, vv) = (ru, rv);
            }
            let mut residual = sof.geometry.point_at(uu, vv).distance_to(&target);
            let scale = samples
                .points
                .first()?
                .distance_to(samples.points.last()?)
                .max(1e-9);
            let gate = 1.0e-3 * scale.max(1.0);
            if residual > gate {
                // The pcurve seed can be DISPLACED (a setback fan's
                // stale rim pcurves sit 1.8e-3 off) and the seeded
                // local projection then converges near the displaced
                // location, reading 1.4e-3 where the true nearest
                // point is 5e-5 away — the face fell back to the
                // uniform-grid path and the mesh cracked. Retry with a COLD
                // global projection under the SAME gate before giving
                // up; the ring-level periodic unwrap downstream
                // restores branch continuity exactly as it does for
                // analytic inversions.
                if let Some((gu, gv)) = popkern_core::geometry::uv_projection::project_point_to_uv(
                    &sof.geometry,
                    &target,
                ) {
                    let g_res = sof.geometry.point_at(gu, gv).distance_to(&target);
                    if g_res < residual {
                        (uu, vv) = (gu, gv);
                        residual = g_res;
                    }
                }
            }
            // NOTE: widening this gate by recorded entity tolerance
            // was tried twice (with and without the cold retry above)
            // and is FALSIFIED both times: it repartitions the
            // boundary-first/uniform-grid face split and the mixed
            // boundaries crack far worse than the fallback it avoids
            // (hexagon6_blunt 0 -> 482; pentagon5_blunt/G2
            // 106 -> 492). Do not reintroduce it.
            if residual > gate {
                if popkern_core::trace::trace_enabled("boundary-first") {
                    eprintln!("    [bf] edge {eh:?}: off-surface residual {residual:.3e}");
                }
                return None;
            }
            uv.push((uu, vv));
            // Normalize negative zeros: two arcs meeting at an anchor
            // evaluate the shared point as +0.0 / −0.0 in a coordinate,
            // the position-weld keys are BIT patterns, and the mismatch
            // left an unpaired sliver at the joint (PK-344).
            pos.push(popkern_core::geometry::Point3::new(
                target.x + 0.0,
                target.y + 0.0,
                target.z + 0.0,
            ));
            let next_idx = if forward { idx + 1 } else { idx - 1 };
            sources.push((eh, f, samples.fractions[next_idx]));
        }
        he = data.next;
        if he == start {
            break;
        }
    }
    // Canonicalize near-duplicate CONSECUTIVE positions (arc joints:
    // the two arcs meeting at an anchor evaluate the shared point
    // ~1e-17 apart, and the position-weld keys are BIT patterns).
    // Both entries snap to the pair's lexicographically smaller bit
    // pattern — a direction-independent choice, so the two faces
    // sharing the edges pick the SAME canonical point even though
    // they traverse in opposite orders (a keep-first dedupe diverged
    // and left an unpaired sliver at a strip-anchored rim, PK-344).
    if !pos.is_empty() {
        let scale = pos
            .iter()
            .zip(pos.iter().cycle().skip(1))
            .map(|(a, b)| a.distance_to(b))
            .fold(0.0_f64, f64::max)
            .max(1e-12);
        let weld = scale * 1e-9;
        let n = pos.len();
        for i in 0..n {
            let j = (i + 1) % n;
            let (a, b) = (pos[i], pos[j]);
            let bits =
                |p: &popkern_core::geometry::Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
            if bits(&a) != bits(&b) && a.distance_to(&b) <= weld {
                let canonical = if bits(&a) <= bits(&b) { a } else { b };
                pos[i] = canonical;
                pos[j] = canonical;
            }
        }
    }
    Some(BoundaryRing { uv, pos, sources })
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
