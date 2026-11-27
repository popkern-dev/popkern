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

use popkern_core::geometry::chart::ChartDescriptor;
use popkern_math::hash::FastHashMap as HashMap;

use popkern_core::geometry::Point3;
use popkern_core::topology::{BRep, EdgeHandle, FaceHandle};

use crate::TessellationOptions;
use crate::mesh::FaceMesh;

#[path = "boundary_first_apex.rs"]
mod boundary_first_apex;
#[path = "boundary_first_apex_sweep.rs"]
mod boundary_first_apex_sweep;
#[path = "boundary_first_cdt.rs"]
mod boundary_first_cdt;
#[path = "boundary_first_cone_development.rs"]
mod boundary_first_cone_development;
#[path = "boundary_first_delaunay.rs"]
mod boundary_first_delaunay;
#[path = "boundary_first_deviation.rs"]
mod boundary_first_deviation;
#[path = "boundary_first_fidelity.rs"]
mod boundary_first_fidelity;
#[path = "boundary_first_fold.rs"]
mod boundary_first_fold;
use boundary_first_fidelity::closed_chart_fidelity_ok;
#[path = "boundary_first_atlas.rs"]
mod boundary_first_atlas;
#[path = "boundary_first_edge_samples.rs"]
mod boundary_first_edge_samples;
#[path = "boundary_first_lattice.rs"]
mod boundary_first_lattice;
#[path = "boundary_first_mesh_rings.rs"]
mod boundary_first_mesh_rings;
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
use boundary_first_apex::close_through_apex;
pub(crate) use boundary_first_edge_samples::{EdgeSamples, sample_edges};
use boundary_first_fold::lift_to_chain;
use boundary_first_mesh_rings::mesh_rings;
use boundary_first_wrap_pair::{join_wrapped_rings, stash_wrapped_ring};

/// One boundary ring: ordered UV polygon with matched 3D positions
/// and, per point, the canonical edge + curve fraction it came from
/// (the handle encroachment refinement needs to split the RIGHT
/// shared edge span).
#[derive(Clone)]
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

/// The atlas's lattice rung for `face_handle`, a face its structured
/// meshing declined on the final edge samples: `None` when the face
/// reaches no cap (not the atlas's), `Some(None)` when the rung declined.
pub(crate) fn tessellate_face_atlas_lattice(
    brep: &BRep,
    face_handle: FaceHandle,
    edge_samples: &HashMap<EdgeHandle, EdgeSamples>,
    options: &TessellationOptions,
) -> Option<Option<FaceMesh>> {
    boundary_first_atlas::tessellate_on_atlas_lattice(brep, face_handle, edge_samples, options)
}

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
    // seam spine (a full sphere's meridian pair): with no holes its
    // chart "polygon" bounds no interior (Delaunay meshed the primitive
    // sphere 0/46), and only the atlas meshes it. With holes, and a torus
    // face's seam SQUARE, it meshes here (a drilled torus, PK-459; a ball
    // with Viviani's figure eight cut out fell back to the grid, the union
    // 1.5e-5 high, PK-1140).
    let periods = ChartDescriptor::evaluation_periods_of(surface);
    // The face's chart artifacts, from the chart service.
    let descriptor = ChartDescriptor::of(surface, &brep.tolerance);
    let seam_square = matches!(periods, (Some(_), Some(_))) && !face.inner_loops.is_empty();
    let bare_spine =
        !seam_square && face.inner_loops.is_empty() && outer_loop_is_spine(brep, face.outer_loop)?;

    // Build boundary rings from canonical samples, then chart them.
    let mut loops = vec![face.outer_loop];
    loops.extend(face.inner_loops.iter().copied());
    // A ring that cannot be built leaves the face to the structured spine
    // path below or the caller's fallback.
    let built: Option<Vec<BoundaryRing>> = loops
        .into_iter()
        .map(|lh| build_ring(brep, lh, edge_samples))
        .collect();
    // A face reaching a chart singularity meshes on the atlas (pieces in
    // charts with no singular row; chart-service.md §6a).
    if let Some(built) = &built
        && let Some(out) = boundary_first_atlas::tessellate_on_atlas(
            brep,
            face_handle,
            (surface, orientation),
            built,
            (
                &descriptor,
                options,
                boundary_first_atlas::PieceMeshing::Structured,
            ),
            split_requests,
        )
    {
        return out;
    }
    if bare_spine {
        // A spine between two apexes (a spindle's points: geometric
        // singularities the atlas leaves to the face chart) meshes as a
        // two-apex revolution (PK-339); one between chart poles is the
        // atlas's, and what it declined falls to the grid.
        let apexes = !descriptor.singular_rows.is_empty()
            && descriptor
                .singular_rows
                .iter()
                .all(|r| r.kind == popkern_core::geometry::chart::SingularKind::Apex);
        if apexes {
            return boundary_first_revolution_bipolar::seam_spine_revolution(
                brep,
                face_handle,
                edge_samples,
                options,
            );
        }
        if popkern_core::trace::trace_enabled("atlas") {
            eprintln!("    [atlas] face {face_handle:?}: a bare seam spine the atlas declined");
        }
        return None;
    }
    let built = built?;
    // A cone face sweeping round its apex close beside it meshes in the
    // cone's development (an isometry; the face chart's azimuth is
    // stretched without bound there). A decline keeps the face chart.
    let curve_at = |e: EdgeHandle, f: f64| {
        let c = brep.edges.get(e)?.curve.as_ref()?;
        let eval: &dyn popkern_core::geometry::curve::CurveEval = &c.geometry;
        Some(eval.point_at(c.t_start + (c.t_end - c.t_start) * f))
    };
    if let Some((dev, rings, azimuths)) =
        boundary_first_cone_development::development_for(surface, &built, &curve_at)
        && let Some(fm) = boundary_first_cone_development::mesh_developed(
            face_handle,
            (surface, orientation),
            (&descriptor, options),
            (&dev, &rings, &azimuths),
            split_requests,
        )
    {
        return Some(fm);
    }
    let rings = match chart_rings(
        (face_handle, surface, orientation),
        built,
        (periods, None, &descriptor),
        options,
    )? {
        Charted::Rings(rings) => rings,
        Charted::Mesh(fm) => return Some(fm),
    };
    mesh_rings(
        (face_handle, surface, orientation),
        rings,
        options,
        split_requests,
        &descriptor,
        None,
    )
}

/// Whether a loop traverses every one of its edges at least twice (a seam
/// spine: a full sphere's meridian pair). `None` on a broken loop.
fn outer_loop_is_spine(
    brep: &BRep,
    loop_handle: popkern_core::topology::LoopHandle,
) -> Option<bool> {
    let mut edge_counts: HashMap<EdgeHandle, usize> = HashMap::default();
    let start = brep.loops.get(loop_handle)?.half_edge;
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
    Some(!edge_counts.is_empty() && edge_counts.values().all(|&c| c >= 2))
}

/// A face's rings charted for the boundary-first pipeline, or a mesh a
/// structured path served outright.
enum Charted {
    Rings(Vec<BoundaryRing>),
    Mesh(FaceMesh),
}

/// Chart a face's built rings (one per loop, ring 0 the outer): slits
/// dissolved, each lifted across the period cut, full wraps joined into
/// the outer band, pinches split. `None` when a wrap cannot be joined.
fn chart_rings(
    (face_handle, surface, orientation): (
        FaceHandle,
        &popkern_core::geometry::SurfaceGeometry,
        bool,
    ),
    built_rings: Vec<BoundaryRing>,
    (periods, slit_azimuth, descriptor): (
        (Option<f64>, Option<f64>),
        Option<f64>,
        &ChartDescriptor,
    ),
    options: &TessellationOptions,
) -> Option<Charted> {
    let mut rings: Vec<BoundaryRing> = Vec::new();
    let (mut wrapped_rings, mut wrapped_period) = (Vec::<BoundaryRing>::new(), None::<(f64, bool)>); // PK-468
    let sole_loop = built_rings.len() == 1;
    for built in built_rings {
        // Slit dissolution (PK-335): coincident double-traversals
        // (shell rim slit arms) are chart artifacts — cut the ring
        // into its genuine boundary components before charting, or
        // Delaunay coverage chases the zero-width arms forever.
        for mut ring in boundary_first_slit::dissolve_slits(built, periods) {
            // Periodic surfaces (analytic cylinders: u = angle) can have
            // pcurves authored in different period branches; a ring that
            // crosses the seam jumps ~2π in the UV chart and Delaunay
            // cannot cover the jump segment. Unwrap u to the branch
            // nearest each previous sample (positions are canonical 3D
            // points and unaffected; cos/sin accept any branch).
            {
                let (period_u, period_v) = periods;
                // Samples at an apex (a geometric singularity: a cone's
                // point), from the chart service. A chart pole never
                // reaches here: a face reaching one meshes on the atlas,
                // its pieces in charts with no singular row.
                let mut u_degenerate =
                    boundary_first_unwrap::at_apex(&ring, surface, period_u, descriptor);
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
                    if sole_loop
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
                        return Some(Charted::Mesh(fm));
                    }
                    match close_through_apex(ring, &mut u_degenerate, surface, orientation, p) {
                        Ok(lobes) => rings.extend(lobes), // PK-931 apex row
                        Err(ring) => {
                            stash_wrapped_ring(ring, &mut wrapped_rings); // PK-468 band slit
                            boundary_first_wrap_pair::record_wrap(&mut wrapped_period, (p, false))?;
                        }
                    }
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
            (&wrapped_rings, &rings),
            p,
            surface,
            orientation,
            (options.chord_height, slit_azimuth),
            transposed,
        )?;
        rings.splice(0..0, joined);
    }
    Some(Charted::Rings(rings))
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
    let mut first_leg_starts: HashMap<EdgeHandle, (f64, f64)> = HashMap::default();
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
                Some(_) => lift_to_chain((uu, vv), uv.last().copied(), ring_periods),
                None => {
                    let shift = boundary_first_fold::fold_shift_for_leg(
                        &sof.geometry,
                        first_leg_starts.get(&eh).copied(),
                        (uu, vv),
                        uv.last().copied(),
                        ring_periods,
                    );
                    leg_shift = Some(shift);
                    if shift == (0.0, 0.0)
                        && boundary_first_fold::starts_unfolded(
                            &sof.geometry,
                            first_leg_starts.get(&eh).copied(),
                            (uu, vv),
                            ring_periods,
                        )
                    {
                        (uu, vv)
                    } else if shift == (0.0, 0.0) {
                        lift_to_chain((uu, vv), uv.last().copied(), ring_periods)
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
