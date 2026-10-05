//! Determinism observability (PK-085): geometry fingerprints
//! for run-to-run bisection. Enable POPKERN_TRACE=determinism
//! for stage-boundary fingerprints and
//! POPKERN_TRACE=determinism-faces for per-face hashes; diffing
//! two runs localizes the first divergent pipeline stage / the
//! faces a divergent builder wrote. Zero cost when unset.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;

/// FNV-fold one geometry value into the hash.
#[inline]
fn fold(hash: &mut u64, value: f64) {
    *hash ^= value.to_bits();
    *hash = hash.wrapping_mul(0x100000001b3);
}

/// FNV-fold a Debug rendering into the hash. Rust's f64 Debug
/// prints the shortest round-trip representation, so this is
/// bit-exact for float-bearing geometry like pcurves and edge
/// curves whose internal structure varies by variant.
#[inline]
fn fold_debug(hash: &mut u64, value: &impl std::fmt::Debug) {
    for byte in format!("{value:?}").bytes() {
        *hash ^= u64::from(byte);
        *hash = hash.wrapping_mul(0x100000001b3);
    }
}

/// Deterministic fingerprint over vertices, bi-nonic Gregory
/// grids (plain and dual-twist, including candidate pairs),
/// half-edge pcurves, face-use pcurve authority, and edge
/// curves — everything the shared-trim measure reads.
pub(crate) fn determinism_fingerprint(brep: &BRep) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for (_, vertex) in &brep.vertices {
        fold(&mut hash, vertex.position.x);
        fold(&mut hash, vertex.position.y);
        fold(&mut hash, vertex.position.z);
    }
    for (_, face) in &brep.faces {
        let Some(surface) = &face.surface else {
            continue;
        };
        let SurfaceGeometry::Gregory(gp) = &surface.geometry else {
            continue;
        };
        match gp {
            GregoryPatch::RectangularBiNonic(r) => {
                for row in &r.grid {
                    for point in row {
                        fold(&mut hash, point.x);
                        fold(&mut hash, point.y);
                        fold(&mut hash, point.z);
                    }
                }
            }
            GregoryPatch::RectangularBiNonicDualTwist(r) => {
                for row in &r.grid {
                    for point in row {
                        fold(&mut hash, point.x);
                        fold(&mut hash, point.y);
                        fold(&mut hash, point.z);
                    }
                }
                for (p_u, p_v) in &r.inner_cps {
                    for point in [p_u, p_v] {
                        fold(&mut hash, point.x);
                        fold(&mut hash, point.y);
                        fold(&mut hash, point.z);
                    }
                }
            }
            _ => {}
        }
    }
    hash
}

/// Pcurve-only fingerprint (half-edge pcurves + face-use
/// authority).
pub(crate) fn determinism_pcurve_fingerprint(brep: &BRep) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for (_, half_edge) in &brep.half_edges {
        if let Some(pcurve) = &half_edge.pcurve {
            fold_debug(&mut hash, pcurve);
        }
    }
    for (key, pcurve) in &brep.face_use_pcurves {
        fold_debug(&mut hash, key);
        fold_debug(&mut hash, pcurve);
    }
    hash
}

/// Edge-curve-only fingerprint.
pub(crate) fn determinism_edge_curve_fingerprint(brep: &BRep) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    for (_, edge) in &brep.edges {
        if let Some(curve) = &edge.curve {
            fold_debug(&mut hash, curve);
        }
    }
    hash
}

/// Hash any Debug-printable value (bit-exact for float-bearing
/// geometry — see `fold_debug`). For fingerprinting pipeline
/// INPUTS (e.g. computed fillet surfaces) rather than BRep state.
pub(crate) fn determinism_value_hash(value: &impl std::fmt::Debug) -> u64 {
    let mut hash = 0xcbf29ce484222325u64;
    fold_debug(&mut hash, value);
    hash
}

/// True when stage-boundary fingerprinting is enabled.
pub(crate) fn determinism_trace_enabled() -> bool {
    popkern_core::trace::trace_enabled("determinism")
}

/// Pre-apply fingerprint of the input BRep plus every computed edge
/// geometry bundle (surface + both 3D trims), printed before the
/// apply pipeline mutates anything. Zero cost when tracing is off.
pub(crate) fn pre_apply_trace(brep: &BRep, edge_geometries: &[crate::construct::EdgeGeometry]) {
    if !determinism_trace_enabled() {
        return;
    }
    determinism_trace("pre_apply_input_brep", brep);
    for (i, eg) in edge_geometries.iter().enumerate() {
        eprintln!(
            "dtrace: pre_apply eg[{i}] edge={:?} f1={:?} f2={:?} surf={:016x} trim1={:016x} trim2={:016x}",
            eg.analysis.edge,
            eg.analysis.face_1,
            eg.analysis.face_2,
            determinism_value_hash(&eg.surface),
            determinism_value_hash(&eg.trim_3d_1),
            determinism_value_hash(&eg.trim_3d_2),
        );
    }
}

/// Stage-boundary fingerprint print (POPKERN_TRACE=determinism).
pub(crate) fn determinism_trace(label: &str, brep: &BRep) {
    if std::env::var("PK531_WATCH").is_ok() {
        let targets = [
            (0.8935211557340824, -0.07736140869204894),
            (0.9264249281536544, -0.10126739868667767),
        ];
        let near: Vec<String> = brep
            .vertices
            .iter()
            .filter(|(_, v)| {
                v.position.z.abs() < 1e-6
                    && targets
                        .iter()
                        .any(|(x, y)| (v.position.x - x).phypot(v.position.y - y) < 5e-2)
            })
            .map(|(h, v)| format!("{h:?}({:.4},{:.4})", v.position.x, v.position.y))
            .collect();
        eprintln!("[watch] {label}: {near:?}");
    }
    // POPKERN_TRACE=pcurve-fidelity: the worst same-fraction pcurve
    // deviation over the body's ANALYTIC-support faces at every stage,
    // so a stage that rewrites an edge range without re-trimming its
    // pcurves names itself (PK-396).
    if popkern_core::trace::trace_enabled("pcurve-fidelity") {
        let supports: Vec<popkern_core::topology::FaceHandle> = brep
            .faces
            .iter()
            .filter(|(_, f)| {
                matches!(
                    f.surface.as_ref().map(|s| &s.geometry),
                    Some(SurfaceGeometry::Plane { .. }) | Some(SurfaceGeometry::Cylinder { .. })
                )
            })
            .map(|(h, _)| h)
            .collect();
        match crate::pcurve_fidelity::measure_pcurve_fidelity(brep, &supports) {
            Some(fid) => eprintln!(
                "[pcurve-fidelity] {label}: worst {:.3e} edge {:?} face {:?} fraction {:.3} (measured {})",
                fid.worst, fid.worst_edge, fid.worst_face, fid.worst_fraction, fid.measured
            ),
            None => eprintln!("[pcurve-fidelity] {label}: nothing measurable"),
        }
    }
    if popkern_core::trace::trace_enabled("determinism") {
        eprintln!(
            "dtrace: {label} geo={:016x} pc={:016x} ec={:016x}",
            determinism_fingerprint(brep),
            determinism_pcurve_fingerprint(brep),
            determinism_edge_curve_fingerprint(brep)
        );
    }
    // Blend-face census at stage boundaries (POPKERN_TRACE=census):
    // names the pass that merges away or re-surfaces a strip face —
    // the region path's serving debts were invisible without a
    // per-stage census.
    if popkern_core::trace::trace_enabled("census") {
        let mut kinds: std::collections::BTreeMap<&str, usize> = std::collections::BTreeMap::new();
        let mut blends: Vec<String> = Vec::new();
        for (fh, f) in brep.faces.iter() {
            let k = match f.surface.as_ref().map(|s| &s.geometry) {
                Some(popkern_core::geometry::surface::SurfaceGeometry::Plane { .. }) => "plane",
                Some(popkern_core::geometry::surface::SurfaceGeometry::Cylinder { .. }) => "cyl",
                Some(popkern_core::geometry::surface::SurfaceGeometry::Sphere { .. }) => "sphere",
                Some(popkern_core::geometry::surface::SurfaceGeometry::Torus { .. }) => "torus",
                Some(popkern_core::geometry::surface::SurfaceGeometry::Nurbs(_)) => "nurbs",
                Some(_) => "other",
                None => "none",
            };
            *kinds.entry(k).or_default() += 1;
            if matches!(k, "nurbs" | "sphere" | "torus") {
                blends.push(format!("{fh:?}:{k}"));
            }
        }
        eprintln!("[census] {label}: {kinds:?} blends [{}]", blends.join(" "));
    }
    // Same-face ("slit") edge census at stage boundaries: an edge
    // whose two loop sides belong to one face is a wrapped-face seam
    // that boundary-first meshing declines on; the stage that mints
    // one names the pass owing a dissolve.
    if popkern_core::trace::trace_enabled("slit") {
        let mut slits = 0usize;
        for (_, e) in brep.edges.iter() {
            let Some(he) = brep.half_edges.get(e.half_edge) else {
                continue;
            };
            let Some(tw) = brep.half_edges.get(he.twin) else {
                continue;
            };
            let (Some(la), Some(lb)) = (brep.loops.get(he.loop_), brep.loops.get(tw.loop_)) else {
                continue;
            };
            if la.face == lb.face {
                slits += 1;
            }
        }
        if slits > 0 {
            eprintln!("[slit] {label}: {slits} same-face edge(s)");
        }
    }
}

/// Per-face fingerprints (POPKERN_TRACE=determinism-faces): the
/// divergent-face set between two runs names the builder that
/// wrote them. Dual-twist faces report grid and pair hashes
/// separately.
pub(crate) fn determinism_face_trace(label: &str, brep: &BRep) {
    if !popkern_core::trace::trace_enabled("determinism-faces") {
        return;
    }
    for (handle, face) in &brep.faces {
        let Some(surface) = &face.surface else {
            continue;
        };
        let SurfaceGeometry::Gregory(gp) = &surface.geometry else {
            continue;
        };
        match gp {
            GregoryPatch::RectangularBiNonic(r) => {
                let mut hash = 0xcbf29ce484222325u64;
                for row in &r.grid {
                    for point in row {
                        fold(&mut hash, point.x);
                        fold(&mut hash, point.y);
                        fold(&mut hash, point.z);
                    }
                }
                eprintln!("ftrace: {label} {handle:?} {hash:016x}");
            }
            GregoryPatch::RectangularBiNonicDualTwist(r) => {
                let mut grid_hash = 0xcbf29ce484222325u64;
                for row in &r.grid {
                    for point in row {
                        fold(&mut grid_hash, point.x);
                        fold(&mut grid_hash, point.y);
                        fold(&mut grid_hash, point.z);
                    }
                }
                let mut pair_hash = 0xcbf29ce484222325u64;
                for (p_u, p_v) in &r.inner_cps {
                    for point in [p_u, p_v] {
                        fold(&mut pair_hash, point.x);
                        fold(&mut pair_hash, point.y);
                        fold(&mut pair_hash, point.z);
                    }
                }
                eprintln!(
                    "ftrace: {label} {handle:?} grid={grid_hash:016x} pairs={pair_hash:016x}"
                );
            }
            _ => {}
        }
    }
}

/// PK-141 campaign observability (POPKERN_TRACE=trim-end): dump every
/// edge geometry's trim-curve endpoints at geometry-phase exit —
/// separator pairing, corner anchoring, and region siding all key
/// on these positions.
pub(crate) fn trace_trim_endpoints(edge_geometries: &[crate::construct::EdgeGeometry]) {
    if !popkern_core::trace::trace_enabled("trim-end") {
        return;
    }
    use popkern_core::geometry::curve::CurveEval;
    for (i, eg) in edge_geometries.iter().enumerate() {
        for (si, trim) in [&eg.trim_3d_1, &eg.trim_3d_2].into_iter().enumerate() {
            let (t0, t1) = trim.parameter_range();
            let a = trim.point_at(t0);
            let b = trim.point_at(t1);
            eprintln!(
                "[trim-end] eg{i} side{si} ({:.4},{:.4},{:.4}) -> ({:.4},{:.4},{:.4})",
                a.x, a.y, a.z, b.x, b.y, b.z
            );
        }
    }
}
