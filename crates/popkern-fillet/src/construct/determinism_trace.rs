//! Determinism observability (PK-085): geometry fingerprints
//! for run-to-run bisection. Enable POPKERN_DETERMINISM_TRACE=1
//! for stage-boundary fingerprints and
//! POPKERN_DETERMINISM_FACE_TRACE=1 for per-face hashes; diffing
//! two runs localizes the first divergent pipeline stage / the
//! faces a divergent builder wrote. Zero cost when unset.

use popkern_core::geometry::gregory::GregoryPatch;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::topology::BRep;

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
    std::env::var_os("POPKERN_DETERMINISM_TRACE").is_some()
}

/// Stage-boundary fingerprint print (POPKERN_DETERMINISM_TRACE=1).
pub(crate) fn determinism_trace(label: &str, brep: &BRep) {
    if std::env::var_os("POPKERN_DETERMINISM_TRACE").is_some() {
        eprintln!(
            "dtrace: {label} geo={:016x} pc={:016x} ec={:016x}",
            determinism_fingerprint(brep),
            determinism_pcurve_fingerprint(brep),
            determinism_edge_curve_fingerprint(brep)
        );
    }
    // Blend-face census at stage boundaries (PK_CENSUS_TRACE=1):
    // names the pass that merges away or re-surfaces a strip face —
    // the region path's serving debts were invisible without a
    // per-stage census.
    if std::env::var_os("PK_CENSUS_TRACE").is_some() {
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
    if std::env::var_os("PK_SLIT_TRACE").is_some() {
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

/// Per-face fingerprints (POPKERN_DETERMINISM_FACE_TRACE=1): the
/// divergent-face set between two runs names the builder that
/// wrote them. Dual-twist faces report grid and pair hashes
/// separately.
pub(crate) fn determinism_face_trace(label: &str, brep: &BRep) {
    if std::env::var_os("POPKERN_DETERMINISM_FACE_TRACE").is_none() {
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
