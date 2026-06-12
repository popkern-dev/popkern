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

/// Deterministic fingerprint over vertices and bi-nonic Gregory
/// grids (plain and dual-twist, including candidate pairs).
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

/// Stage-boundary fingerprint print (POPKERN_DETERMINISM_TRACE=1).
pub(crate) fn determinism_trace(label: &str, brep: &BRep) {
    if std::env::var_os("POPKERN_DETERMINISM_TRACE").is_some() {
        eprintln!("dtrace: {label} {:016x}", determinism_fingerprint(brep));
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
