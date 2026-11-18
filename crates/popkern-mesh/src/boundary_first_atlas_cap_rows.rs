//! The atlas's row bookkeeping between pieces: band row samples carried
//! into the caps, and the merge's closure check along the rows.

use popkern_math::hash::FastHashMap as HashMap;

use popkern_core::geometry::Point3;
use popkern_core::topology::FaceHandle;

use super::super::super::BoundaryRing;
use super::super::Cut;
use crate::mesh::FaceMesh;

/// A row sample: face parameters and position.
pub(super) type RowSample = ((f64, f64), Point3);

/// A boundary edge of `mesh` (used by one triangle) whose ends both lie on
/// a cap row, if any.
pub(super) fn exposed_row_edge(mesh: &FaceMesh, cuts: &[Cut]) -> Option<(usize, usize)> {
    let mut uses: HashMap<(usize, usize), usize> = HashMap::default();
    for t in &mesh.triangles {
        for k in 0..3 {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            *uses.entry((a.min(b), a.max(b))).or_default() += 1;
        }
    }
    let on_row = |i: usize| cuts.iter().any(|c| c.on_row(mesh.vertices[i].uv));
    uses.into_iter()
        .filter(|&(_, n)| n == 1)
        .map(|(e, _)| e)
        .find(|&(a, b)| on_row(a) && on_row(b))
}

/// One face mesh from its pieces' meshes: vertices at one position
/// welded (the row samples two pieces share, the crossings).
pub(super) fn merge(face: FaceHandle, meshes: Vec<FaceMesh>) -> FaceMesh {
    let mut out = FaceMesh {
        face,
        vertices: Vec::new(),
        triangles: Vec::new(),
    };
    let mut at: HashMap<[u64; 3], usize> = HashMap::default();
    for m in meshes {
        let map: Vec<usize> = m
            .vertices
            .iter()
            .map(|v| {
                let key = [
                    v.position.x.to_bits(),
                    v.position.y.to_bits(),
                    v.position.z.to_bits(),
                ];
                *at.entry(key).or_insert_with(|| {
                    out.vertices.push(v.clone());
                    out.vertices.len() - 1
                })
            })
            .collect();
        out.triangles
            .extend(m.triangles.iter().map(|t| t.map(|i| map[i])));
    }
    out
}

/// `rings` with every band vertex lying on `cut`'s row inserted, in order
/// along the row, into the ring segment between two row samples that
/// spans it (the band's slit cuts its row circle at a point the cap's
/// copy never had).
pub(super) fn with_band_row_samples(
    rings: Vec<BoundaryRing>,
    cut: &Cut,
    band: &[RowSample],
) -> Vec<BoundaryRing> {
    let key = |p: &Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    let on_row: Vec<&((f64, f64), Point3)> =
        band.iter().filter(|(uv, _)| cut.on_row(*uv)).collect();
    if on_row.is_empty() {
        return rings;
    }
    rings
        .into_iter()
        .map(|ring| {
            let have: std::collections::HashSet<[u64; 3]> = ring.pos.iter().map(key).collect();
            let missing: Vec<&RowSample> = on_row
                .iter()
                .copied()
                .filter(|(_, p)| !have.contains(&key(p)))
                .collect();
            if missing.is_empty() {
                return ring;
            }
            let n = ring.uv.len();
            let mut out = BoundaryRing {
                uv: Vec::with_capacity(n + missing.len()),
                pos: Vec::with_capacity(n + missing.len()),
                sources: Vec::with_capacity(n + missing.len()),
            };
            for i in 0..n {
                out.uv.push(ring.uv[i]);
                out.pos.push(ring.pos[i]);
                out.sources.push(ring.sources[i]);
                let j = (i + 1) % n;
                if !(cut.on_row(ring.uv[i]) && cut.on_row(ring.uv[j])) {
                    continue;
                }
                // Missing samples between these two row samples in space:
                // on the chord's span, ordered from i.
                let (a, b) = (ring.pos[i], ring.pos[j]);
                let ab = b - a;
                let len2 = ab.dot(&ab);
                if len2 <= 0.0 {
                    continue;
                }
                let mut between: Vec<(f64, &RowSample)> = missing
                    .iter()
                    .map(|m| ((m.1 - a).dot(&ab) / len2, *m))
                    .filter(|(t, m)| {
                        *t > 0.0 && *t < 1.0 && {
                            let foot = a + ab * *t;
                            foot.distance_to(&m.1) <= 0.25 * len2.sqrt()
                        }
                    })
                    .collect();
                between.sort_by(|x, y| x.0.total_cmp(&y.0));
                // Face parameters on any branch serve a cap (it reads
                // positions; parameters are looked up by position).
                for (_, m) in between {
                    out.uv.push(m.0);
                    out.pos.push(m.1);
                    out.sources
                        .push((popkern_core::topology::EdgeHandle::default(), 0.0, 0.0));
                }
            }
            out
        })
        .collect()
}
