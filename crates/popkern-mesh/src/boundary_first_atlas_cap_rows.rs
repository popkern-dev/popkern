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

/// The band's rings with every pair sharing a seam traversal joined into
/// one ring without it. The face's ring held both uses of its seam; split
/// at the cap rows, a band straddling the seam came as two rings each
/// holding one use, sharing their corners on the rows, and the
/// triangulation met each corner twice (a rod beside a sphere's poles: its
/// band declined every round). Where ring `a` runs `p0 … pk` and ring `b`
/// runs `pk … p0` (positions to the bit), the seam lies inside the band:
/// the joined ring is `p0`, `b` past the run, `pk`, `a` past the run. The
/// run's ends stay (the rows' points, which the caps share); its interior
/// samples go (the band is meshed across the seam, as each cap is).
pub(super) fn join_seam_rings(mut rings: Vec<BoundaryRing>) -> Vec<BoundaryRing> {
    let bits = |p: &Point3| [p.x.to_bits(), p.y.to_bits(), p.z.to_bits()];
    'join: loop {
        for a in 0..rings.len() {
            for b in (a + 1)..rings.len() {
                let (ra, rb) = (&rings[a], &rings[b]);
                let (na, nb) = (ra.pos.len(), rb.pos.len());
                let at_b: HashMap<[u64; 3], usize> = rb
                    .pos
                    .iter()
                    .enumerate()
                    .map(|(j, p)| (bits(p), j))
                    .collect();
                // A shared run's start in `a`: shared, its successor in `a`
                // is its predecessor in `b`, and its predecessor in `a` is
                // not `b`'s successor (so the run starts here).
                let reversed_step = |i: usize| -> Option<usize> {
                    let &j = at_b.get(&bits(&ra.pos[i]))?;
                    (bits(&ra.pos[(i + 1) % na]) == bits(&rb.pos[(j + nb - 1) % nb])).then_some(j)
                };
                for i in 0..na {
                    if reversed_step(i).is_none() || reversed_step((i + na - 1) % na).is_some() {
                        continue;
                    }
                    // The run p0 = a[i] … pk = a[i + k].
                    let mut k = 1;
                    while k + 1 < na && reversed_step((i + k) % na).is_some() {
                        k += 1;
                    }
                    let (Some(&j_end), Some(&j_start)) = (
                        at_b.get(&bits(&ra.pos[(i + k) % na])),
                        at_b.get(&bits(&ra.pos[i])),
                    ) else {
                        continue;
                    };
                    // `b` past the run: from after p0 (j_start) to before
                    // pk (j_end), cyclic; `a` past the run likewise.
                    let b_len = (j_end + nb - j_start) % nb;
                    let a_len = (i + na - (i + k) % na) % na;
                    let mut order: Vec<(usize, usize)> = vec![(a, i)];
                    order.extend((1..b_len).map(|t| (b, (j_start + t) % nb)));
                    order.push((a, (i + k) % na));
                    order.extend((1..a_len).map(|t| (a, (i + k + t) % na)));
                    // Each kept sample keeps the source of the segment it
                    // starts: the ends start segments into `b` and `a`.
                    let source = |&(r, idx): &(usize, usize), slot: usize| match slot {
                        0 => rings[b].sources[j_start],
                        s if s == b_len => rings[a].sources[(i + k) % na],
                        _ => rings[r].sources[idx],
                    };
                    let joined = BoundaryRing {
                        uv: order.iter().map(|&(r, idx)| rings[r].uv[idx]).collect(),
                        pos: order.iter().map(|&(r, idx)| rings[r].pos[idx]).collect(),
                        sources: order
                            .iter()
                            .enumerate()
                            .map(|(s, o)| source(o, s))
                            .collect(),
                    };
                    if joined.uv.len() < 3 {
                        continue;
                    }
                    rings.remove(b);
                    rings[a] = joined;
                    continue 'join;
                }
            }
        }
        return rings;
    }
}
