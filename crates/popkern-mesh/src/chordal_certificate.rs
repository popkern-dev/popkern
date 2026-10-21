//! Certified chordal deviation of an emitted mesh (PK-505, design
//! `docs/architecture/certified-chordal-bound.md`).
//!
//! Per face: every triangle is bounded against the surface patch over
//! its parameter triangle (`triangle_chord_bound`), and every boundary
//! segment against the edge curve it lies on (`curve_piece_bound` plus
//! the endpoints' distance from the curve). A face is `Certified` when
//! both stay within the tolerance, `Exceeds` when a bound is proven but
//! larger, and `Inconclusive` when no bound could be proven — never
//! accepted on sampled evidence.

use std::collections::HashMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::certify::CertifyDecline;
use popkern_core::geometry::certify::chord::{
    bent_triangle_chord_bound, quad_pair_offset, rectangle_chord_bound, triangle_chord_bound,
};
use popkern_core::geometry::certify::surface_jet::SurfaceJetSource;
use popkern_core::topology::{BRep, FaceHandle};

use crate::mesh::TriMesh;

/// Why a face's bound could not be proven.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChordalInconclusive {
    /// The face carries no surface geometry.
    NoSurface,
    /// The surface or an edge curve has no certified enclosure.
    Decline(CertifyDecline),
    /// A boundary segment has no edge curve of the face to be bounded
    /// against (the face's edges are all curveless).
    SegmentOnNoEdge {
        /// The segment's first end.
        at: Point3,
    },
}

/// A face's verdict.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ChordalStatus {
    /// Every triangle and boundary segment is proven within the
    /// tolerance.
    Certified,
    /// Proven bounds, at least one above the tolerance.
    Exceeds,
    /// No bound could be proven.
    Inconclusive(ChordalInconclusive),
}

/// One face's certificate.
#[derive(Debug, Clone, PartialEq)]
pub struct FaceChordal {
    /// The face.
    pub face: FaceHandle,
    /// Triangles of this face in the mesh.
    pub triangles: usize,
    /// The largest proven triangle bound.
    pub surface_bound: f64,
    /// Index (into `TriMesh::triangles`) of the triangle with it.
    pub worst_triangle: Option<usize>,
    /// A sampled lower bound of the face's chordal deviation: the
    /// largest distance from the surface, on a grid over that
    /// triangle's parameter piece, to the face's nearest emitted
    /// triangle. It shows how tight the proven bound is, and whether an
    /// `Exceeds` is the mesh's (witness past the tolerance too) or the
    /// bound's.
    pub witness: f64,
    /// The largest proven boundary-segment bound.
    pub boundary_bound: f64,
    /// The verdict.
    pub status: ChordalStatus,
}

/// The certificate of a whole mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct ChordalCertificate {
    /// The tolerance certified against.
    pub tolerance: f64,
    /// One entry per face with triangles, in mesh order.
    pub faces: Vec<FaceChordal>,
}

impl ChordalCertificate {
    /// Whether every face is certified within the tolerance.
    pub fn certified(&self) -> bool {
        self.faces
            .iter()
            .all(|f| f.status == ChordalStatus::Certified)
    }

    /// The largest proven bound over every face (surface or boundary).
    pub fn worst(&self) -> f64 {
        self.faces
            .iter()
            .map(|f| f.surface_bound.max(f.boundary_bound))
            .fold(0.0, f64::max)
    }

    /// Faces whose bound could not be proven.
    pub fn inconclusive(&self) -> Vec<&FaceChordal> {
        self.faces
            .iter()
            .filter(|f| matches!(f.status, ChordalStatus::Inconclusive(_)))
            .collect()
    }

    /// Faces proven above the tolerance.
    pub fn exceeding(&self) -> Vec<&FaceChordal> {
        self.faces
            .iter()
            .filter(|f| f.status == ChordalStatus::Exceeds)
            .collect()
    }
}

/// Certify `mesh` (tessellated from `brep`) against `tolerance`.
pub fn chordal_certificate(brep: &BRep, mesh: &TriMesh, tolerance: f64) -> ChordalCertificate {
    let mut by_face: Vec<(FaceHandle, Vec<usize>)> = Vec::new();
    let mut slot: HashMap<FaceHandle, usize> = HashMap::new();
    for (t, &face) in mesh.face_map.iter().enumerate() {
        let k = *slot.entry(face).or_insert_with(|| {
            by_face.push((face, Vec::new()));
            by_face.len() - 1
        });
        by_face[k].1.push(t);
    }
    // Faces are independent (each prepares its own jet source): they
    // run like tessellation's face pass, in parallel and collected in
    // order, so the certificate is the same at any thread count.
    let faces = crate::face_pass::dispatch(by_face, false, |(face, tris)| {
        certify_face(brep, mesh, face, &tris, tolerance)
    });
    ChordalCertificate { tolerance, faces }
}

fn certify_face(
    brep: &BRep,
    mesh: &TriMesh,
    face: FaceHandle,
    tris: &[usize],
    tolerance: f64,
) -> FaceChordal {
    let mut out = FaceChordal {
        face,
        triangles: tris.len(),
        surface_bound: 0.0,
        worst_triangle: None,
        witness: 0.0,
        boundary_bound: 0.0,
        status: ChordalStatus::Certified,
    };
    let inconclusive = |mut out: FaceChordal, why| {
        out.status = ChordalStatus::Inconclusive(why);
        out
    };
    let Some(surface) = brep
        .faces
        .get(face)
        .and_then(|f| f.surface.as_ref())
        .map(|s| &s.geometry)
    else {
        return inconclusive(out, ChordalInconclusive::NoSurface);
    };
    let mut source = match SurfaceJetSource::new(surface) {
        Ok(s) => s,
        Err(e) => return inconclusive(out, ChordalInconclusive::Decline(e)),
    };
    let breaks = c0_breaks(surface);
    let periods = surface.periodicity();
    let mut bound_piece = |piece: &Piece| match piece {
        Piece::Triangle(w, p) => triangle_chord_bound(&mut source, &breaks, *w, *p, tolerance),
        Piece::Rectangle(u, v, c) => {
            rectangle_chord_bound(&mut source, &breaks, *u, *v, *c, tolerance)
        }
        Piece::Quad(u, v, c, main) => match quad_pair_offset(*c, *main) {
            Some(h) => {
                rectangle_chord_bound(&mut source, &breaks, *u, *v, *c, tolerance).map(|b| b + h)
            }
            None => Ok(f64::INFINITY),
        },
    };
    // Each triangle's bound and the piece holding it; `plain` keeps the
    // unwrapped parameters of triangles bounded as one affine piece.
    let mut best: Vec<(f64, Piece)> = Vec::with_capacity(tris.len());
    let mut plain: Vec<Option<[(f64, f64); 3]>> = Vec::with_capacity(tris.len());
    for &t in tris {
        let tri = mesh.triangles[t];
        let p = tri.map(|i| mesh.vertices[i].position);
        let w = unwrap(tri.map(|i| mesh.vertices[i].uv), periods);
        let pieces = pole_pieces(surface, w, p, tolerance);
        plain.push(matches!(pieces[..], [Piece::Triangle(..)]).then_some(w));
        let mut worst = (0.0_f64, Piece::Triangle(w, p));
        for piece in pieces {
            match bound_piece(&piece) {
                Ok(b) if b >= worst.0 => worst = (b, piece),
                Ok(_) => {}
                Err(e) => return inconclusive(out, ChordalInconclusive::Decline(e)),
            }
        }
        best.push(worst);
    }
    // An affine correspondence charges a trapezoid's in-plane shear (a
    // cone strip between two rulings: 2.6e-2 parametric, 2e-3 apart).
    // Two triangles tiling a parameter rectangle are bounded together,
    // through the bilinear patch, and the smaller bound is kept.
    if best.iter().any(|(b, _)| *b > tolerance) {
        let mut by_edge: HashMap<(usize, usize), Vec<usize>> = HashMap::new();
        for (k, &t) in tris.iter().enumerate() {
            let tri = mesh.triangles[t];
            for e in 0..3 {
                let (a, b) = (tri[e], tri[(e + 1) % 3]);
                by_edge.entry((a.min(b), a.max(b))).or_default().push(k);
            }
        }
        for k in 0..tris.len() {
            let Some(w) = plain[k].filter(|_| best[k].0 > tolerance) else {
                continue;
            };
            let vt = mesh.triangles[tris[k]];
            let pt = vt.map(|i| mesh.vertices[i].position);
            for e in 0..3 {
                let (a, b) = (vt[e], vt[(e + 1) % 3]);
                let Some(partners) = by_edge.get(&(a.min(b), a.max(b))) else {
                    continue;
                };
                for &j in partners.iter().filter(|&&j| j != k && plain[j].is_some()) {
                    let vs = mesh.triangles[tris[j]];
                    let Some(&o) = vs.iter().find(|&&i| i != a && i != b) else {
                        continue;
                    };
                    let off = (mesh.vertices[o].uv, mesh.vertices[o].position);
                    let Some(quad) = quad_piece((vt, w, pt), vs, off, periods) else {
                        continue;
                    };
                    let b = match bound_piece(&quad) {
                        Ok(b) => b,
                        Err(e) => return inconclusive(out, ChordalInconclusive::Decline(e)),
                    };
                    // Either quad tiling a triangle bounds it.
                    for i in [k, j] {
                        if b < best[i].0 {
                            best[i] = (b, quad);
                        }
                    }
                }
            }
        }
    }
    // Triangles still over: pair them through bent edges, which do not
    // charge a surface's sag inside its tangent plane (PK-908).
    if best.iter().any(|(b, _)| *b > tolerance) {
        let mut weld: HashMap<[u64; 3], usize> = HashMap::new();
        let mut id = |p: Point3| {
            let next = weld.len();
            *weld
                .entry([p.x, p.y, p.z].map(|x| (x + 0.0).to_bits()))
                .or_insert(next)
        };
        let input: Vec<Option<BendTriangle>> = tris
            .iter()
            .zip(&plain)
            .map(|(&t, w)| {
                let p = mesh.triangles[t].map(|i| mesh.vertices[i].position);
                w.map(|w| BendTriangle {
                    weld: p.map(&mut id),
                    w,
                    p,
                })
            })
            .collect();
        let mut bounds: Vec<f64> = best.iter().map(|b| b.0).collect();
        bend_pass(surface, &input, &mut bounds, tolerance, |w, p, c| {
            bent_triangle_chord_bound(&mut source, &breaks, w, p, c, tolerance)
                .ok()
                .flatten()
        });
        for (b, x) in best.iter_mut().zip(bounds) {
            b.0 = x;
        }
    }
    if let Some((k, (b, piece))) = best
        .iter()
        .enumerate()
        .max_by(|x, y| x.1.0.total_cmp(&y.1.0))
    {
        out.surface_bound = *b;
        out.worst_triangle = Some(tris[k]);
        let all: Vec<[Point3; 3]> = tris
            .iter()
            .map(|&t| mesh.triangles[t].map(|i| mesh.vertices[i].position))
            .collect();
        out.witness = piece.sampled_against(surface, &all);
    }
    match boundary_bound(brep, mesh, face, tris, tolerance) {
        Ok(b) => out.boundary_bound = b,
        Err(why) => return inconclusive(out, why),
    }
    if out.surface_bound > tolerance || out.boundary_bound > tolerance {
        out.status = ChordalStatus::Exceeds;
    }
    out
}

#[path = "chordal_bend.rs"]
mod chordal_bend;
#[path = "chordal_boundary.rs"]
mod chordal_boundary;
#[path = "chordal_pieces.rs"]
mod chordal_pieces;
use chordal_bend::{BendTriangle, bend_pass};
use chordal_boundary::boundary_bound;
pub(crate) use chordal_pieces::{Piece, c0_breaks, pole_pieces, quad_piece, unwrap};

#[cfg(test)]
#[path = "chordal_certificate_tests.rs"]
mod tests;
