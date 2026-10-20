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
use popkern_core::geometry::certify::chord::{rectangle_chord_bound, triangle_chord_bound};
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
    /// A sampled lower bound of that triangle's deviation (the same
    /// correspondence, on a grid): how tight the proven bound is, and
    /// a witness that an `Exceeds` is the mesh's and not the bound's.
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
    let faces = by_face
        .into_iter()
        .map(|(face, tris)| certify_face(brep, mesh, face, &tris, tolerance))
        .collect();
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
    for &t in tris {
        let tri = mesh.triangles[t];
        let p = tri.map(|i| mesh.vertices[i].position);
        let w = unwrap(tri.map(|i| mesh.vertices[i].uv), periods);
        for piece in pole_pieces(surface, w, p, tolerance) {
            let bound = match &piece {
                Piece::Triangle(w, p) => {
                    triangle_chord_bound(&mut source, &breaks, *w, *p, tolerance)
                }
                Piece::Rectangle(u, v, c) => {
                    rectangle_chord_bound(&mut source, &breaks, *u, *v, *c, tolerance)
                }
            };
            match bound {
                Ok(b) if b > out.surface_bound => {
                    out.surface_bound = b;
                    out.worst_triangle = Some(t);
                    out.witness = piece.sampled(surface);
                }
                Ok(_) => {}
                Err(e) => return inconclusive(out, ChordalInconclusive::Decline(e)),
            }
        }
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

#[path = "chordal_boundary.rs"]
mod chordal_boundary;
#[path = "chordal_pieces.rs"]
mod chordal_pieces;
use chordal_boundary::boundary_bound;
pub(crate) use chordal_pieces::{Piece, c0_breaks, pole_pieces, unwrap};

#[cfg(test)]
#[path = "chordal_certificate_tests.rs"]
mod tests;
