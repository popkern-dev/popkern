//! The period rule (plan §4.3a), measured along the face: two
//! candidates on one face are one point only when the face itself
//! joins them within their zones.
//!
//! Each candidate knows where it lies on the face, in the face's chart
//! and period branch (`Incidence::Face::uv`). Their separation on the
//! face is the length of the surface path over the uv segment between
//! them, after moving one by whole periods to the nearest branch: two
//! points either side of a seam are one period apart in uv and none on
//! the face, and two points beside a pole are far apart in u and
//! nowhere on the surface; but two points of a face that folds back on
//! itself (adjacent turns of a coil touching) are close in space and a
//! whole flank apart on the face, like two points on an edge curling
//! back near itself.
//!
//! A uv is evidence only when the surface maps it back onto the
//! candidate within its radius: a record's own uv curve misses its 3D
//! ends by 3e-3..2e-2 on cones and spheres (PK-1031), and trusting it
//! separated 36 same-point pairs in 11 booleans. Measured with that
//! check (bool, examples, mesh, blend; table on): 2786/2786 pass and
//! the rule separates nothing outside its seam, pole and fold fixtures.

use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::topology::BRep;

use super::candidates::{CandidateVertex, Incidence};

/// Samples along the uv segment for its surface length.
const PATH_SAMPLES: usize = 32;

/// `d` moved by whole periods `period` to the nearest branch.
fn nearest_branch(d: f64, period: Option<f64>) -> f64 {
    match period {
        Some(p) if p > 0.0 => d - p * (d / p).round(),
        _ => d,
    }
}

/// The surface length between `p` and `q` along a face both lie on,
/// when it exceeds their zones together: they are a whole stretch of
/// the face apart however close in space.
pub(crate) fn separation_along_face(
    operands: [&BRep; 2],
    p: &CandidateVertex,
    q: &CandidateVertex,
) -> Option<f64> {
    for ip in &p.incidence {
        let &Incidence::Face {
            operand,
            face,
            uv: Some(a),
        } = ip
        else {
            continue;
        };
        let Some(b) = q.incidence.iter().find_map(|iq| match *iq {
            Incidence::Face {
                operand: o,
                face: f,
                uv: Some(b),
            } if o == operand && f == face => Some(b),
            _ => None,
        }) else {
            continue;
        };
        let Some(surface) = operands[operand]
            .faces
            .get(face)
            .and_then(|f| f.surface.as_ref())
        else {
            continue;
        };
        let (pu, pv) = surface.geometry.periodicity();
        let (du, dv) = (nearest_branch(b.0 - a.0, pu), nearest_branch(b.1 - a.1, pv));
        let at = |k: usize| {
            let w = k as f64 / PATH_SAMPLES as f64;
            surface.geometry.point_at(a.0 + du * w, a.1 + dv * w)
        };
        let length: f64 = (1..=PATH_SAMPLES)
            .map(|k| at(k - 1).distance_to(&at(k)))
            .sum();
        if length > p.zone + q.zone {
            return Some(length);
        }
    }
    None
}
