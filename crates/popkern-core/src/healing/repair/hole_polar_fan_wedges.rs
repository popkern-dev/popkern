//! Polar-wedge geometry construction for the N-sided G1 hole cap
//! (PK-248). See [`super::hole_polar_fan`] for the serving contract;
//! this module only builds the wedge patches and spoke curves.
//!
//! # Wedge layout
//!
//! Wedge `k` covers rim side `k` (ring vertex `V_k → V_{k+1}`) and
//! the pole `C`. As a degree-5 rectangular Gregory patch (CCW
//! boundary convention of `from_boundary_beziers`):
//!
//! ```text
//!   boundaries[0]  (u = 0 edge)  rim side k, V_k → V_{k+1}
//!   boundaries[1]  (v = 1 edge)  spoke k+1 reversed, V_{k+1} → C
//!   boundaries[2]  (u = 1 edge)  collapsed at C (all CPs = C)
//!   boundaries[3]  (v = 0 edge)  spoke k, C → V_k
//! ```
//!
//! Adjacent wedges share each spoke's quintic control polygon
//! verbatim, so the fan is G0 across spokes exactly.
//!
//! # G1 mechanisms, by seam
//!
//! * **Rim (cap↔wall)**: the single-patch rung's zero-end-leg rim
//!   quintic plus wall-oriented ribbon inner members
//!   ([`super::hole_rim_ribbon`]).
//! * **Spokes (wedge↔wedge)**: the midpoint-offset pairing from the
//!   fillet cap's cross-spoke rule — both wedges' spoke-adjacent
//!   inner members are built from the SAME base point on the spoke
//!   row and the SAME transverse direction `D = n̂ × t̂` with opposite
//!   signs, so every spoke row CP is exactly the midpoint of the two
//!   members and the cross fields are coplanar with the spoke.
//! * **Pole**: spoke tangents at `C` are projected into the ring
//!   tangent plane (vertex-enclosure precondition — straight chords
//!   to a non-planar rim have non-coplanar tangents at `C`, which no
//!   inner-CP choice can repair), and the collapsed side's inner
//!   members come from the Shirman-Séquin enclosure solve.

use crate::geometry::gregory::RectangularGregoryDeg5;
use crate::geometry::gregory::vertex_enclosure::shirman_sequin_vertex_enclosure_solve;
use crate::geometry::{FillBoundaryEdge, Point3, Vec3};
use crate::topology::BRep;

use super::hole_fill::DetectedHole;
use super::hole_polar_fan::PolarFanCapError;
use super::hole_rim_ribbon::{fit_quintic_zero_leg, rim_ribbon_directions};

/// The constructed fan: one wedge per rim side plus the shared spoke
/// control polygons (spoke `k` runs pole → ring vertex `k`).
pub(super) struct PolarWedgeSet {
    /// The pole all wedges collapse to.
    pub pole: Point3,
    /// One degree-5 Gregory wedge per rim side, in ring order.
    pub wedges: Vec<RectangularGregoryDeg5>,
    /// Quintic spoke control polygons, pole → `frames[k].start`.
    pub spokes: Vec<[Point3; 6]>,
}

/// Build the polar wedges for a hole whose ring frames are already
/// oriented (`frames[k].end == frames[k + 1].start`).
pub(super) fn build_polar_wedges(
    brep: &BRep,
    hole: &DetectedHole,
    frames: &[FillBoundaryEdge],
    tolerance: f64,
) -> Result<PolarWedgeSet, PolarFanCapError> {
    let n = frames.len();
    let pole = ring_centroid(frames);
    let ring_normal =
        ring_newell_normal(frames).ok_or_else(|| PolarFanCapError::DegenerateRing {
            detail: "ring polygon has no usable Newell normal".into(),
        })?;

    // Spokes: pole → V_k quintics whose tangent AT the pole lies in
    // the ring tangent plane (Hermite cubic, degree-elevated). On a
    // planar ring this reproduces the straight chord exactly.
    let mut spokes = Vec::with_capacity(n);
    let mut pole_tangents = Vec::with_capacity(n);
    for (k, frame) in frames.iter().enumerate() {
        let chord = frame.start - pole;
        let length = chord.length();
        if length <= tolerance {
            return Err(PolarFanCapError::DegenerateRing {
                detail: format!("ring vertex {k} coincides with the pole ({length:.3e})"),
            });
        }
        let radial = chord * (1.0 / length);
        let in_plane = chord - ring_normal * chord.dot(&ring_normal);
        let Some(pole_tangent) = in_plane.normalized() else {
            return Err(PolarFanCapError::DegenerateRing {
                detail: format!("ring vertex {k} sits on the pole axis; no in-plane spoke tangent"),
            });
        };
        // Zero end leg at the RIM vertex (same corner reality as the
        // rim quintics: at a ring vertex two walls disagree, so the
        // spoke's rim leg must not feed the wall-facing cross field —
        // measured 1.3e-1 mid-band wall G1 with a full Hermite leg,
        // pure corner-leg pollution). The pole leg stays: enclosure
        // coplanarity at the pole needs the in-plane tangent.
        let mut spoke = hermite_quintic(pole, frame.start, pole_tangent * length, radial * length);
        spoke[4] = frame.start;
        spokes.push(spoke);
        pole_tangents.push(pole_tangent);
    }

    // Rim rows and wall ribbons (shared with the single-patch rung).
    let mut rims = Vec::with_capacity(n);
    let mut rim_members = Vec::with_capacity(n);
    for (k, frame) in frames.iter().enumerate() {
        let row = fit_quintic_zero_leg(frame, k, tolerance)?;
        let ribbon = rim_ribbon_directions(brep, hole, frame, k)?;
        let h = frame.start.distance_to(&frame.end) / 3.0;
        rim_members.push([
            midpoint(row[1], row[2]) + ribbon[0] * h,
            midpoint(row[3], row[4]) + ribbon[1] * h,
        ]);
        rims.push(row);
    }

    // Spoke members: shared base + shared transverse direction with
    // opposite signs, so the spoke row CP is the exact midpoint of
    // the two wedges' members (the fillet cross-spoke rule). `plus`
    // is the CCW side of the spoke (wedge k for spoke k), `minus`
    // the CW side (wedge k−1).
    let mut spoke_members = Vec::with_capacity(n);
    for (k, spoke) in spokes.iter().enumerate() {
        let h = pole.distance_to(&frames[k].start) / 3.0;
        let station = |a: Point3, b: Point3, fallback: Vec3| -> Result<_, PolarFanCapError> {
            let tangent = (b - a).normalized().unwrap_or(fallback);
            let transverse = ring_normal.cross(&tangent).normalized().ok_or_else(|| {
                PolarFanCapError::DegenerateRing {
                    detail: format!("spoke {k} runs along the ring normal; no transverse field"),
                }
            })?;
            let base = midpoint(a, b);
            Ok((base + transverse * h, base - transverse * h))
        };
        let near_pole = station(spoke[1], spoke[2], pole_tangents[k])?;
        let near_rim = station(spoke[3], spoke[4], pole_tangents[k])?;
        spoke_members.push(SpokeMembers {
            near_pole,
            near_rim,
        });
    }

    // Pole enclosure: Shirman-Séquin solve over the parallelogram
    // raws; `spoke_directions[k]` is the spoke shared between wedges
    // k and k+1, i.e. spoke k+1.
    let raws: Vec<Point3> = (0..n)
        .map(|k| pole + (spokes[k][1] - pole) + (spokes[(k + 1) % n][1] - pole))
        .collect();
    let dirs: Vec<Vec3> = (0..n).map(|k| pole_tangents[(k + 1) % n]).collect();
    let enclosure = shirman_sequin_vertex_enclosure_solve(pole, ring_normal, &raws, &dirs)
        .ok_or(PolarFanCapError::EnclosureUnsolvable)?;

    // Assemble the wedges.
    let mut wedges = Vec::with_capacity(n);
    for k in 0..n {
        let next = (k + 1) % n;
        let boundaries = [rims[k], reversed(&spokes[next]), [pole; 6], spokes[k]];
        let inner_cps = [
            (rim_members[k][0], spoke_members[k].near_rim.0),
            (rim_members[k][1], spoke_members[next].near_rim.1),
            (enclosure[k], spoke_members[next].near_pole.1),
            (enclosure[k], spoke_members[k].near_pole.0),
        ];
        let wedge = RectangularGregoryDeg5::from_boundary_beziers(
            boundaries,
            inner_cps,
            tolerance.max(1e-9),
        )
        .map_err(|error| PolarFanCapError::Assembly {
            wedge: k,
            detail: error.to_string(),
        })?;
        wedges.push(wedge);
    }

    Ok(PolarWedgeSet {
        pole,
        wedges,
        spokes,
    })
}

/// Members flanking one spoke at its two interior stations, as
/// `(ccw_side, cw_side)` pairs.
struct SpokeMembers {
    /// Members at the station nearer the pole (t ≈ 0.3).
    near_pole: (Point3, Point3),
    /// Members at the station nearer the rim vertex (t ≈ 0.7).
    near_rim: (Point3, Point3),
}

/// Centroid of the ring vertices (`frames[k].start` covers each
/// vertex exactly once on an oriented ring).
fn ring_centroid(frames: &[FillBoundaryEdge]) -> Point3 {
    let inverse = 1.0 / frames.len() as f64;
    let mut c = Point3::new(0.0, 0.0, 0.0);
    for frame in frames {
        c = Point3::new(
            c.x + frame.start.x * inverse,
            c.y + frame.start.y * inverse,
            c.z + frame.start.z * inverse,
        );
    }
    c
}

/// Unit Newell normal of the ring polygon sampled at each side's
/// start and midpoint (mids matter when the rim curves).
fn ring_newell_normal(frames: &[FillBoundaryEdge]) -> Option<Vec3> {
    let points: Vec<Point3> = frames
        .iter()
        .flat_map(|frame| [frame.start, frame.mid])
        .collect();
    let mut normal = Vec3::new(0.0, 0.0, 0.0);
    for (index, current) in points.iter().enumerate() {
        let next = &points[(index + 1) % points.len()];
        normal.x += (current.y - next.y) * (current.z + next.z);
        normal.y += (current.z - next.z) * (current.x + next.x);
        normal.z += (current.x - next.x) * (current.y + next.y);
    }
    normal.normalized()
}

/// Quintic Bezier interpolating `start → end` with the given end
/// derivatives (cubic Hermite, degree-elevated twice — exact).
fn hermite_quintic(start: Point3, end: Point3, d_start: Vec3, d_end: Vec3) -> [Point3; 6] {
    let cubic = [
        start,
        start + d_start * (1.0 / 3.0),
        end - d_end * (1.0 / 3.0),
        end,
    ];
    let mut quartic = [Point3::new(0.0, 0.0, 0.0); 5];
    for (i, slot) in quartic.iter_mut().enumerate() {
        let w = i as f64 / 4.0;
        let lower = cubic[i.saturating_sub(1)];
        let same = cubic[i.min(3)];
        *slot = lerp(same, lower, w);
    }
    let mut quintic = [Point3::new(0.0, 0.0, 0.0); 6];
    for (i, slot) in quintic.iter_mut().enumerate() {
        let w = i as f64 / 5.0;
        let lower = quartic[i.saturating_sub(1)];
        let same = quartic[i.min(4)];
        *slot = lerp(same, lower, w);
    }
    quintic
}

/// `a + w · (b − a)`.
fn lerp(a: Point3, b: Point3, w: f64) -> Point3 {
    Point3::new(
        a.x + w * (b.x - a.x),
        a.y + w * (b.y - a.y),
        a.z + w * (b.z - a.z),
    )
}

/// Midpoint of two points.
fn midpoint(a: Point3, b: Point3) -> Point3 {
    Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z))
}

/// The control polygon walked in the opposite direction.
fn reversed(cps: &[Point3; 6]) -> [Point3; 6] {
    let mut out = *cps;
    out.reverse();
    out
}
