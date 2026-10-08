//! A hinge grown on two existing bodies along their seam (PK-702,
//! design `docs/architecture/hinges.md` §5.5): a box and its lid, say.
//!
//! The seam is a straight line where the bodies meet (the box's top back
//! edge). The pin axis runs parallel to it, `outward` of the seam by the
//! knuckle radius plus the clearance, so each body's knuckles clear the
//! other body as it turns. Each body owns alternate knuckle segments; a
//! knuckle is a disc round the axis joined to its body's wall by a web
//! that stays on its own side of the seam — one line–arc–line profile,
//! extruded over the segment, bored, then unioned onto the body. The
//! bodies are only unioned (never cut), so every original face keeps its
//! name.

use popkern_core::geometry::analytic_surfaces::orthonormal_basis;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_math::PortableMath;
use popkern_sketch::{PathSegment, ProfileShape, SketchPlane};

use crate::HingeError;
use crate::construct::{construction, extrude, subtract, union};
use crate::params::{HingeComponent, HingeKind, HingeParams, hinge_layout};

/// Where a seam hinge sits.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Seam {
    /// A point on the seam line where the hinge starts.
    pub origin: Point3,
    /// Along the seam (unit); the hinge runs `length` this way.
    pub axis: Vec3,
    /// Away from both bodies, perpendicular to the seam (unit).
    pub outward: Vec3,
    /// From body A across the seam into body B, perpendicular to both
    /// (unit).
    pub up: Vec3,
}

impl Seam {
    fn check(&self) -> Result<(), HingeError> {
        let unit = |v: Vec3| (v.length() - 1.0).abs() < 1e-9;
        let square = self.axis.dot(&self.outward).abs() < 1e-9
            && self.axis.dot(&self.up).abs() < 1e-9
            && self.outward.dot(&self.up).abs() < 1e-9;
        if !(unit(self.axis) && unit(self.outward) && unit(self.up) && square) {
            return Err(HingeError::InvalidSpec {
                field: "seam",
                value: 0.0,
                requirement: "axis, outward and up are orthonormal",
            });
        }
        Ok(())
    }

    /// The world point at seam coordinates `(n, v, s)`: `n` outward, `v`
    /// up, `s` along the seam.
    fn at(&self, n: f64, v: f64, s: f64) -> Point3 {
        self.origin + self.outward * n + self.up * v + self.axis * s
    }
}

/// The web reaches this fraction of the knuckle radius into its body's
/// wall, and spans from the clearance off the seam to this fraction of
/// the radius below (above) it.
const WEB_DEPTH: f64 = 0.5;
const WEB_SPAN: f64 = 0.75;

/// Build one component of a seam hinge: `LeafA` returns `body` (body A)
/// with its knuckles, `LeafB` body B with its knuckles, `Pin` the pin of
/// a butt hinge (no body). `params.leaf_width` and `params.angle_deg`
/// are unused: the bodies are the leaves.
///
/// # Errors
///
/// [`hinge_layout`]'s declines; [`HingeError::InvalidSpec`] for a seam
/// that is not orthonormal, a missing body, or a clearance too large for
/// the web (≥ half the knuckle radius); construction variants when the
/// kernel declines (a web that misses the wall shows as a union that
/// changes nothing — refused as a construction error).
pub fn seam_hinge_component(
    kind: HingeKind,
    component: HingeComponent,
    seam: &Seam,
    params: &HingeParams,
    body: Option<&BRep>,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    seam.check()?;
    // The leaf width only matters to plate hinges; give the layout one
    // wide enough not to refuse.
    let mut layout_params = params.clone();
    layout_params.leaf_width = layout_params.leaf_width.max(4.0 * params.length + 1.0);
    let layout = hinge_layout(kind, &layout_params)?;
    let r = layout.knuckle_radius;
    let c = params.clearance;
    if c >= 0.5 * r {
        return Err(HingeError::InvalidSpec {
            field: "clearance",
            value: c,
            requirement: "less than half the knuckle radius on a seam hinge",
        });
    }
    let a = r + c.max(params.length * 1e-4);
    if component == HingeComponent::Pin {
        if kind == HingeKind::PrintInPlace {
            return Err(HingeError::InvalidSpec {
                field: "component",
                value: 0.0,
                requirement: "leaf_a or leaf_b: a print-in-place hinge's pin is part of body A",
            });
        }
        return segment_cylinder(
            seam,
            a,
            params.pin_diameter / 2.0,
            0.0,
            params.length,
            tolerance,
        );
    }
    let Some(body) = body else {
        return Err(HingeError::InvalidSpec {
            field: "body",
            value: 0.0,
            requirement: "a body to grow the knuckles on",
        });
    };
    let side = if component == HingeComponent::LeafA {
        -1.0
    } else {
        1.0
    };
    let integral_pin = kind == HingeKind::PrintInPlace && component == HingeComponent::LeafA;
    let margin = 0.5 * r;
    // Each web buries WEB_DEPTH·R into the wall: the wall must be there.
    // Ray-voting classification (the winding number is exact only on
    // straight-edged polyhedra, and a box with rounded corners is not).
    let (near, far) = (c.max(r * 1e-3), WEB_SPAN * r);
    let aabbs: Vec<_> = body
        .faces
        .keys()
        .filter_map(|fh| {
            let pts = body.face_boundary_sample_points(fh, 32);
            popkern_core::geometry::intersection::surface_subdivision::Aabb::from_points(&pts)
                .map(|a| (fh, a))
        })
        .collect();
    for (z0, z1) in layout.owned(component) {
        // A quarter beyond the web end: the wall must not end where the
        // web does (the web end face would be coplanar with its inside).
        let probe = seam.at(
            -WEB_DEPTH * r * 1.25,
            side * 0.5 * (near + far),
            0.5 * (z0 + z1),
        );
        let inside = popkern_bool::point_in_solid_final(&probe, body, &aabbs, &tolerance)
            .is_ok_and(|c| c == popkern_bool::PointClassification::Inside);
        if !inside {
            return Err(HingeError::SeamWall {
                body: if side < 0.0 { "A" } else { "B" },
                depth: WEB_DEPTH * r,
                band: [near, far],
            });
        }
    }
    let mut result = body.clone();
    for (z0, z1) in layout.owned(component) {
        let mut knuckle = extrude_knuckle(seam, side, a, r, c, z0, z1, tolerance)?;
        if !integral_pin {
            let bore = segment_cylinder(
                seam,
                a,
                layout.bore_radius,
                z0 - margin,
                z1 + margin,
                tolerance,
            )?;
            knuckle = subtract(&knuckle, &bore, "seam knuckle bore", tolerance)?;
        }
        let shells = result.shells.len();
        result = union(&result, &knuckle, "seam knuckle", tolerance)?;
        if result.shells.len() != shells {
            return Err(construction(
                "seam knuckle",
                "the knuckle did not join the body (its web misses the wall)".into(),
            ));
        }
    }
    if integral_pin {
        // One pin through every B zone, as in the plate hinge.
        if let Some((lo, hi)) = crate::leaf::integral_pin_span(&layout.segments, params.length) {
            let pin = segment_cylinder(seam, a, params.pin_diameter / 2.0, lo, hi, tolerance)?;
            result = union(&result, &pin, "seam integral pin", tolerance)?;
        }
    }
    Ok(result)
}

/// A sketch plane at seam station `s` spanned by outward and up, its
/// normal along the seam (flipping the station when the frame is
/// left-handed, so the extrusion always runs toward +s).
fn station_plane(seam: &Seam, s0: f64, s1: f64) -> (SketchPlane, f64) {
    let normal = seam.outward.cross(&seam.up);
    let forward = normal.dot(&seam.axis) > 0.0;
    let origin = seam.at(0.0, 0.0, if forward { s0 } else { s1 });
    (
        SketchPlane {
            origin,
            u_axis: seam.outward,
            v_axis: seam.up,
            normal,
        },
        s1 - s0,
    )
}

/// A cylinder of radius `radius` on the pin axis from station `s0` to
/// `s1`.
fn segment_cylinder(
    seam: &Seam,
    a: f64,
    radius: f64,
    s0: f64,
    s1: f64,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    let (plane, depth) = station_plane(seam, s0, s1);
    let shape = ProfileShape::Path {
        start: [a + radius, 0.0],
        segments: arc_pieces(
            seam,
            [a, 0.0],
            radius,
            0.0,
            std::f64::consts::TAU,
            [a + radius, 0.0],
        ),
    };
    extrude(shape, plane, depth, "seam pin", tolerance)
}

/// One knuckle: the disc round the axis at `(a, 0)` and its web into the
/// body on `side` of the seam (−1 below, +1 above), over `[s0, s1]`.
#[allow(clippy::too_many_arguments)]
fn extrude_knuckle(
    seam: &Seam,
    side: f64,
    a: f64,
    r: f64,
    c: f64,
    s0: f64,
    s1: f64,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    let (near, far) = (side * c.max(r * 1e-3), side * WEB_SPAN * r);
    let depth = -WEB_DEPTH * r;
    let n_near = a - (r * r - near * near).sqrt();
    let n_far = a - (r * r - far * far).sqrt();
    // Angles of the web's two meeting points on the disc, and the long
    // way round between them (over the far side of the axis).
    let start = near.patan2(n_near - a);
    let end = far.patan2(n_far - a);
    let (from, to) = if side < 0.0 {
        // Below the seam: from the near point clockwise over the top.
        (
            start.rem_euclid(std::f64::consts::TAU),
            end.rem_euclid(std::f64::consts::TAU) - std::f64::consts::TAU,
        )
    } else {
        // Above: the mirror, counter-clockwise under the bottom.
        let s = start.rem_euclid(std::f64::consts::TAU) - std::f64::consts::TAU;
        (s, end.rem_euclid(std::f64::consts::TAU))
    };
    let mut segments = vec![PathSegment::Line { to: [n_near, near] }];
    segments.extend(arc_pieces(seam, [a, 0.0], r, from, to, [n_far, far]));
    segments.push(PathSegment::Line { to: [depth, far] });
    let shape = ProfileShape::Path {
        start: [depth, near],
        segments,
    };
    let (plane, length) = station_plane(seam, s0, s1);
    extrude(shape, plane, length, "seam knuckle profile", tolerance)
}

/// A circular arc about `centre` (seam-plane coordinates) from angle
/// `from` to `to` (radians, either direction), ending exactly at `end`,
/// split where the extruded cylinder's chart seam and its quarter turns
/// fall, so no piece straddles the seam (PK-704) and none reaches 90°.
fn arc_pieces(
    seam: &Seam,
    centre: [f64; 2],
    r: f64,
    from: f64,
    to: f64,
    end: [f64; 2],
) -> Vec<PathSegment<f64>> {
    use std::f64::consts::{FRAC_PI_2, TAU};
    // The chart's u = 0 direction for this axis, in (outward, up).
    let (perp, _) = orthonormal_basis(seam.outward.cross(&seam.up));
    let psi = perp.dot(&seam.up).patan2(perp.dot(&seam.outward));
    let dir = if to >= from { 1.0 } else { -1.0 };
    let mut cuts = vec![from];
    // Every psi + k·π/2 strictly between from and to.
    let k0 = ((from - psi) / FRAC_PI_2).floor() as i64 - 1;
    let k1 = ((to - psi) / FRAC_PI_2).ceil() as i64 + 1;
    let mut inner: Vec<f64> = (k0.min(k1)..=k0.max(k1))
        .map(|k| psi + k as f64 * FRAC_PI_2)
        .filter(|&t| (t - from) * dir > 1e-9 && (to - t) * dir > 1e-9)
        .collect();
    inner.sort_by(|x, y| ((x - from) * dir).total_cmp(&((y - from) * dir)));
    cuts.extend(inner);
    cuts.push(to);
    debug_assert!((to - from).abs() <= TAU + 1e-9);
    let point = |t: f64| [centre[0] + r * t.pcos(), centre[1] + r * t.psin()];
    cuts.windows(2)
        .enumerate()
        .map(|(k, w)| {
            let last = k + 2 == cuts.len();
            PathSegment::Arc {
                via: point(0.5 * (w[0] + w[1])),
                to: if last { end } else { point(w[1]) },
            }
        })
        .collect()
}
