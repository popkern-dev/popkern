//! Living hinges (design §5.4): two panels joined by a flexure, each
//! style ONE profile extruded once — no booleans.
//!
//! Frame: the bend axis runs along +Z from `z = 0` to `z = length`;
//! panel A lies toward −X, panel B toward +X; the mounting (bottom) face
//! is the plane `y = 0` and the panels' top faces `y = thickness`.
//!
//! - **Straight:** a thin web of thickness `web` on the mounting side
//!   spans the flex zone, its roots rounded with radius `root_radius`
//!   (tangent arcs in the profile). For flexible materials (PP, PETG,
//!   TPU).
//! - **Slotted:** a serpentine kerf pattern for rigid materials — slots
//!   across the flex zone, parallel to the bend axis, open alternately to
//!   the `z = 0` and `z = length` edges, each closed end a semicircle.
//!   The panel stays one simple outline, so the whole hinge is one
//!   extrusion through the thickness.

use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{PathSegment, ProfileShape};

use crate::HingeError;
use crate::construct::{extrude, plane};

/// Which flexure.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LivingHingeStyle {
    /// A thin web with rounded roots.
    Straight,
    /// A serpentine of slots open alternately to each edge.
    Slotted,
}

impl LivingHingeStyle {
    /// The skin spelling (`straight`, `slotted`).
    pub fn spelling(self) -> &'static str {
        match self {
            LivingHingeStyle::Straight => "straight",
            LivingHingeStyle::Slotted => "slotted",
        }
    }

    /// The style a skin spelling names.
    pub fn parse(text: &str) -> Option<LivingHingeStyle> {
        match text {
            "straight" => Some(LivingHingeStyle::Straight),
            "slotted" => Some(LivingHingeStyle::Slotted),
            _ => None,
        }
    }
}

/// A living hinge's dimensions, in model units.
#[derive(Debug, Clone, PartialEq)]
pub struct LivingHingeParams {
    /// The flexure.
    pub style: LivingHingeStyle,
    /// Along the bend axis.
    pub length: f64,
    /// Width of the flex zone (across the bend axis).
    pub span: f64,
    /// Each panel's width beyond the flex zone.
    pub panel_width: f64,
    /// Panel thickness.
    pub thickness: f64,
    /// Straight: web thickness.
    pub web: f64,
    /// Straight: root fillet radius (0 = sharp roots).
    pub root_radius: f64,
    /// Slotted: distance between slot centre lines.
    pub slot_pitch: f64,
    /// Slotted: slot width.
    pub slot_width: f64,
    /// Slotted: material left between a slot's closed end and the edge.
    pub bridge: f64,
}

impl LivingHingeParams {
    /// A straight living hinge with a flex zone of `span` and length
    /// `length`; panels 2·span wide, thickness span/5, web a fifth of
    /// the thickness, roots rounded at half the room above the web.
    pub fn straight(span: f64, length: f64) -> LivingHingeParams {
        let thickness = span / 5.0;
        let web = thickness / 5.0;
        LivingHingeParams {
            style: LivingHingeStyle::Straight,
            length,
            span,
            panel_width: 2.0 * span,
            thickness,
            web,
            root_radius: ((thickness - web) / 2.0).min(span / 4.0),
            slot_pitch: span / 8.0,
            slot_width: span / 24.0,
            bridge: length / 10.0,
        }
    }

    /// A slotted (serpentine) living hinge: as [`LivingHingeParams::straight`]
    /// with eight slots across the span, a third of the pitch wide, and
    /// bridges of a tenth of the length.
    pub fn slotted(span: f64, length: f64) -> LivingHingeParams {
        LivingHingeParams {
            style: LivingHingeStyle::Slotted,
            ..LivingHingeParams::straight(span, length)
        }
    }

    /// Slot centre lines across the flex zone (slotted style).
    pub fn slot_centres(&self) -> Vec<f64> {
        let count = (self.span / self.slot_pitch).floor().max(0.0) as usize;
        let first = -0.5 * self.slot_pitch * (count as f64 - 1.0);
        (0..count)
            .map(|k| first + self.slot_pitch * k as f64)
            .collect()
    }
}

fn positive(field: &'static str, value: f64) -> Result<(), HingeError> {
    if value.is_finite() && value > 0.0 {
        Ok(())
    } else {
        Err(HingeError::InvalidSpec {
            field,
            value,
            requirement: "finite and > 0",
        })
    }
}

fn fail(field: &'static str, value: f64, requirement: &'static str) -> HingeError {
    HingeError::InvalidSpec {
        field,
        value,
        requirement,
    }
}

/// Check `p` for its style.
///
/// # Errors
///
/// [`HingeError::InvalidSpec`] naming the field and its limit.
pub fn check_living_hinge(p: &LivingHingeParams) -> Result<(), HingeError> {
    positive("length", p.length)?;
    positive("span", p.span)?;
    positive("panel_width", p.panel_width)?;
    positive("thickness", p.thickness)?;
    match p.style {
        LivingHingeStyle::Straight => {
            positive("web", p.web)?;
            if !(p.root_radius.is_finite() && p.root_radius >= 0.0) {
                return Err(fail("root_radius", p.root_radius, "finite and ≥ 0"));
            }
            if p.web + p.root_radius >= p.thickness {
                return Err(fail(
                    "root_radius",
                    p.root_radius,
                    "less than the thickness above the web",
                ));
            }
            if 2.0 * p.root_radius >= p.span {
                return Err(fail(
                    "root_radius",
                    p.root_radius,
                    "less than half the span",
                ));
            }
        }
        LivingHingeStyle::Slotted => {
            positive("slot_pitch", p.slot_pitch)?;
            positive("slot_width", p.slot_width)?;
            positive("bridge", p.bridge)?;
            if p.slot_width >= p.slot_pitch {
                return Err(fail(
                    "slot_width",
                    p.slot_width,
                    "less than the slot pitch (ribs need material)",
                ));
            }
            if p.slot_pitch > p.span {
                return Err(fail("slot_pitch", p.slot_pitch, "at most the span"));
            }
            if p.bridge + p.slot_width / 2.0 >= p.length {
                return Err(fail(
                    "bridge",
                    p.bridge,
                    "short enough to leave a slot (bridge + slot_width/2 < length)",
                ));
            }
        }
    }
    Ok(())
}

/// Build a living hinge: one body, both panels and the flexure.
///
/// # Errors
///
/// [`check_living_hinge`]'s declines; construction variants when the
/// kernel declines the extrusion.
pub fn living_hinge(p: &LivingHingeParams, tolerance: Tolerance) -> Result<BRep, HingeError> {
    check_living_hinge(p)?;
    match p.style {
        LivingHingeStyle::Straight => extrude(
            straight_profile(p),
            plane(
                Point3::new(0.0, 0.0, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 1.0, 0.0),
            ),
            p.length,
            "straight living hinge",
            tolerance,
        ),
        // The outline lies in the XZ plane; the plane's normal X × Z is
        // −Y, so start on the top face and extrude down through it.
        LivingHingeStyle::Slotted => extrude(
            slotted_profile(p),
            plane(
                Point3::new(0.0, p.thickness, 0.0),
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 1.0),
            ),
            p.thickness,
            "slotted living hinge",
            tolerance,
        ),
    }
}

/// A quarter arc about `centre` from angle `a` to `b` (radians).
fn quarter(centre: [f64; 2], r: f64, a: f64, b: f64) -> PathSegment<f64> {
    use popkern_math::PortableMath;
    let m = 0.5 * (a + b);
    PathSegment::Arc {
        via: [centre[0] + r * m.pcos(), centre[1] + r * m.psin()],
        to: [centre[0] + r * b.pcos(), centre[1] + r * b.psin()],
    }
}

/// The straight hinge's cross-section in (x, y).
fn straight_profile(p: &LivingHingeParams) -> ProfileShape<f64> {
    use std::f64::consts::{FRAC_PI_2, PI};
    let (h, x_out, t, w, r) = (
        p.span / 2.0,
        p.span / 2.0 + p.panel_width,
        p.thickness,
        p.web,
        p.root_radius,
    );
    let line = |x: f64, y: f64| PathSegment::Line { to: [x, y] };
    let mut segments = vec![line(x_out, 0.0), line(x_out, t), line(h, t)];
    if r > 0.0 {
        // Down panel B's inner face, round the concave root (centre
        // inside the gap above the web), along the web, round the other
        // root, up panel A's inner face.
        segments.push(line(h, w + r));
        segments.push(quarter([h - r, w + r], r, 0.0, -FRAC_PI_2));
        segments.push(line(-h + r, w));
        segments.push(quarter([-h + r, w + r], r, -FRAC_PI_2, -PI));
    } else {
        segments.push(line(h, w));
        segments.push(line(-h, w));
    }
    segments.push(line(-h, t));
    segments.push(line(-x_out, t));
    ProfileShape::Path {
        start: [-x_out, 0.0],
        segments,
    }
}

/// The slotted hinge's outline in (x, z): along the `z = 0` edge left
/// to right, dipping into each slot open there; back along `z = length`
/// right to left, dipping into each slot open there.
fn slotted_profile(p: &LivingHingeParams) -> ProfileShape<f64> {
    use std::f64::consts::{FRAC_PI_2, PI};
    let (x_out, l, half, b) = (
        p.span / 2.0 + p.panel_width,
        p.length,
        p.slot_width / 2.0,
        p.bridge,
    );
    let centres = p.slot_centres();
    let line = |x: f64, z: f64| PathSegment::Line { to: [x, z] };
    let mut segments = Vec::new();
    // Slot k is open at z = 0 when k is even, at z = length when odd.
    for (k, &c) in centres.iter().enumerate() {
        if k % 2 == 0 {
            let top = l - b - half;
            segments.push(line(c - half, 0.0));
            segments.push(line(c - half, top));
            segments.push(quarter([c, top], half, PI, FRAC_PI_2));
            segments.push(quarter([c, top], half, FRAC_PI_2, 0.0));
            segments.push(line(c + half, 0.0));
        }
    }
    segments.push(line(x_out, 0.0));
    segments.push(line(x_out, l));
    for (k, &c) in centres.iter().enumerate().rev() {
        if k % 2 == 1 {
            let bottom = b + half;
            segments.push(line(c + half, l));
            segments.push(line(c + half, bottom));
            segments.push(quarter([c, bottom], half, 0.0, -FRAC_PI_2));
            segments.push(quarter([c, bottom], half, -FRAC_PI_2, -PI));
            segments.push(line(c - half, l));
        }
    }
    segments.push(line(-x_out, l));
    ProfileShape::Path {
        start: [-x_out, 0.0],
        segments,
    }
}
