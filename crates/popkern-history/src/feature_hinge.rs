//! Parameter types for the hinge features (PK-693, design
//! `docs/architecture/hinges.md`):
//! [`Feature::Hinge`](crate::feature::Feature::Hinge) components and
//! [`Feature::LivingHinge`](crate::feature::Feature::LivingHinge) bodies.
//!
//! The tree builds one body per feature, so a hinge records one feature
//! per component (leaf A, leaf B, the pin). They share one
//! [`HingeDesc`]: its parameters are owned by the first component
//! recorded (its `owner` flag; only its parameter list names them), so
//! a hinge shows one slider set and an edit replays every component.

use serde::{Deserialize, Serialize};

use crate::feature_fastener::ThreadDesc;
use crate::parameter::ParameterId;

/// Which hinge (mirrors `popkern_hinge::HingeKind`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum HingeKindDesc {
    /// Two leaves and a separate pin.
    Butt,
    /// Two leaves printed together, the pin part of leaf A.
    PrintInPlace,
}

/// Which body of a hinge a feature builds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum HingeComponentDesc {
    /// The leaf toward −X (owns the hinge's parameters).
    LeafA,
    /// The leaf toward +X (turned to the open angle).
    LeafB,
    /// The separate pin of a butt hinge.
    Pin,
}

impl HingeComponentDesc {
    /// The skin spelling (`leaf_a`, `leaf_b`, `pin`).
    pub fn spelling(self) -> &'static str {
        match self {
            HingeComponentDesc::LeafA => "leaf_a",
            HingeComponentDesc::LeafB => "leaf_b",
            HingeComponentDesc::Pin => "pin",
        }
    }
}

/// Screw holes through each leaf.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct HingeHolesDesc {
    /// Holes per leaf (rounded to a whole count at replay).
    pub per_leaf: ParameterId,
    /// The screw the holes clear.
    pub thread: ThreadDesc,
    /// Sink a 90° countersink.
    pub countersink: bool,
}

/// A hinge's dimensions, shared by its components.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct HingeDesc {
    /// Along the pin axis.
    pub length: ParameterId,
    /// Each leaf, from the pin axis to its outer edge.
    pub leaf_width: ParameterId,
    /// Leaf plate thickness.
    pub thickness: ParameterId,
    /// Total knuckle count (rounded to a whole count at replay, like
    /// pattern counts).
    pub knuckles: ParameterId,
    /// Pin diameter.
    pub pin_diameter: ParameterId,
    /// Knuckle outer diameter (`None`: a leaf rolled round the pin).
    #[serde(default)]
    pub knuckle_diameter: Option<ParameterId>,
    /// Radial clearance (bore over pin, notch over knuckle).
    pub clearance: ParameterId,
    /// Axial gap between knuckles.
    pub gap: ParameterId,
    /// Open angle (an angle parameter: 180° flat).
    pub angle: ParameterId,
    /// Screw holes, if any.
    #[serde(default)]
    pub holes: Option<HingeHolesDesc>,
}

impl HingeDesc {
    /// Every parameter with its role, in a fixed order.
    pub fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        let mut v = vec![
            ("length", self.length),
            ("leaf_width", self.leaf_width),
            ("thickness", self.thickness),
            ("knuckles", self.knuckles),
            ("pin_diameter", self.pin_diameter),
            ("clearance", self.clearance),
            ("gap", self.gap),
            ("angle", self.angle),
        ];
        if let Some(id) = self.knuckle_diameter {
            v.push(("knuckle_diameter", id));
        }
        if let Some(h) = &self.holes {
            v.push(("holes", h.per_leaf));
            v.push(("screw_major_diameter", h.thread.major_diameter));
            v.push(("screw_pitch", h.thread.pitch));
        }
        v
    }
}

/// Which flexure (mirrors `popkern_hinge::LivingHingeStyle`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum LivingHingeStyleDesc {
    /// A thin web with rounded roots.
    Straight,
    /// A serpentine of slots open alternately to each edge.
    Slotted,
}

/// A living hinge's dimensions; the style's own parameters are set,
/// the other style's are `None` (so no unused slider shows).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct LivingHingeDesc {
    /// Along the bend axis.
    pub length: ParameterId,
    /// Width of the flex zone.
    pub span: ParameterId,
    /// Each panel's width beyond the flex zone.
    pub panel_width: ParameterId,
    /// Panel thickness.
    pub thickness: ParameterId,
    /// Straight: web thickness.
    #[serde(default)]
    pub web: Option<ParameterId>,
    /// Straight: root fillet radius.
    #[serde(default)]
    pub root_radius: Option<ParameterId>,
    /// Slotted: slot pitch.
    #[serde(default)]
    pub slot_pitch: Option<ParameterId>,
    /// Slotted: slot width.
    #[serde(default)]
    pub slot_width: Option<ParameterId>,
    /// Slotted: bridge at each slot's closed end.
    #[serde(default)]
    pub bridge: Option<ParameterId>,
}

impl LivingHingeDesc {
    /// Every parameter with its role, in a fixed order.
    pub fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        let mut v = vec![
            ("length", self.length),
            ("span", self.span),
            ("panel_width", self.panel_width),
            ("thickness", self.thickness),
        ];
        for (role, id) in [
            ("web", self.web),
            ("root_radius", self.root_radius),
            ("slot_pitch", self.slot_pitch),
            ("slot_width", self.slot_width),
            ("bridge", self.bridge),
        ] {
            if let Some(id) = id {
                v.push((role, id));
            }
        }
        v
    }
}

/// Where a seam hinge sits, resolved when it was recorded (the seam
/// spans two bodies, and a feature resolves selections against its own
/// body only): a point on the seam, the seam direction, and the
/// outward and up directions (unit, orthonormal).
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct SeamDesc {
    /// Where the hinge starts on the seam.
    pub origin: [f64; 3],
    /// Along the seam.
    pub axis: [f64; 3],
    /// Away from both bodies.
    pub outward: [f64; 3],
    /// From body A into body B.
    pub up: [f64; 3],
}

/// A seam hinge's dimensions, shared by its components.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SeamHingeDesc {
    /// Along the seam.
    pub length: ParameterId,
    /// Total knuckle count (rounded at replay).
    pub knuckles: ParameterId,
    /// Pin diameter.
    pub pin_diameter: ParameterId,
    /// Knuckle wall round the bore (`None` with `knuckle_diameter`).
    pub wall: ParameterId,
    /// Knuckle outer diameter (`None`: pin + 2·clearance + 2·wall).
    #[serde(default)]
    pub knuckle_diameter: Option<ParameterId>,
    /// Radial clearance.
    pub clearance: ParameterId,
    /// Axial gap between knuckles.
    pub gap: ParameterId,
}

impl SeamHingeDesc {
    /// Every parameter with its role, in a fixed order.
    pub fn param_ids(&self) -> Vec<(&'static str, ParameterId)> {
        let mut v = vec![
            ("length", self.length),
            ("knuckles", self.knuckles),
            ("pin_diameter", self.pin_diameter),
            ("wall", self.wall),
            ("clearance", self.clearance),
            ("gap", self.gap),
        ];
        if let Some(id) = self.knuckle_diameter {
            v.push(("knuckle_diameter", id));
        }
        v
    }
}
