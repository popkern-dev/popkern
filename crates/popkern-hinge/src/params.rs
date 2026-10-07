//! Hinge parameters, their checks, and the layout every component is
//! built from (the numbers a GUI can preview without building).

use popkern_fastener::{ClearanceFit, ThreadSpec};

use crate::HingeError;

/// Which hinge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HingeKind {
    /// Two leaves and a separate pin (a piano hinge is a butt hinge
    /// with many knuckles).
    Butt,
    /// Two leaves printed together: the pin is part of leaf A's
    /// knuckles and runs through leaf B's bored knuckles.
    PrintInPlace,
}

impl HingeKind {
    /// The skin spelling (`butt`, `print_in_place`).
    pub fn spelling(self) -> &'static str {
        match self {
            HingeKind::Butt => "butt",
            HingeKind::PrintInPlace => "print_in_place",
        }
    }

    /// The kind a skin spelling names.
    pub fn parse(text: &str) -> Option<HingeKind> {
        match text {
            "butt" => Some(HingeKind::Butt),
            "print_in_place" => Some(HingeKind::PrintInPlace),
            _ => None,
        }
    }
}

/// One body of a hinge.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HingeComponent {
    /// The leaf toward −X (owns the first knuckle; carries the pin in a
    /// print-in-place hinge).
    LeafA,
    /// The leaf toward +X (turned by the open angle).
    LeafB,
    /// The separate pin of a butt hinge.
    Pin,
}

impl HingeComponent {
    /// The skin spelling (`leaf_a`, `leaf_b`, `pin`).
    pub fn spelling(self) -> &'static str {
        match self {
            HingeComponent::LeafA => "leaf_a",
            HingeComponent::LeafB => "leaf_b",
            HingeComponent::Pin => "pin",
        }
    }

    /// The component a skin spelling names.
    pub fn parse(text: &str) -> Option<HingeComponent> {
        match text {
            "leaf_a" => Some(HingeComponent::LeafA),
            "leaf_b" => Some(HingeComponent::LeafB),
            "pin" => Some(HingeComponent::Pin),
            _ => None,
        }
    }
}

/// Screw holes through each leaf.
#[derive(Debug, Clone)]
pub struct HingeHoles {
    /// Holes per leaf, spaced evenly along it.
    pub per_leaf: u32,
    /// The screw the holes clear (ISO 273 normal fit).
    pub thread: ThreadSpec,
    /// Sink a 90° countersink for a flush countersunk head.
    pub countersink: bool,
}

/// A hinge's dimensions, in model units.
#[derive(Debug, Clone)]
pub struct HingeParams {
    /// Along the pin axis.
    pub length: f64,
    /// Each leaf, from the pin axis to its outer edge.
    pub leaf_width: f64,
    /// Leaf plate thickness.
    pub thickness: f64,
    /// Total knuckle count (leaf A owns the even ones, from z = 0).
    pub knuckles: u32,
    /// Pin diameter.
    pub pin_diameter: f64,
    /// Knuckle outer diameter (`None`: pin + 2·clearance + 2·thickness,
    /// a leaf rolled round the pin).
    pub knuckle_diameter: Option<f64>,
    /// Radial clearance: the bore over the pin, and each leaf's notch
    /// over the other leaf's knuckles.
    pub clearance: f64,
    /// Axial gap between adjacent knuckles.
    pub gap: f64,
    /// Open angle in degrees: 180 flat, 0 folded shut.
    pub angle_deg: f64,
    /// Screw holes, if any.
    pub holes: Option<HingeHoles>,
}

impl HingeParams {
    /// A butt hinge of `length` with proportional defaults (any unit):
    /// leaf width L/4, thickness L/30, pin L/15, five knuckles, an exact
    /// pin fit, gaps of L/400, flat.
    pub fn new(length: f64) -> HingeParams {
        HingeParams {
            length,
            leaf_width: length / 4.0,
            thickness: length / 30.0,
            knuckles: 5,
            pin_diameter: length / 15.0,
            knuckle_diameter: None,
            clearance: 0.0,
            gap: length / 400.0,
            angle_deg: 180.0,
            holes: None,
        }
    }

    /// A print-in-place hinge of `length`: as [`HingeParams::new`] with
    /// printed clearances and gaps of L/200.
    pub fn print_in_place(length: f64) -> HingeParams {
        HingeParams {
            clearance: length / 200.0,
            gap: length / 200.0,
            ..HingeParams::new(length)
        }
    }
}

/// Where one screw hole sits on leaf A (leaf B's mirror at `x → −x`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HolePosition {
    /// Distance of the hole axis from the pin axis.
    pub offset: f64,
    /// Position along the pin axis.
    pub z: f64,
}

/// The numbers every component is built from.
#[derive(Debug, Clone, PartialEq)]
pub struct HingeLayout {
    /// Knuckle radius R.
    pub knuckle_radius: f64,
    /// Bore radius (pin radius + clearance).
    pub bore_radius: f64,
    /// Where a leaf stops inside the other leaf's knuckle zones
    /// (distance from the axis).
    pub notch: f64,
    /// Knuckle z ranges, in order; leaf A owns the even ones.
    pub segments: Vec<(f64, f64)>,
    /// Screw hole positions on each leaf.
    pub holes: Vec<HolePosition>,
}

impl HingeLayout {
    /// The segments `component` owns (leaf A: even, leaf B: odd).
    pub fn owned(&self, component: HingeComponent) -> Vec<(f64, f64)> {
        let parity = usize::from(component == HingeComponent::LeafB);
        self.segments
            .iter()
            .enumerate()
            .filter(|(i, _)| i % 2 == parity)
            .map(|(_, s)| *s)
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

fn non_negative(field: &'static str, value: f64) -> Result<(), HingeError> {
    if value.is_finite() && value >= 0.0 {
        Ok(())
    } else {
        Err(HingeError::InvalidSpec {
            field,
            value,
            requirement: "finite and ≥ 0",
        })
    }
}

/// The most knuckles a hinge may have (a piano hinge stays well under).
pub const MAX_KNUCKLES: u32 = 400;

/// Check `params` for `kind` and lay the hinge out.
///
/// # Errors
///
/// [`HingeError`] specification variants naming the field and limit.
pub fn hinge_layout(kind: HingeKind, params: &HingeParams) -> Result<HingeLayout, HingeError> {
    let p = params;
    positive("length", p.length)?;
    positive("leaf_width", p.leaf_width)?;
    positive("thickness", p.thickness)?;
    positive("pin_diameter", p.pin_diameter)?;
    non_negative("clearance", p.clearance)?;
    non_negative("gap", p.gap)?;
    if !(p.angle_deg.is_finite() && (0.0..360.0).contains(&p.angle_deg)) {
        return Err(HingeError::InvalidSpec {
            field: "angle",
            value: p.angle_deg,
            requirement: "in [0, 360) degrees",
        });
    }
    let min_knuckles = match kind {
        HingeKind::Butt => 2,
        // Leaf B must have a knuckle between two of leaf A's to be
        // captured on the integral pin.
        HingeKind::PrintInPlace => 3,
    };
    if p.knuckles < min_knuckles || p.knuckles > MAX_KNUCKLES {
        return Err(HingeError::InvalidSpec {
            field: "knuckles",
            value: f64::from(p.knuckles),
            requirement: if kind == HingeKind::Butt {
                "between 2 and 400"
            } else {
                "between 3 and 400 (leaf B must sit between two of leaf A's knuckles)"
            },
        });
    }
    if kind == HingeKind::PrintInPlace {
        if p.clearance <= 0.0 {
            return Err(HingeError::PrintClearance {
                field: "clearance",
                value: p.clearance,
            });
        }
        if p.gap <= 0.0 {
            return Err(HingeError::PrintClearance {
                field: "gap",
                value: p.gap,
            });
        }
    }
    let n = f64::from(p.knuckles);
    let segment = (p.length - (n - 1.0) * p.gap) / n;
    if segment <= p.gap || segment <= 0.0 {
        return Err(HingeError::KnucklesDoNotFit {
            knuckles: p.knuckles,
            gap: p.gap,
            length: p.length,
            segment,
        });
    }
    let bore_radius = p.pin_diameter / 2.0 + p.clearance;
    let knuckle_radius = match p.knuckle_diameter {
        Some(d) => {
            positive("knuckle_diameter", d)?;
            d / 2.0
        }
        None => bore_radius + p.thickness,
    };
    if knuckle_radius <= bore_radius {
        return Err(HingeError::NoKnuckleWall {
            knuckle_radius,
            bore_radius,
        });
    }
    if p.thickness >= 2.0 * knuckle_radius {
        return Err(HingeError::InvalidSpec {
            field: "thickness",
            value: p.thickness,
            requirement: "less than the knuckle diameter",
        });
    }
    // An exact-fit hinge still keeps the notch off the knuckle cylinder
    // (a notch face tangent to it is a graze contact).
    let notch = knuckle_radius + p.clearance.max(p.length * 1e-4);
    if p.leaf_width <= notch {
        return Err(HingeError::PlateTooNarrow {
            leaf_width: p.leaf_width,
            needed: notch,
        });
    }
    let segments = (0..p.knuckles)
        .map(|i| {
            let z0 = f64::from(i) * (segment + p.gap);
            (z0, z0 + segment)
        })
        .collect();
    let holes = match &p.holes {
        Some(h) => hole_positions(h, p, notch)?,
        None => Vec::new(),
    };
    Ok(HingeLayout {
        knuckle_radius,
        bore_radius,
        notch,
        segments,
        holes,
    })
}

/// Holes evenly along the leaf, centred on the plate beyond the notch.
fn hole_positions(
    holes: &HingeHoles,
    p: &HingeParams,
    notch: f64,
) -> Result<Vec<HolePosition>, HingeError> {
    if holes.per_leaf == 0 {
        return Ok(Vec::new());
    }
    holes.thread.validate().map_err(HingeError::from)?;
    let bore = holes.thread.iso_clearance_diameter(ClearanceFit::Normal);
    let head = if holes.countersink {
        holes.thread.iso_countersunk_head_diameter()
    } else {
        bore
    };
    // Room around the head: a tenth of it, so the hole never grazes an
    // edge of the plate.
    let needed = head * 1.1;
    let plate = p.leaf_width - notch;
    let spacing = p.length / f64::from(holes.per_leaf);
    let fail = |dimension, needed, available| HingeError::HoleDoesNotFit {
        hole: 1,
        dimension,
        needed,
        available,
    };
    if needed >= plate {
        return Err(fail("plate width", needed, plate));
    }
    if needed >= spacing {
        return Err(fail("spacing along the leaf", needed, spacing));
    }
    if holes.countersink {
        let depth = (head - bore) / 2.0;
        if depth >= p.thickness * 0.9 {
            return Err(fail("thickness for the countersink", depth, p.thickness));
        }
    }
    let offset = notch + plate / 2.0;
    Ok((0..holes.per_leaf)
        .map(|k| HolePosition {
            offset,
            z: spacing * (f64::from(k) + 0.5),
        })
        .collect())
}
