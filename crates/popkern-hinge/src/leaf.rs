//! Building one hinge component (design §5.2–§5.3): a leaf is its
//! knuckle-and-plate profile extruded along the axis, minus one comb
//! over the other leaf's knuckle zones, minus the bore (or plus the
//! integral pin), minus the screw holes; leaf B is the mirror image,
//! turned to the open angle.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fastener::{ClearanceFit, ClearanceHoleOptions, HoleDepth};
use popkern_math::PortableMath;
use popkern_sketch::{PathSegment, ProfileShape};

use crate::HingeError;
use crate::construct::{extrude, placed, plane, subtract, union, z_cylinder};
use crate::params::{HingeComponent, HingeKind, HingeLayout, HingeParams, hinge_layout};

/// Build one component of a `kind` hinge, in the hinge frame.
///
/// # Errors
///
/// The specification variants of [`HingeError`] from
/// [`hinge_layout`]; `InvalidSpec` for a pin of a print-in-place hinge
/// (it has none); construction variants naming the stage when the
/// kernel declines.
pub fn hinge_component(
    kind: HingeKind,
    component: HingeComponent,
    params: &HingeParams,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    let layout = hinge_layout(kind, params)?;
    match component {
        HingeComponent::Pin => match kind {
            HingeKind::Butt => z_cylinder(
                params.pin_diameter / 2.0,
                0.0,
                params.length,
                "pin",
                tolerance,
            ),
            HingeKind::PrintInPlace => Err(HingeError::InvalidSpec {
                field: "component",
                value: 0.0,
                requirement: "leaf_a or leaf_b: a print-in-place hinge's pin is part of leaf A",
            }),
        },
        leaf => build_leaf(kind, leaf, params, &layout, tolerance),
    }
}

/// +1 for leaf B (plate toward +X), −1 for leaf A.
fn side(component: HingeComponent) -> f64 {
    if component == HingeComponent::LeafB {
        1.0
    } else {
        -1.0
    }
}

fn build_leaf(
    kind: HingeKind,
    component: HingeComponent,
    p: &HingeParams,
    layout: &HingeLayout,
    tolerance: Tolerance,
) -> Result<BRep, HingeError> {
    let sx = side(component);
    let r = layout.knuckle_radius;
    // Tool overshoot past the body: half a knuckle radius.
    let margin = 0.5 * r;
    let stage = if component == HingeComponent::LeafA {
        "leaf A profile"
    } else {
        "leaf B profile"
    };
    let mut leaf = extrude(
        leaf_profile(sx, r, p.thickness, p.leaf_width),
        plane(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
        ),
        p.length,
        stage,
        tolerance,
    )?;
    let other = if component == HingeComponent::LeafA {
        HingeComponent::LeafB
    } else {
        HingeComponent::LeafA
    };
    let comb = extrude(
        comb_profile(sx, layout, &layout.owned(other), p.length, margin, p.gap),
        plane(
            Point3::new(0.0, r + margin, 0.0),
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ),
        2.0 * (r + margin),
        "comb notch tool",
        tolerance,
    )?;
    leaf = subtract(&leaf, &comb, "comb notch", tolerance)?;
    if kind == HingeKind::PrintInPlace && component == HingeComponent::LeafA {
        let own = layout.owned(component);
        let (first, last) = match (own.first(), own.last()) {
            (Some(f), Some(l)) => (*f, *l),
            _ => {
                return Err(crate::construct::construction(
                    "integral pin",
                    "leaf A owns no knuckle".into(),
                ));
            }
        };
        // From inside A's first knuckle to inside its last, or out to
        // the end when B owns the last knuckle — never flush with an A
        // knuckle face.
        let a_owns_last = layout.segments.len() % 2 == 1;
        let z0 = 0.5 * (first.0 + first.1);
        let z1 = if a_owns_last {
            0.5 * (last.0 + last.1)
        } else {
            p.length
        };
        let pin = z_cylinder(p.pin_diameter / 2.0, z0, z1, "integral pin tool", tolerance)?;
        leaf = union(&leaf, &pin, "integral pin", tolerance)?;
    } else {
        let bore = z_cylinder(
            layout.bore_radius,
            -margin,
            p.length + margin,
            "bore tool",
            tolerance,
        )?;
        leaf = subtract(&leaf, &bore, "bore", tolerance)?;
    }
    if let Some(holes) = &p.holes {
        let options = ClearanceHoleOptions {
            fit: ClearanceFit::Normal,
            countersink: holes.countersink.then_some(0.0),
            counterbore: None,
            entry_chamfer: false,
        };
        for hole in &layout.holes {
            leaf = popkern_fastener::holes::clearance_hole(
                &leaf,
                &holes.thread,
                Point3::new(sx * hole.offset, -r + p.thickness, hole.z),
                Vec3::new(0.0, -1.0, 0.0),
                HoleDepth::Through,
                &options,
                tolerance,
            )?;
        }
    }
    let turn = 180.0 - p.angle_deg;
    if component == HingeComponent::LeafB && turn != 0.0 {
        leaf = placed(
            &leaf,
            &Transform3::rotation_z(turn.to_radians()),
            "leaf B pose",
        )?;
    }
    Ok(leaf)
}

/// The leaf's cross-section: the mounting line from the outer edge to
/// below the axis, the knuckle arc round the axis (split under 90° per
/// piece), the front line back out. `sx` = −1 for leaf A (plate toward
/// −X), +1 for leaf B (the mirror image).
fn leaf_profile(sx: f64, r: f64, t: f64, w: f64) -> ProfileShape<f64> {
    // Authored as leaf A (plate toward −X); leaf B is its mirror.
    let pt = |x: f64, y: f64| [-sx * x, y];
    // Leaf A coordinates (plate toward −X); the arc runs from straight
    // below the axis (−90°) counter-clockwise round the far side to the
    // front line y = −R + t on the plate side.
    let y_front = -r + t;
    let x_front = -(r * r - y_front * y_front).max(0.0).sqrt();
    let start = -std::f64::consts::FRAC_PI_2;
    let mut end = y_front.patan2(x_front);
    if end < 0.0 {
        end += 2.0 * std::f64::consts::PI;
    }
    // Piece boundaries at the quarter turns: every piece stays inside
    // one quadrant, so none straddles the cylinder chart's u seam.
    let quarter = std::f64::consts::FRAC_PI_2;
    let mut cuts = vec![start];
    let mut a = start + quarter;
    while a < end - 1e-9 {
        cuts.push(a);
        a += quarter;
    }
    cuts.push(end);
    let mut segments = vec![PathSegment::Line { to: pt(0.0, -r) }];
    for (k, pair) in cuts.windows(2).enumerate() {
        let mid = 0.5 * (pair[0] + pair[1]);
        let to = if k + 2 == cuts.len() {
            pt(x_front, y_front)
        } else {
            pt(r * pair[1].pcos(), r * pair[1].psin())
        };
        segments.push(PathSegment::Arc {
            via: pt(r * mid.pcos(), r * mid.psin()),
            to,
        });
    }
    segments.push(PathSegment::Line {
        to: pt(-w, y_front),
    });
    ProfileShape::Path {
        start: pt(-w, -r),
        segments,
    }
}

/// The comb removing the other leaf's knuckle zones from this leaf, in
/// the (x, z) plane: one tooth per zone, from the notch line across the
/// axis to a spine beyond the knuckles (outside the leaf). Teeth at the
/// hinge ends run past the body so no tooth face sits on an end face.
fn comb_profile(
    sx: f64,
    layout: &HingeLayout,
    zones: &[(f64, f64)],
    length: f64,
    margin: f64,
    gap: f64,
) -> ProfileShape<f64> {
    let r = layout.knuckle_radius;
    // Leaf A coordinates: teeth from x = −notch to the spine at
    // x ∈ [R + margin, R + 2·margin]; mirrored for leaf B.
    let x_tooth = -layout.notch;
    let x_spine_in = r + margin;
    let x_spine_out = r + 2.0 * margin;
    let teeth: Vec<(f64, f64)> = zones
        .iter()
        .map(|&(z0, z1)| {
            let end = length * 1e-9;
            let lo = if z0 <= end { -margin } else { z0 - gap };
            let hi = if z1 >= length - end {
                length + margin
            } else {
                z1 + gap
            };
            (lo, hi)
        })
        .collect();
    let mut points = Vec::new();
    if let (Some(first), Some(last)) = (teeth.first(), teeth.last()) {
        points.push([sx * -x_spine_out, first.0]);
        points.push([sx * -x_spine_out, last.1]);
        for j in (0..teeth.len()).rev() {
            points.push([sx * -x_tooth, teeth[j].1]);
            points.push([sx * -x_tooth, teeth[j].0]);
            if j > 0 {
                points.push([sx * -x_spine_in, teeth[j].0]);
                points.push([sx * -x_spine_in, teeth[j - 1].1]);
            }
        }
    }
    ProfileShape::Polygon { points }
}
