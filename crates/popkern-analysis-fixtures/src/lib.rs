//! Generates the frozen STEP fixtures from the pinned popkern.
//!
//! Fixtures are regenerated only on a pin bump (`cargo run -p
//! popkern-analysis-fixtures`), and the diff is recorded in
//! `PINNED.md`. The drift test regenerates every fixture in memory and
//! requires the bytes to equal the frozen file, so a popkern
//! construction change shows up at the bump instead of hiding behind a
//! frozen file.

#![warn(missing_docs)]
#![forbid(unsafe_code)]

use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_io::{
    AssemblyChild, OpaqueBag, ProductMetadata, RepresentationContext, ShapeRep, ShapeRepKind,
    StepDocument, StepExportOptions, StepProductNode, UnitContext, export_step,
};
use popkern_part::Part;
use std::sync::Arc;

/// Why a fixture could not be built.
#[derive(Debug, thiserror::Error)]
pub enum FixtureError {
    /// A kernel constructor or transform failed.
    #[error("{fixture}: kernel: {detail}")]
    Kernel {
        /// The fixture.
        fixture: &'static str,
        /// The kernel's message.
        detail: String,
    },
    /// The part API refused the body or its STEP export.
    #[error("{fixture}: part: {detail}")]
    Part {
        /// The fixture.
        fixture: &'static str,
        /// The part API's message.
        detail: String,
    },
}

/// A fixture: its file stem and how to produce its STEP bytes.
pub struct Fixture {
    /// File stem under `fixtures/` (and the STEP product name).
    pub name: &'static str,
    /// Produces the STEP file.
    pub step: fn() -> Result<Vec<u8>, FixtureError>,
}

/// Every generated fixture, in a fixed order.
pub const FIXTURES: [Fixture; 10] = [
    Fixture {
        name: "block_two_bores",
        step: || part_step("block_two_bores", block_two_bores()?),
    },
    Fixture {
        name: "block_dimple",
        step: || part_step("block_dimple", block_dimple()?),
    },
    Fixture {
        name: "assembly_two_dimples",
        step: assembly_two_dimples,
    },
    Fixture {
        name: "two_boxes_one_body",
        step: two_boxes_one_body,
    },
    Fixture {
        name: "hinge_print_in_place_60",
        step: hinge_print_in_place_60,
    },
    Fixture {
        name: "stepped_bar_r1",
        step: || part_step("stepped_bar_r1", stepped_bar("stepped_bar_r1", 1.0)?),
    },
    Fixture {
        name: "stepped_bar_r2",
        step: || part_step("stepped_bar_r2", stepped_bar("stepped_bar_r2", 2.0)?),
    },
    Fixture {
        name: "stepped_bar_r3",
        step: || part_step("stepped_bar_r3", stepped_bar("stepped_bar_r3", 3.0)?),
    },
    Fixture {
        name: "stepped_bar_r5",
        step: || part_step("stepped_bar_r5", stepped_bar("stepped_bar_r5", 5.0)?),
    },
    Fixture {
        name: "l_bracket",
        step: || part_step("l_bracket", l_bracket()?),
    },
];

/// `HingeParams::print_in_place(60)`: two leaves, 0.3 mm clearance, as
/// two solids of one product (cross-cutting Thread 3, gate E5).
pub fn hinge_print_in_place_60() -> Result<Vec<u8>, FixtureError> {
    use popkern_hinge::{HingeComponent, HingeKind, HingeParams, hinge_component};
    const NAME: &str = "hinge_print_in_place_60";
    let params = HingeParams::print_in_place(60.0);
    let leaf = |c| {
        hinge_component(HingeKind::PrintInPlace, c, &params, Tolerance::default()).map_err(|e| {
            FixtureError::Kernel {
                fixture: NAME,
                detail: e.to_string(),
            }
        })
    };
    let leaves = vec![leaf(HingeComponent::LeafA)?, leaf(HingeComponent::LeafB)?];
    let doc = StepDocument {
        products: vec![product(NAME, Vec::new(), leaves)],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    export_step(&doc, &StepExportOptions::default()).map_err(|e| FixtureError::Part {
        fixture: NAME,
        detail: e.to_string(),
    })
}

fn part_step(name: &'static str, part: Part) -> Result<Vec<u8>, FixtureError> {
    part.step_bytes(name).map_err(|e| FixtureError::Part {
        fixture: name,
        detail: e.to_string(),
    })
}

fn product(name: &str, children: Vec<AssemblyChild>, solids: Vec<BRep>) -> StepProductNode {
    StepProductNode {
        name: name.to_string(),
        id: name.to_string(),
        children,
        shapes: if solids.is_empty() {
            Vec::new()
        } else {
            vec![Arc::new(ShapeRep {
                kind: ShapeRepKind::AdvancedBrep,
                solids,
                context: RepresentationContext::default(),
                non_brep_items: Vec::new(),
            })]
        },
        metadata: ProductMetadata::default(),
    }
}

/// `block_dimple` placed twice: at the origin, and rotated 90° about z
/// then moved to (30, 5, 0). One product definition, two instances.
pub fn assembly_two_dimples() -> Result<Vec<u8>, FixtureError> {
    const NAME: &str = "assembly_two_dimples";
    let part = block_dimple()?;
    let turned = Transform3::translation(30.0, 5.0, 0.0)
        .compose(&Transform3::rotation_z(std::f64::consts::FRAC_PI_2));
    let child = |transform, name: &str| AssemblyChild {
        product_index: 1,
        transform,
        instance_name: Some(name.to_string()),
    };
    let doc = StepDocument {
        products: vec![
            product(
                NAME,
                vec![
                    child(Transform3::identity(), "left"),
                    child(turned, "right"),
                ],
                Vec::new(),
            ),
            product("block_dimple", Vec::new(), vec![part.brep().clone()]),
        ],
        roots: vec![0],
        units: UnitContext::default(),
        opaque_bag: OpaqueBag::default(),
    };
    export_step(&doc, &StepExportOptions::default()).map_err(|e| FixtureError::Part {
        fixture: NAME,
        detail: e.to_string(),
    })
}

/// Two disjoint boxes (a 10 mm cube and a 2 mm cube at x = 15) as two
/// solids of one body.
pub fn two_boxes_one_body() -> Result<Vec<u8>, FixtureError> {
    const NAME: &str = "two_boxes_one_body";
    let a = cube_10(NAME)?;
    let b = kernel(NAME, BRep::unit_cube(Tolerance::default()))?;
    let b = kernel(NAME, b.transform(&Transform3::scale_uniform(2.0)))?;
    let b = Part::from_brep(b, Tolerance::default(), NAME).map_err(|e| FixtureError::Part {
        fixture: NAME,
        detail: e.to_string(),
    })?;
    let b = kernel(
        NAME,
        b.brep().transform(&Transform3::translation(15.0, 0.0, 0.0)),
    )?;
    let body = kernel(NAME, popkern_bool::concatenate_solids(a, b))?;
    let part =
        Part::from_brep(body, Tolerance::default(), NAME).map_err(|e| FixtureError::Part {
            fixture: NAME,
            detail: e.to_string(),
        })?;
    part_step(NAME, part)
}

fn kernel<T, E: std::fmt::Display>(
    fixture: &'static str,
    r: Result<T, E>,
) -> Result<T, FixtureError> {
    r.map_err(|e| FixtureError::Kernel {
        fixture,
        detail: e.to_string(),
    })
}

/// A 10 mm cube at the origin.
fn cube_10(fixture: &'static str) -> Result<BRep, FixtureError> {
    let cube = kernel(fixture, BRep::unit_cube(Tolerance::default()))?;
    kernel(fixture, cube.transform(&Transform3::scale_uniform(10.0)))
}

fn subtract(fixture: &'static str, a: BRep, b: BRep) -> Result<Part, FixtureError> {
    let out = kernel(
        fixture,
        popkern_bool::boolean(
            popkern_bool::BoolOperation::Subtraction,
            a,
            b,
            Tolerance::default(),
            None,
        ),
    )?;
    Part::from_brep(out, Tolerance::default(), fixture).map_err(|e| FixtureError::Part {
        fixture,
        detail: e.to_string(),
    })
}

/// A 10 mm cube with two vertical through-bores of radius 1.5 at
/// (2.5, 5) and (7.5, 5): planes and two equal cylinders, so a
/// radius-only cylinder selector matches two faces.
pub fn block_two_bores() -> Result<Part, FixtureError> {
    const NAME: &str = "block_two_bores";
    let mut body = cube_10(NAME)?;
    for x in [2.5, 7.5] {
        let drill = kernel(NAME, BRep::cylinder(1.5, 20.0, Tolerance::default()))?;
        let drill = kernel(
            NAME,
            drill.transform(&Transform3::translation(x, 5.0, -5.0)),
        )?;
        body = subtract(NAME, body, drill)?.brep().clone();
    }
    Part::from_brep(body, Tolerance::default(), NAME).map_err(|e| FixtureError::Part {
        fixture: NAME,
        detail: e.to_string(),
    })
}

/// A 10 mm cube with a spherical dimple: a sphere of radius 3 centred
/// at (5, 5, 11) removed from the top face, leaving a concave
/// spherical face.
pub fn block_dimple() -> Result<Part, FixtureError> {
    const NAME: &str = "block_dimple";
    let sphere = kernel(NAME, BRep::sphere(3.0, Tolerance::default()))?;
    let sphere = kernel(
        NAME,
        sphere.transform(&Transform3::translation(5.0, 5.0, 11.0)),
    )?;
    subtract(NAME, cube_10(NAME)?, sphere)
}

/// The polygon `outline` on the XY plane extruded by `depth` along +z.
fn prism(fixture: &'static str, outline: &[(f64, f64)], depth: f64) -> Result<Part, FixtureError> {
    use popkern_part::{Plane, Profile};
    part_op(
        fixture,
        Profile::on(Plane::XY).polygon(outline).extrude(depth),
    )
}

fn part_op(
    fixture: &'static str,
    r: Result<Part, popkern_part::PartError>,
) -> Result<Part, FixtureError> {
    r.map_err(|e| FixtureError::Part {
        fixture,
        detail: e.to_string(),
    })
}

/// A flat stepped bar (Peterson's shoulder fillet in tension, D/d = 3):
/// a 30 × 30 × 5 wide end (x ∈ [0, 30], y ∈ [−15, 15]) and a 30 × 10 × 5
/// narrow end (x ∈ [30, 60], y ∈ [−5, 5]), with both shoulder corners
/// at (30, ±5) filleted at `radius` (the shoulder is 10 wide, so any
/// radius below 10 leaves a flat shoulder face).
pub fn stepped_bar(fixture: &'static str, radius: f64) -> Result<Part, FixtureError> {
    stepped_bar_with(
        fixture,
        &SteppedBar {
            big_d: 30.0,
            step: 10.0,
            radius,
            wide_length: 30.0,
            narrow_length: 30.0,
            thickness: 5.0,
            split_end: false,
        },
    )
}

/// An L-bracket extruded 20 along z: a 40 × 5 base leg and a 5 × 30
/// upright leg, so the inner corner is the concave edge along z at
/// (5, 5) with a dihedral of π/2 between outward normals.
pub fn l_bracket() -> Result<Part, FixtureError> {
    let outline = [
        (0.0, 0.0),
        (40.0, 0.0),
        (40.0, 5.0),
        (5.0, 5.0),
        (5.0, 30.0),
        (0.0, 30.0),
    ];
    prism("l_bracket", &outline, 20.0)
}

/// Plan dimensions of a flat stepped bar (Roark Table 17.1 case 5):
/// wide width `big_d`, narrow width `big_d − 2·step`, both ends centred
/// on y = 0, the shoulder at x = `wide_length`, extruded `thickness`
/// along +z from z = 0.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SteppedBar {
    /// Wide width D.
    pub big_d: f64,
    /// Shoulder height h = (D − d)/2.
    pub step: f64,
    /// Fillet radius r.
    pub radius: f64,
    /// Length of the wide end.
    pub wide_length: f64,
    /// Length of the narrow end.
    pub narrow_length: f64,
    /// Thickness.
    pub thickness: f64,
    /// Split the narrow end face at y = 0 into two faces by extending
    /// the upper half 1 further along x, so opposite tractions on the
    /// halves apply a couple (in-plane bending). The kernel merges a
    /// collinear split, hence the step.
    pub split_end: bool,
}

/// The stepped bar `bar`, with both shoulder corners filleted by
/// popkern's fillet operator.
pub fn stepped_bar_with(fixture: &'static str, bar: &SteppedBar) -> Result<Part, FixtureError> {
    let (w, n) = (0.5 * bar.big_d, 0.5 * bar.big_d - bar.step);
    let (x0, x1) = (bar.wide_length, bar.wide_length + bar.narrow_length);
    let mut outline = vec![(0.0, -w), (x0, -w), (x0, -n), (x1, -n)];
    let x2 = if bar.split_end {
        outline.extend([(x1, 0.0), (x1 + 1.0, 0.0)]);
        x1 + 1.0
    } else {
        x1
    };
    outline.extend([(x2, n), (x0, n), (x0, w), (0.0, w)]);
    let part = prism(fixture, &outline, bar.thickness)?;
    part_op(
        fixture,
        part.fillet(bar.radius, |e| {
            e.straight()
                .direction(popkern_part::Dir::Z)
                .filter(move |i| {
                    (i.midpoint.x - x0).abs() < 1e-6 && (i.midpoint.y.abs() - n).abs() < 1e-6
                })
        }),
    )
}

/// A bar of width `big_d` (y ∈ [−D/2, D/2]) and length `length`
/// (x from 0) with one semicircular notch of radius `radius` in the
/// +y edge at mid-length (Roark Table 17.1 case 3, h = r), extruded
/// `thickness` along +z. The x = `length` end is split at y = 0 into two
/// faces (the upper half extends 1 further), so tractions on the halves
/// can place the load line or apply a couple.
pub fn u_notch_bar(
    fixture: &'static str,
    big_d: f64,
    radius: f64,
    length: f64,
    thickness: f64,
) -> Result<Part, FixtureError> {
    use popkern_part::{Plane, Profile};
    let (w, c) = (0.5 * big_d, 0.5 * length);
    // Sketched with an arc (an exact cylinder after extrusion): a
    // boolean cut leaves the cylinder's seam on the bar face, which the
    // pinned merge declines.
    let profile = Profile::on(Plane::XY)
        .path(0.0, -w)
        .line_to(length, -w)
        .line_to(length, 0.0)
        .line_to(length + 1.0, 0.0)
        .line_to(length + 1.0, w)
        .line_to(c + radius, w)
        .arc_to((c, w - radius), (c - radius, w))
        .line_to(0.0, w)
        .close();
    part_op(fixture, profile.extrude(thickness))
}

/// STEP bytes of fixture `f`.
pub fn step_bytes(f: &Fixture) -> Result<Vec<u8>, FixtureError> {
    (f.step)()
}
