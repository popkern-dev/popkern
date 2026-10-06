//! Edge-hit discovery on long thread coils finds every crossing
//! (PK-679). The boolean splits a body's edges where they pierce the
//! tool's walls; a missed pierce on a face's seam edge leaves an
//! intersection record crossing the chart seam without a node, and
//! the chart walk declines. On the M2.5 x 50 rod the tuned 48-cell
//! grid found 56 of the ~110 pierces of the cylinder's seam line per
//! flank (its 3x3-station cell boxes alias on cells ~2.3 turns wide)
//! and none of the tip cone's seam pierce; the boolean's edge-hit
//! discovery now uses the dense budget.

use popkern_core::geometry::curve::CurveGeometry;
use popkern_core::geometry::intersection::surface_curve::intersect_surface_curve_dense;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;

use crate::engine::{helical_tool, plan_turns};
use crate::spec::ThreadSpec;
use crate::thread_form;

#[test]
fn rod_seam_line_pierces_every_turn_of_a_long_tool_flank() {
    let tolerance = Tolerance::default();
    let spec = match ThreadSpec::parse("M2.5") {
        Ok(s) => s,
        Err(e) => panic!("{e}"),
    };
    let profile = match thread_form::external_groove(&spec) {
        Ok(p) => p,
        Err(e) => panic!("{e}"),
    };
    let length = 49.7;
    let (z0, turns) = plan_turns(0.0, 50.0, spec.lead());
    let tool = match helical_tool(&spec, &profile, z0, turns, 0, "pierce pin", tolerance) {
        Ok(t) => t,
        Err(e) => panic!("{e}"),
    };
    // The rod cylinder's seam line: +x, at the major radius.
    let r = spec.major_diameter() / 2.0;
    let seam = CurveGeometry::Line {
        origin: Point3::new(r, 0.0, 0.0),
        direction: Vec3::new(0.0, 0.0, 1.0),
    };
    let expected = (length / spec.lead()).floor() as usize;
    let mut flanks = 0;
    for (_, face) in tool.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let result = match intersect_surface_curve_dense(
            &surface.geometry,
            (0.0, 1.0),
            (0.0, 1.0),
            &seam,
            (0.0, length),
            &tolerance,
        ) {
            Ok(r) => r,
            Err(e) => panic!("{e}"),
        };
        let n = result.intersections.len();
        // Walls the line actually crosses (the two flanks) cross it
        // once per turn; the line misses the root and crest walls.
        if n > 2 {
            flanks += 1;
            assert!(
                n >= expected && n <= expected + 1,
                "flank pierced {n} times over {length} (one per turn: {expected})"
            );
        }
    }
    assert_eq!(flanks, 2, "both flanks cross the seam line once per turn");
}
