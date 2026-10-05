//! The built thread tool walls are the helicoids they model (PK-675).
//! An iso-profile line of a swept tool wall carries one profile point
//! along the helix, so unscrewing it (rotating back by its angle and
//! dropping the matching share of the lead) must leave one fixed
//! meridian point. Its axial spread is the wall's deviation from the
//! ideal helicoid: path, frame and skinning error together.
//!
//! Measured (peak to peak, M3/M6/M10, both hands, both tools): axial
//! 2.0e-6 to 6.2e-6; radial 4.1e-5 to 1.4e-4, the skinning floor of 24
//! frames per turn (it scales with the radius). Before PK-675 the
//! quadratic-arc helix and a two-per-span arc-length table left
//! 5.7e-3 to 1.8e-2 axially and 1.3e-3 to 4.7e-3 radially.

use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;

use crate::engine::{helical_tool, plan_turns};
use crate::spec::{Handedness, ThreadSpec};
use crate::thread_form::{self, ToolProfile};

/// The largest axial and radial spread of unscrewed iso-profile lines
/// over the tool's swept walls.
fn helicoid_spread(spec: &ThreadSpec, profile: &ToolProfile) -> (f64, f64) {
    let (z0, turns) = plan_turns(0.0, 4.0 * spec.lead(), spec.lead());
    let tool = match helical_tool(
        spec,
        profile,
        z0,
        turns,
        0,
        "fidelity",
        Tolerance::default(),
    ) {
        Ok(t) => t,
        Err(e) => panic!("tool build failed: {e}"),
    };
    let winding = if spec.handedness() == Handedness::Left {
        -1.0
    } else {
        1.0
    };
    let (mut axial, mut radial) = (0.0f64, 0.0f64);
    let mut walls = 0;
    for (_, face) in tool.faces.iter() {
        let Some(surface) = face.surface.as_ref() else {
            continue;
        };
        let SurfaceGeometry::Nurbs(_) = &surface.geometry else {
            continue;
        };
        walls += 1;
        let (ur, vr) = surface.geometry.parameter_range();
        let nv = 2000;
        for i in 0..=6 {
            let u = ur.0 + (ur.1 - ur.0) * i as f64 / 6.0;
            let (mut z_lo, mut z_hi, mut r_lo, mut r_hi) = (f64::MAX, f64::MIN, f64::MAX, f64::MIN);
            let mut unwrapped = 0.0;
            let mut previous = None;
            for j in 0..=nv {
                let v = vr.0 + (vr.1 - vr.0) * j as f64 / nv as f64;
                let p = surface.geometry.point_at(u, v);
                let angle = p.y.atan2(p.x);
                if let Some(prev) = previous {
                    let mut d: f64 = angle - prev;
                    d -= (d / std::f64::consts::TAU).round() * std::f64::consts::TAU;
                    unwrapped += d;
                }
                previous = Some(angle);
                let z = p.z - winding * spec.lead() * unwrapped / std::f64::consts::TAU;
                let r = p.x.hypot(p.y);
                z_lo = z_lo.min(z);
                z_hi = z_hi.max(z);
                r_lo = r_lo.min(r);
                r_hi = r_hi.max(r);
            }
            axial = axial.max(z_hi - z_lo);
            radial = radial.max(r_hi - r_lo);
        }
    }
    assert!(walls >= 4, "expected swept walls, found {walls}");
    (axial, radial)
}

#[test]
fn thread_tool_walls_are_helicoids() {
    for designation in ["M3", "M6", "M10"] {
        let base = match ThreadSpec::parse(designation) {
            Ok(s) => s,
            Err(e) => panic!("{designation}: {e}"),
        };
        for spec in [base.clone(), base.left_handed()] {
            let profiles = [
                ("external", thread_form::external_groove(&spec)),
                ("internal", thread_form::internal_ridge(&spec, 0.0)),
            ];
            for (kind, profile) in profiles {
                let profile = match profile {
                    Ok(p) => p,
                    Err(e) => panic!("{designation} {kind}: {e}"),
                };
                let (axial, radial) = helicoid_spread(&spec, &profile);
                assert!(
                    axial < 1e-5,
                    "{designation} {kind}: wall leaves the helicoid by {axial:.3e} axially"
                );
                assert!(
                    radial < 2e-4,
                    "{designation} {kind}: wall leaves the helicoid by {radial:.3e} radially"
                );
            }
        }
    }
}
