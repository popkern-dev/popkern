//! The minimizer (design §8): a failing scene reduced while it still
//! fails the same oracle. The descriptor holds its placements, so the
//! reduction is the scene similarity to the identity at unit scale (the
//! contact and curve mapped back to A's frame): a failure that survives
//! it is not placement luck.

use popkern_core::geometry::{Point3, Transform3};

use super::judge::judge;
use super::scene::Scene;

/// The reduced scene and whether it still fails `oracle`.
pub fn minimize(scene: &Scene, oracle: &str) -> (Scene, bool) {
    let sim = Transform3 {
        matrix: scene.scene,
    };
    let Ok(inv) = sim.inverse() else {
        return (scene.clone(), false);
    };
    let back = |p: [f64; 3]| {
        let q = inv.apply_point(&Point3::new(p[0], p[1], p[2]));
        [q.x, q.y, q.z]
    };
    let reduced = Scene {
        scene: Transform3::identity().matrix,
        scale: 1.0,
        offset: scene.offset / scene.scale,
        contact: scene.contact.map(back),
        curve: scene.curve.iter().map(|p| back(*p)).collect(),
        ..scene.clone()
    };
    let still = judge(&reduced).failures.iter().any(|f| f.oracle == oracle);
    if still {
        (reduced, true)
    } else {
        (scene.clone(), false)
    }
}
