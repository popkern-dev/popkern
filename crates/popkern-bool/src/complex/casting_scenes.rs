#![cfg(test)]
//! The six-axis arm's end castings, rebuilt for complex tests: a ring
//! (a solid disc about y) with a web through its rim at the +x extreme,
//! and a boss (a cylinder along +x) meeting both. The operands are
//! built as a script builds them: one periodic wall face per cylinder,
//! its seam on the sketch plane's u axis (the ring's on -x), the web a
//! sketched rectangle.

use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::BoolOperation;

/// One casting.
#[derive(Debug, Clone, Copy)]
pub(crate) struct Casting {
    /// Ring disc centre (x, z) and radius; the disc spans y from
    /// `ring_y` down by 10.
    pub ring: (f64, f64, f64),
    /// The ring's top face y.
    pub ring_y: f64,
    /// Web x and y extents (its z spans `web_z`).
    pub web_x: (f64, f64),
    /// Web y extent.
    pub web_y: (f64, f64),
    /// Web z extent.
    pub web_z: (f64, f64),
    /// Boss axis (y, z), radius, x extent.
    pub boss: (f64, f64, f64),
    /// Boss x extent.
    pub boss_x: (f64, f64),
    /// The boss wall's seam angle about its axis from +x of the
    /// unplaced cylinder (−π/2 puts it on −y, where the sketch plane
    /// facing +x puts it).
    pub boss_seam: f64,
}

impl Casting {
    /// The tangent boss (a): ring r 49 at x -817, web y -25..-17, boss
    /// r 33 at height `boss_y` (its top generator at boss_y + 33) from
    /// x -770 for `boss_length`.
    pub(crate) fn tangent_boss(boss_y: f64, boss_length: f64) -> Self {
        Casting {
            ring: (-817.0, 162.5, 49.0),
            ring_y: -16.0,
            web_x: (-782.95, -740.0),
            web_y: (-25.0, -17.0),
            web_z: (137.5, 187.5),
            boss: (boss_y, 162.5, 33.0),
            boss_x: (-770.0, -770.0 + boss_length),
            boss_seam: -std::f64::consts::FRAC_PI_2,
        }
    }

    /// The near-tangent cap (b): ring r 58 at x -425 (its extreme at
    /// x -367), web y -11..-3, boss r 42.5 at y -45 from `cap_x` to
    /// x -335.
    pub(crate) fn near_tangent_cap(cap_x: f64) -> Self {
        Casting {
            ring: (-425.0, 162.5, 58.0),
            ring_y: -2.0,
            web_x: (-383.95, -335.0),
            web_y: (-11.0, -3.0),
            web_z: (132.5, 192.5),
            boss: (-45.0, 162.5, 42.5),
            boss_x: (cap_x, -335.0),
            boss_seam: -std::f64::consts::FRAC_PI_2,
        }
    }

    /// The ring with its web, and the boss.
    pub(crate) fn operands(&self) -> (BRep, BRep) {
        use std::f64::consts::{FRAC_PI_2, PI};
        let tol = Tolerance::default();
        let (rx, rz, rr) = self.ring;
        let ring = cylinder(
            rr,
            10.0,
            PI,
            [
                Transform3::rotation_x(FRAC_PI_2),
                Transform3::translation(rx, self.ring_y, rz),
            ],
        );
        let web = sketched_prism(
            popkern_sketch::ProfileShape::Rect {
                width: self.web_x.1 - self.web_x.0,
                height: self.web_y.1 - self.web_y.0,
            },
            [
                0.5 * (self.web_x.0 + self.web_x.1),
                0.5 * (self.web_y.0 + self.web_y.1),
            ],
            popkern_sketch::SketchPlane {
                origin: Point3::new(0.0, 0.0, self.web_z.0),
                u_axis: Vec3::new(1.0, 0.0, 0.0),
                v_axis: Vec3::new(0.0, 1.0, 0.0),
                normal: Vec3::new(0.0, 0.0, 1.0),
            },
            self.web_z.1 - self.web_z.0,
        );
        let a = crate::boolean(BoolOperation::Union, ring, web, tol, None).expect("ring and web");
        let (by, bz, br) = self.boss;
        let boss = cylinder(
            br,
            self.boss_x.1 - self.boss_x.0,
            self.boss_seam,
            [
                Transform3::rotation_y(FRAC_PI_2),
                Transform3::translation(self.boss_x.0, by, bz),
            ],
        );
        (a, boss)
    }

    /// The union served by the complex alone (no scoping, no fallback):
    /// valid, with material in the ring, the web and the boss and none
    /// past them. Panics naming the casting otherwise.
    pub(crate) fn complex_union_serves(&self) {
        let tol = Tolerance::default();
        let (a, boss) = self.operands();
        let ix = crate::intersect::find_intersections(&a, &boss, &tol).expect("intersections");
        let (mut body, _) = crate::cellular::boolean::cellular_graph_merge(
            &a,
            &boss,
            &ix,
            BoolOperation::Union,
            &tol,
            None,
            crate::complex::serve::SelectionSource::Complex,
        )
        .unwrap_or_else(|d| panic!("the complex serves {self:?}: {}", d.into_error()));
        body.validate()
            .unwrap_or_else(|e| panic!("{self:?}: invalid: {e}"));
        // As the serve finishes it: the classifier below reads pcurves.
        let _ = body.recompute_pcurves_where_missing();
        assert!(
            crate::operation::first_missing_pcurve_face(&body).is_none(),
            "{self:?}: pcurve-incomplete"
        );
        let boxes = crate::intersect::compute_face_aabbs(&body);
        let (rx, rz, rr) = self.ring;
        let (by, bz, br) = self.boss;
        let ring_mid_y = self.ring_y - 5.0;
        let boss_mid_x = 0.5 * (self.boss_x.0 + self.boss_x.1);
        for (p, inside) in [
            ((rx, ring_mid_y, rz), true),
            (
                (
                    0.5 * (self.web_x.0 + self.web_x.1),
                    0.5 * (self.web_y.0 + self.web_y.1),
                    0.5 * (self.web_z.0 + self.web_z.1),
                ),
                true,
            ),
            ((boss_mid_x, by, bz), true),
            ((rx, ring_mid_y, rz + rr + 20.0), false),
            ((rx, self.ring_y + 6.0, rz), false),
            ((boss_mid_x, by, bz + br + 5.0), false),
            ((self.boss_x.1 + 5.0, by, bz), false),
        ] {
            let point = Point3::new(p.0, p.1, p.2);
            let got = crate::classify::point_in_solid_final(&point, &body, &boxes, &tol)
                .expect("classifies");
            let want = if inside {
                crate::PointClassification::Inside
            } else {
                crate::PointClassification::Outside
            };
            assert_eq!(got, want, "{self:?} at {p:?}");
        }
    }
}

/// `shape` at `center` sketched on `plane`, extruded along its normal
/// by `height`.
fn sketched_prism(
    shape: popkern_sketch::ProfileShape<f64>,
    center: [f64; 2],
    plane: popkern_sketch::SketchPlane,
    height: f64,
) -> BRep {
    let normal = plane.normal;
    let sheet = popkern_sketch::AnalyticProfile { shape, center }
        .to_brep(plane, Tolerance::default(), None)
        .expect("sheet");
    let face = sheet.regions[0].face;
    let mut b = sheet.brep;
    b.extrude(face, normal, height, None).expect("extrude");
    b
}

/// A solid cylinder of `radius` and `height` along z, its seam turned
/// `seam_turn` about the axis from +x, then `placed`.
fn cylinder(radius: f64, height: f64, seam_turn: f64, placed: [Transform3; 2]) -> BRep {
    let tol = Tolerance::default();
    BRep::cylinder(radius, height, tol)
        .and_then(|c| c.transform(&Transform3::rotation_z(seam_turn)))
        .and_then(|c| c.transform(&placed[0]))
        .and_then(|c| c.transform(&placed[1]))
        .expect("cylinder")
}

#[cfg(test)]
#[path = "casting_scenes_tests.rs"]
mod tests;
