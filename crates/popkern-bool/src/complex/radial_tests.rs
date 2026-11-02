//! Radial ordering on hand-built fans (plan §5, Stage 4 gate).
//!
//! Each fixture is one complex edge with a fan of faces around it. The
//! faces' surfaces sit on the faces of a placeholder operand (radial
//! ordering reads only a face's surface), and each fin's sense is
//! chosen so its face leaves the edge toward a given side.

use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::nurbs_surface::NurbsSurface;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use super::super::ids::{EdgeId, FinId};
use super::super::model::{CellComplex, FaceOwner};
use super::super::radial_ladder::TieLadder;
use super::super::radial_stations::same_cycle;
use super::{FaceGeometry, RadialDecline, order_radially_with, station};

const TIE: f64 = 1e-4;
const LINEAR: f64 = 1e-6;

/// One face of a fan: its surface and the side it leaves the edge
/// toward at the edge's midpoint.
struct Sheet {
    surface: SurfaceGeometry,
    side: Vec3,
}

fn plane_z0() -> SurfaceGeometry {
    SurfaceGeometry::Plane {
        origin: Point3::new(0.0, 0.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
    }
}

/// A cylinder along x resting on z = 0 along the x axis.
fn cylinder_on_x_axis(radius: f64) -> SurfaceGeometry {
    SurfaceGeometry::Cylinder {
        origin: Point3::new(0.0, 0.0, radius),
        axis: Vec3::new(1.0, 0.0, 0.0),
        radius,
    }
}

/// The fan's operand and complex: one edge `a → b` (closed when they
/// coincide, along `curve`), one face per sheet.
fn fan(
    curve: Option<CurveOnEdge>,
    a: Point3,
    b: Point3,
    sheets: &[Sheet],
) -> (BRep, CellComplex, Vec<FinId>) {
    let mut brep = BRep::unit_cube(Tolerance::default()).expect("cube");
    let handles: Vec<_> = brep.faces.keys().collect();
    assert!(sheets.len() <= handles.len());
    let mut c = CellComplex::new(1);
    let va = c.add_vertex(a, LINEAR, None);
    let vb = if a.distance_to(&b) > 0.0 {
        c.add_vertex(b, LINEAR, None)
    } else {
        va
    };
    let edge = c.add_edge(va, vb, curve, LINEAR, None);
    let mut fins = Vec::new();
    for (sheet, &fh) in sheets.iter().zip(&handles) {
        let surf = SurfaceOnFace {
            geometry: sheet.surface.clone(),
            orientation: true,
        };
        brep.faces[fh].surface = Some(surf);
        let face = c.add_face(
            vec![FaceOwner {
                operand: 0,
                face: fh,
                same_sense: true,
            }],
            true,
            Vec::new(),
        );
        let (p, t) = station(&c, edge, 0.5).expect("station");
        let surf = brep.faces[fh].surface.as_ref().expect("surface");
        let (u, v) = surf.geometry.closest_point(&p).expect("on the surface");
        let n = surf.folded_normal(u, v).expect("normal");
        let forward = n.cross(&t).dot(&sheet.side) > 0.0;
        let lp = c
            .add_loop(face, &[(edge, forward, None)], true)
            .expect("loop");
        fins.push(c.loops[lp.index()].first);
    }
    (brep, c, fins)
}

fn order(brep: &BRep, c: &mut CellComplex, ladder: TieLadder) -> Result<Vec<FinId>, RadialDecline> {
    let geo = FaceGeometry {
        operands: [brep, brep],
        graph: None,
    };
    order_radially_with(c, &geo, TIE, LINEAR, ladder)?;
    Ok(c.edges[0].radial.clone())
}

fn assert_cycle(got: &[FinId], want: &[FinId]) {
    let tag = |o: &[FinId]| o.iter().map(|&f| (f, false)).collect::<Vec<_>>();
    assert!(
        same_cycle(&tag(want), &tag(got)),
        "got {got:?}, want {want:?}"
    );
}

/// A plane and a cylinder resting on it along the x axis, both leaving
/// the line toward ±y: [plane +y, cylinder +y, cylinder −y, plane −y].
fn line_contact() -> (BRep, CellComplex, Vec<FinId>) {
    let (y, ny) = (Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
    fan(
        None,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        &[
            Sheet {
                surface: plane_z0(),
                side: y,
            },
            Sheet {
                surface: cylinder_on_x_axis(0.5),
                side: y,
            },
            Sheet {
                surface: cylinder_on_x_axis(0.5),
                side: ny,
            },
            Sheet {
                surface: plane_z0(),
                side: ny,
            },
        ],
    )
}

#[test]
fn a_line_contact_orders_by_curvature() {
    // About +x, counter-clockwise runs from +y up toward +z, where the
    // cylinder bends.
    let (brep, mut c, f) = line_contact();
    let got = order(&brep, &mut c, TieLadder::default()).expect("ordered");
    assert_cycle(&got, &[f[0], f[1], f[2], f[3]]);
}

#[test]
fn each_rung_alone_orders_the_line_contact_the_same() {
    for ladder in [
        TieLadder {
            curvature: true,
            offset_probe: false,
        },
        TieLadder {
            curvature: false,
            offset_probe: true,
        },
    ] {
        let (brep, mut c, f) = line_contact();
        let got = order(&brep, &mut c, ladder).expect("ordered");
        assert_cycle(&got, &[f[0], f[1], f[2], f[3]]);
    }
}

/// Negative control: with no rung the line contact is a tie, and the
/// decline names the tied fins, the station and what was measured.
#[test]
fn without_the_ladder_a_line_contact_declines_ambiguous() {
    let (brep, mut c, f) = line_contact();
    let off = TieLadder {
        curvature: false,
        offset_probe: false,
    };
    match order(&brep, &mut c, off) {
        Err(RadialDecline::Ambiguous {
            edge,
            fins,
            station,
            evidence,
        }) => {
            assert_eq!(edge, EdgeId::from_index(0));
            assert_eq!(fins.len(), 2);
            assert!(fins.iter().all(|x| f.contains(x)));
            assert_eq!(station, 0.5);
            assert!(
                evidence
                    .iter()
                    .all(|e| e.curvature.is_none() && e.probe.is_none())
            );
        }
        other => panic!("expected Ambiguous, got {other:?}"),
    }
}

#[test]
fn internally_tangent_cylinders_order_by_curvature() {
    // Both bend toward +z; the smaller (curvature 2) lies above the
    // larger (curvature 1) next to the line.
    let (y, ny) = (Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
    let (brep, mut c, f) = fan(
        None,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        &[
            Sheet {
                surface: cylinder_on_x_axis(0.5),
                side: y,
            },
            Sheet {
                surface: cylinder_on_x_axis(1.0),
                side: y,
            },
            Sheet {
                surface: cylinder_on_x_axis(0.5),
                side: ny,
            },
            Sheet {
                surface: cylinder_on_x_axis(1.0),
                side: ny,
            },
        ],
    );
    let got = order(&brep, &mut c, TieLadder::default()).expect("ordered");
    assert_cycle(&got, &[f[1], f[0], f[2], f[3]]);
}

#[test]
fn a_torus_lying_on_a_plane_orders_around_its_contact_circle() {
    let (major, minor) = (2.0, 0.5);
    let z = Vec3::new(0.0, 0.0, 1.0);
    let circle = CurveOnEdge {
        geometry: CurveGeometry::Circle {
            center: Point3::new(0.0, 0.0, 0.0),
            normal: z,
            radius: major,
        },
        t_start: 0.0,
        t_end: std::f64::consts::TAU,
    };
    let torus = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, minor),
        axis: z,
        major_radius: major,
        minor_radius: minor,
    };
    // The midpoint station's outward radial direction.
    let start = circle.geometry.point_at(0.0);
    let probe = fan(Some(circle.clone()), start, start, &[]);
    let (p, t) = station(&probe.1, EdgeId::from_index(0), 0.5).expect("station");
    let out = Vec3::new(p.x, p.y, 0.0).normalized().expect("off axis");
    let inward = out * -1.0;
    let (brep, mut c, f) = fan(
        Some(circle),
        start,
        start,
        &[
            Sheet {
                surface: plane_z0(),
                side: out,
            },
            Sheet {
                surface: torus.clone(),
                side: out,
            },
            Sheet {
                surface: torus,
                side: inward,
            },
            Sheet {
                surface: plane_z0(),
                side: inward,
            },
        ],
    );
    let got = order(&brep, &mut c, TieLadder::default()).expect("ordered");
    // Counter-clockwise from `out` heads up (the torus side) or down.
    if t.cross(&out).dot(&z) > 0.0 {
        assert_cycle(&got, &[f[0], f[1], f[2], f[3]]);
    } else {
        assert_cycle(&got, &[f[0], f[3], f[2], f[1]]);
    }
}

/// Negative control: two coincident planes tie on every rung.
#[test]
fn coincident_faces_decline_with_their_measured_keys() {
    let (y, ny) = (Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
    let (brep, mut c, _) = fan(
        None,
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.0, 0.0, 0.0),
        &[
            Sheet {
                surface: plane_z0(),
                side: y,
            },
            Sheet {
                surface: plane_z0(),
                side: y,
            },
            Sheet {
                surface: plane_z0(),
                side: ny,
            },
        ],
    );
    match order(&brep, &mut c, TieLadder::default()) {
        Err(RadialDecline::Ambiguous { evidence, .. }) => {
            assert_eq!(evidence.len(), 2);
            for e in evidence {
                assert!(e.curvature.is_some_and(|k| k.abs() < 1e-12), "{e:?}");
                assert!(e.probe.is_some_and(|d| d.abs() < 1e-9), "{e:?}");
            }
        }
        other => panic!("expected Ambiguous, got {other:?}"),
    }
}

/// Negative control: a face crossing another inside the edge. The
/// saddle z = xy holds the x axis and leaves it toward (0, 1, x):
/// below the plane z = 0 for x < 0, above it for x > 0. On x ∈ [−1, 2]
/// the stations sit at x = 0.5 and −0.25, either side of the crossing.
#[test]
fn faces_crossing_inside_the_edge_are_inconsistent() {
    let saddle = SurfaceGeometry::Nurbs(NurbsSurface {
        degree_u: 1,
        degree_v: 1,
        knots_u: vec![0.0, 0.0, 1.0, 1.0],
        knots_v: vec![0.0, 0.0, 1.0, 1.0],
        control_points: vec![
            Point3::new(-1.0, -1.0, 1.0),
            Point3::new(-1.0, 1.0, -1.0),
            Point3::new(2.0, -1.0, -2.0),
            Point3::new(2.0, 1.0, 2.0),
        ],
        num_u: 2,
        num_v: 2,
        weights: vec![1.0; 4],
    });
    let (y, ny) = (Vec3::new(0.0, 1.0, 0.0), Vec3::new(0.0, -1.0, 0.0));
    let (brep, mut c, f) = fan(
        None,
        Point3::new(-1.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        &[
            Sheet {
                surface: plane_z0(),
                side: y,
            },
            Sheet {
                surface: saddle,
                side: y,
            },
            Sheet {
                surface: plane_z0(),
                side: ny,
            },
        ],
    );
    match order(&brep, &mut c, TieLadder::default()) {
        Err(RadialDecline::Inconsistent { edge, stations }) => {
            assert_eq!(edge, EdgeId::from_index(0));
            assert_eq!(stations.len(), 2);
            assert_eq!(stations[0].0, 0.5);
            assert_eq!(stations[1].0, 0.25);
            assert!(stations.iter().all(|(_, o)| o.len() == f.len()));
        }
        other => panic!("expected Inconsistent, got {other:?}"),
    }
}
