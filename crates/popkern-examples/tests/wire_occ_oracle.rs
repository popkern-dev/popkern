//! Wire × solid booleans against OpenCASCADE (PK-1068, plan 4.12 item
//! 1): Common(wire, solid) and Cut(wire, solid) lengths over a corpus of
//! segments, a bent polyline, a wire on a face, and circles against a
//! sphere, each in four placements (rigid motions and scales 10 and
//! 0.01). `occ_wire.py` builds the same operands from their parameters.
//! The OCC leg SKIPS when the OCP venv is absent; the popkern lengths
//! are still checked against each other (inside + outside = the wire).

use popkern_bool::{BoolOperation, boolean_with_wire};
use popkern_core::constructors::BRepBuilder;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry, CurveOnEdge};
use popkern_core::geometry::{Point3, Transform3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde_json::{Value, json};

struct Placement {
    axis: [f64; 3],
    angle: f64,
    translation: [f64; 3],
    scale: f64,
}

impl Placement {
    fn transform(&self) -> Transform3 {
        let [ax, ay, az] = self.axis;
        let [tx, ty, tz] = self.translation;
        let rotation = Transform3::rotation_around_axis(
            Point3::new(0.0, 0.0, 0.0),
            Vec3::new(ax, ay, az),
            self.angle,
        )
        .expect("rotation");
        Transform3::translation(tx, ty, tz)
            .compose(&rotation.compose(&Transform3::scale_uniform(self.scale)))
    }

    fn json(&self) -> Value {
        json!({"axis": self.axis, "angle": self.angle, "translation": self.translation, "scale": self.scale})
    }
}

const PLACEMENTS: [Placement; 4] = [
    Placement {
        axis: [0.0, 0.0, 1.0],
        angle: 0.0,
        translation: [0.0, 0.0, 0.0],
        scale: 1.0,
    },
    Placement {
        axis: [1.0, 1.0, 0.0],
        angle: 0.7,
        translation: [3.0, -2.0, 1.0],
        scale: 1.0,
    },
    Placement {
        axis: [0.0, 1.0, 1.0],
        angle: 2.1,
        translation: [-1.0, 4.0, 2.0],
        scale: 10.0,
    },
    Placement {
        axis: [1.0, 0.0, 1.0],
        angle: 1.3,
        translation: [0.02, 0.0, -0.01],
        scale: 0.01,
    },
];

fn solid(spec: &Value, tol: Tolerance) -> BRep {
    match spec["kind"].as_str() {
        Some("box") => BRep::unit_cube(tol).expect("cube"),
        Some("sphere") => BRep::sphere(spec["radius"].as_f64().expect("r"), tol).expect("sphere"),
        other => panic!("solid {other:?}"),
    }
}

fn point(v: &Value) -> Point3 {
    let f = |i: usize| v[i].as_f64().expect("number");
    Point3::new(f(0), f(1), f(2))
}

fn wire(spec: &Value, tol: Tolerance) -> BRep {
    let mut b = BRepBuilder::new(tol);
    match spec["kind"].as_str() {
        Some("polyline") => {
            let points: Vec<Point3> = spec["points"]
                .as_array()
                .expect("points")
                .iter()
                .map(point)
                .collect();
            let vs: Vec<_> = points.iter().map(|&p| b.add_vertex(p)).collect();
            for k in 0..points.len() - 1 {
                let (eh, _, _) = b.add_edge(vs[k], vs[k + 1]);
                let segment =
                    CurveOnEdge::line_segment(points[k], points[k + 1], 1e-12).expect("segment");
                b.set_edge_curve(eh, segment);
            }
        }
        Some("circle") => {
            let geometry = CurveGeometry::Circle {
                center: point(&spec["centre"]),
                normal: Vec3::new(0.0, 0.0, 1.0),
                radius: spec["radius"].as_f64().expect("r"),
            };
            let v = b.add_vertex(geometry.point_at(0.0));
            let (eh, _, _) = b.add_edge(v, v);
            b.set_edge_curve(
                eh,
                CurveOnEdge {
                    geometry,
                    t_start: 0.0,
                    t_end: std::f64::consts::TAU,
                },
            );
        }
        other => panic!("wire {other:?}"),
    }
    b.build_unchecked()
}

/// Total length of a wire body's edges, by fine chords (2e4 per edge:
/// a circle's chord sum is within 1e-8 of its length).
fn length(w: &BRep) -> f64 {
    w.edges
        .values()
        .map(|e| {
            let c = e.curve.as_ref().expect("curve");
            let n = 20_000;
            (0..n)
                .map(|k| {
                    let t = |k: usize| c.t_start + (c.t_end - c.t_start) * k as f64 / n as f64;
                    c.geometry
                        .point_at(t(k))
                        .distance_to(&c.geometry.point_at(t(k + 1)))
                })
                .sum::<f64>()
        })
        .sum()
}

fn oracle(spec: &Value) -> Option<Value> {
    let python = std::path::Path::new(&std::env::var("HOME").unwrap_or_default())
        .join("venvs/ocp/bin/python");
    let script =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tools/occ-oracle/occ_wire.py");
    if !python.exists() || !script.exists() {
        return None;
    }
    let out = std::process::Command::new(&python)
        .arg(&script)
        .arg(spec.to_string())
        .output()
        .expect("oracle");
    Some(serde_json::from_slice(&out.stdout).expect("oracle json"))
}

#[test]
fn wire_booleans_match_occ() {
    let tol = Tolerance::default();
    let cube = json!({"kind": "box", "corners": [0, 0, 0, 1, 1, 1]});
    let ball = json!({"kind": "sphere", "centre": [0, 0, 0], "radius": 1.0});
    let poly = |points: Value| json!({"kind": "polyline", "points": points});
    let corpus = [
        (
            "segment_through",
            &cube,
            poly(json!([[-1, 0.5, 0.5], [2, 0.5, 0.5]])),
        ),
        (
            "segment_ending_inside",
            &cube,
            poly(json!([[-1, 0.5, 0.5], [0.5, 0.5, 0.5]])),
        ),
        (
            "bent",
            &cube,
            poly(json!([[-1, 0.5, 0.5], [0.5, 0.5, 0.5], [0.5, 0.5, 2]])),
        ),
        ("on_a_face", &cube, poly(json!([[-1, 0.5, 1], [2, 0.5, 1]]))),
        (
            "circle_crossing_a_sphere",
            &ball,
            json!({"kind": "circle", "centre": [0.5, 0, 0.3], "radius": 0.8}),
        ),
        (
            "circle_outside_a_sphere",
            &ball,
            json!({"kind": "circle", "centre": [0, 0, 0.3], "radius": 1.5}),
        ),
        (
            "segment_through_the_poles",
            &ball,
            poly(json!([[0, 0, -2], [0, 0, 2]])),
        ),
    ];
    let mut worst = (0.0f64, String::new());
    let mut checked = 0;
    for (name, solid_spec, wire_spec) in &corpus {
        for (k, placement) in PLACEMENTS.iter().enumerate() {
            let row = format!("{name}@{k}");
            let t = placement.transform();
            let s = solid(solid_spec, tol).transform(&t).expect("place solid");
            let w = wire(wire_spec, tol).transform(&t).expect("place wire");
            let inside = boolean_with_wire(BoolOperation::Intersection, &w, &s, tol)
                .unwrap_or_else(|e| panic!("{row}: {e}"));
            let outside = boolean_with_wire(BoolOperation::Subtraction, &w, &s, tol)
                .unwrap_or_else(|e| panic!("{row}: {e}"));
            let (li, lo, lw) = (length(&inside), length(&outside), length(&w));
            let floor = 1e-9 * placement.scale;
            assert!(
                (li + lo - lw).abs() <= 1e-7 * lw.max(floor),
                "{row}: inside {li} + outside {lo} != wire {lw}"
            );
            let spec =
                json!({"solid": solid_spec, "wire": wire_spec, "placement": placement.json()});
            let Some(occ) = oracle(&spec) else {
                continue;
            };
            assert_eq!(occ["done"], true, "{row}: OCC {occ}");
            for (ours, key) in [(li, "common_length"), (lo, "cut_length")] {
                let theirs = occ[key].as_f64().expect(key);
                let g = (ours - theirs).abs() / ours.abs().max(theirs.abs()).max(floor);
                assert!(g <= 1e-6, "{row} {key}: {ours} vs OCC {theirs} ({g:.2e})");
                if g > worst.0 {
                    worst = (g, format!("{row} {key}"));
                }
            }
            checked += 1;
        }
    }
    eprintln!(
        "wire OCC oracle: {checked} rows, worst {:.3e} at {}",
        worst.0, worst.1
    );
}
