//! `corner_cap: "rolling_ball"` (PK-861): where three rounds of one
//! radius meet on planar supports, the corner is the exact rolling-ball
//! sphere, G1 to every round.
//!
//! The default corner is the setback G2 cap (KP20 patches). The opt-in
//! keeps the setback topology, which at such a corner already ends each
//! round on the ball's cross-section through its centre, and replaces
//! the cap by one sphere face. Each body is checked against Steiner's
//! formula for a convex polyhedron rounded by `r`: the inner polyhedron
//! (faces moved in by `r`) swept by the ball,
//! `V' + S' r + sum(L (pi - a) / 2) r^2 + 4/3 pi r^3`.

use popkern_core::geometry::{SurfaceGeometry, Vec3};
use popkern_core::topology::BRep;
use popkern_rhai::ScriptHost;

const R: f64 = 4.2;
const RADIUS: f64 = 0.7;

fn rounded(normals: &str, corner_cap: &str) -> BRep {
    let source = format!(
        r#"fn dot(a, b) {{ a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }}
fn scale(a, s) {{ [a[0] * s, a[1] * s, a[2] * s] }}
fn norm(v) {{ scale(v, 1.0 / sqrt(dot(v, v))) }}
let m = model();
let big = 12.0;
let die = m.profile(xy_at(-big)).rect(2.0 * big, 2.0 * big).sketch("block").extrude(blind(2.0 * big));
let i = 0;
for raw in {normals} {{
    let n = norm(raw);
    let deep = 3.0 * big;
    let slab = m.tool(|t| t.profile(plane_at(scale(n, -deep), n)).rect(6.0 * big, 6.0 * big)
        .sketch("slab_" + i).extrude(blind(deep + {R})));
    die = die.intersect(slab, #{{ label: "cut_" + i }});
    i += 1;
}}
die = die.fillet(edges(), {RADIUS}, #{{ label: "rounds", corner_cap: "{corner_cap}" }});
m
"#
    );
    let mut run = ScriptHost::new()
        .run(&source, "rounded.rhai")
        .unwrap_or_else(|e| panic!("{corner_cap}: {}", e.render()));
    run.model.take().expect("model").to_brep().expect("brep")
}

fn spheres(brep: &BRep) -> usize {
    brep.faces
        .values()
        .filter(|f| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Sphere { .. })
            )
        })
        .count()
}

fn check(brep: &BRep, corners: usize, steiner: f64, what: &str) {
    assert_eq!(spheres(brep), corners, "{what}: one ball per corner");
    let issues = brep.validate_all();
    assert!(issues.is_empty(), "{what}: {issues:?}");
    let volume = popkern_measure::body_volume(brep, &brep.tolerance).expect("volume");
    assert!(
        ((volume - steiner) / steiner).abs() < 1e-5,
        "{what}: volume {volume} vs Steiner {steiner}"
    );
}

#[test]
fn a_rounded_tetrahedron_takes_rolling_ball_corners() {
    let rho = R - RADIUS;
    let a = rho * 2.0 * 6f64.sqrt();
    let dihedral = (1.0f64 / 3.0).acos();
    let steiner = a.powi(3) / (6.0 * 2f64.sqrt())
        + 3f64.sqrt() * a * a * RADIUS
        + 6.0 * a * (std::f64::consts::PI - dihedral) / 2.0 * RADIUS * RADIUS
        + 4.0 / 3.0 * std::f64::consts::PI * RADIUS.powi(3);
    let brep = rounded(
        "[[-1.0, -1.0, -1.0], [-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]]",
        "rolling_ball",
    );
    check(&brep, 4, steiner, "tetrahedron");
}

#[test]
fn a_rounded_dodecahedron_takes_rolling_ball_corners() {
    let phi = (1.0 + 5f64.sqrt()) / 2.0;
    let normals = format!(
        "[[0.0, 1.0, {phi}], [0.0, -1.0, {phi}], [0.0, 1.0, -{phi}], [0.0, -1.0, -{phi}], \
         [1.0, {phi}, 0.0], [-1.0, {phi}, 0.0], [1.0, -{phi}, 0.0], [-1.0, -{phi}, 0.0], \
         [{phi}, 0.0, 1.0], [-{phi}, 0.0, 1.0], [{phi}, 0.0, -1.0], [-{phi}, 0.0, -1.0]]"
    );
    let s5 = 5f64.sqrt();
    let a = (R - RADIUS) / ((250.0 + 110.0 * s5).sqrt() / 20.0);
    let dihedral = (-1.0 / s5).acos();
    let steiner = (15.0 + 7.0 * s5) / 4.0 * a.powi(3)
        + 3.0 * (25.0 + 10.0 * s5).sqrt() * a * a * RADIUS
        + 30.0 * a * (std::f64::consts::PI - dihedral) / 2.0 * RADIUS * RADIUS
        + 4.0 / 3.0 * std::f64::consts::PI * RADIUS.powi(3);
    // PK-864: the carve takes the ball directly (4.7 s here; the KP20
    // cap it once fitted and threw away took the corner to 15-25 s).
    let watch = std::time::Instant::now();
    let brep = rounded(&normals, "rolling_ball");
    let took = watch.elapsed().as_secs_f64();
    assert!(
        took < 10.0,
        "dodecahedron rolling-ball fillet took {took:.1} s"
    );
    check(&brep, 20, steiner, "dodecahedron");
    volume_does_not_depend_on_the_ball_chart(&brep, steiner);
}

/// Cube corners are already the exact ball on every path: the opt-in
/// leaves the served body as the default serves it.
#[test]
fn a_rounded_cube_is_unchanged_by_the_opt_in() {
    let normals = "[[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, -1.0, 0.0], \
                   [0.0, 0.0, 1.0], [0.0, 0.0, -1.0]]";
    let (auto, ball) = (rounded(normals, "auto"), rounded(normals, "rolling_ball"));
    let a = 2.0 * (R - RADIUS);
    let steiner = a.powi(3)
        + 6.0 * a * a * RADIUS
        + 12.0 * a * std::f64::consts::FRAC_PI_4 * RADIUS * RADIUS
        + 4.0 / 3.0 * std::f64::consts::PI * RADIUS.powi(3);
    check(&ball, 8, steiner, "cube");
    let volume = |b: &BRep| popkern_measure::body_volume(b, &b.tolerance).expect("volume");
    assert_eq!(
        volume(&auto).to_bits(),
        volume(&ball).to_bits(),
        "cube: identical serve"
    );
}

/// A cap's chart pole from its direction.
type ChartAxis = Box<dyn Fn(Vec3) -> Vec3>;

/// Every ball cap re-charted with its pole along `axis(cap direction)`
/// and fresh pcurves (the geometry is unchanged: only the chart moves).
fn recharted(brep: &BRep, axis: impl Fn(Vec3) -> Vec3) -> BRep {
    let mut out = brep.clone();
    let caps: Vec<_> = out
        .faces
        .iter()
        .filter(|(_, f)| {
            matches!(
                f.surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Sphere { .. })
            )
        })
        .map(|(h, _)| h)
        .collect();
    for &face in &caps {
        let Some(SurfaceGeometry::Sphere { center, .. }) =
            out.faces[face].surface.as_ref().map(|s| s.geometry.clone())
        else {
            continue;
        };
        let start = out.loops[out.faces[face].outer_loop].half_edge;
        let mut toward = Vec3::new(0.0, 0.0, 0.0);
        let mut he = start;
        loop {
            let h = &out.half_edges[he];
            toward = toward + (out.vertices[h.origin].position - center);
            he = h.next;
            if he == start {
                break;
            }
        }
        let pole = axis(toward.normalized().expect("cap direction"));
        if let Some(SurfaceGeometry::Sphere { axis, .. }) =
            out.faces[face].surface.as_mut().map(|s| &mut s.geometry)
        {
            *axis = pole;
        }
        let mut he = start;
        loop {
            out.half_edges[he].pcurve = None;
            he = out.half_edges[he].next;
            if he == start {
                break;
            }
        }
        out.face_use_pcurves.retain(|key, _| key.face != face);
    }
    out.recompute_pcurves_for_lenient(&caps);
    out
}

/// PK-872: the volume quadrature integrates each boundary chord's
/// sliver up to the trim, so the caps' volume no longer depends on how
/// the chart samples their rims (straight chords read -1.13e-5 of
/// Steiner with the pole perpendicular to each cap and +2.7e-6 with
/// +Z): every chart reads within 1e-7.
fn volume_does_not_depend_on_the_ball_chart(brep: &BRep, steiner: f64) {
    let perpendicular = |d: Vec3| {
        let other = if d.x.abs() < 0.9 {
            Vec3::new(1.0, 0.0, 0.0)
        } else {
            Vec3::new(0.0, 1.0, 0.0)
        };
        d.cross(&other).normalized().expect("perpendicular")
    };
    let charts: [(&str, ChartAxis); 4] = [
        ("+Z", Box::new(|_| Vec3::new(0.0, 0.0, 1.0))),
        ("+X", Box::new(|_| Vec3::new(1.0, 0.0, 0.0))),
        ("perpendicular to each cap", Box::new(perpendicular)),
        (
            "tilted 50 degrees off each cap",
            Box::new(move |d: Vec3| {
                let p = perpendicular(d);
                (d * 50f64.to_radians().cos() + p * 50f64.to_radians().sin())
                    .normalized()
                    .expect("tilted")
            }),
        ),
    ];
    let mut report = Vec::new();
    for (name, axis) in &charts {
        let body = recharted(brep, axis);
        let volume = popkern_measure::body_volume(&body, &body.tolerance).expect("volume");
        report.push((name, (volume - steiner) / steiner));
    }
    eprintln!("[pk872] chart -> relative volume error: {report:?}");
    for (name, error) in &report {
        assert!(
            error.abs() < 1e-7,
            "{name}: relative volume error {error:e}"
        );
    }
}
