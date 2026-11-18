//! A NURBS ball (a clamped net closed round its axis, both rows
//! collapsed to poles) meshes on the atlas, alone and as a cavity
//! (PK-1267). Its seam pair's second leg starts at a pole the point
//! inversion leaves on the far border (u = 1, the seed's): the chain lift
//! took it back to u = 0 and the ring ran down and up one border. Alone
//! the band declined and the ball fell back to the uniform grid; as a
//! box's cavity (its face reversed) the band meshed inside out, the body
//! 2.4e-3 of the ball's volume high (the generated corpus's PxB-N scene).

use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::TessellationOptions;

const BALL: &str = r#"fn ball_net(c, r, s) {
    let w = sqrt(0.5);
    let ci = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-1.0, 1.0], [-1.0, 0.0], [-1.0, -1.0], [0.0, -1.0], [1.0, -1.0], [1.0, 0.0]];
    let cw = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
    let mer = [[0.0, -r], [r, -r], [r, 0.0], [r, r], [0.0, r]];
    let mw = [1.0, w, 1.0, w, 1.0];
    let net = [];
    let wts = [];
    for i in 0..9 {
        let row = [];
        let wr = [];
        for j in 0..5 {
            row.push([c[0] + s[0] * mer[j][0] * ci[i][0], c[1] + s[1] * mer[j][0] * ci[i][1], c[2] + s[2] * mer[j][1]]);
            wr.push(cw[i] * mw[j]);
        }
        net.push(row);
        wts.push(wr);
    }
    #{ net: net, weights: wts }
}
fn ball(m, c, r, s) {
    let b = ball_net(c, r, s);
    m.sheet_nurbs(b.net, #{ degree: [2, 2], weights: b.weights,
        knots: [[0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]] })
}"#;

const CENTRE: [f64; 3] = [0.4, -0.3, 0.2];
const RADIUS: f64 = 1.1;
const SCALE: [f64; 3] = [1.0, 1.15, 0.95];

fn body(script: &str) -> BRep {
    let source = format!("{BALL}\nlet m = model();\n{script}\nm");
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host.run(&source, "nurbs_ball_atlas").expect("runs");
    run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone()
}

fn ellipsoid_volume() -> f64 {
    4.0 / 3.0 * std::f64::consts::PI * RADIUS.powi(3) * SCALE.iter().product::<f64>()
}

/// Meshed on the atlas (no face left to the grid), watertight, every
/// facet of the ball's face turned outward from its centre (`cavity`:
/// toward it); returns the mesh's enclosed volume.
fn assert_meshes(name: &str, body: &BRep, cavity: bool) -> f64 {
    let options = TessellationOptions::default();
    let mesh = popkern_mesh::tessellate(body, &options).expect("mesh");
    assert!(
        mesh.grid_faces.is_empty(),
        "{name}: {:?} fell back to the grid",
        mesh.grid_faces
    );
    assert_eq!(
        popkern_mesh::welded_open_edges(&mesh),
        0,
        "{name}: watertight"
    );
    let c = Point3::new(CENTRE[0], CENTRE[1], CENTRE[2]);
    let mut turned = 0.0;
    let mut volume = 0.0;
    for (t, face) in mesh.triangles.iter().zip(&mesh.face_map) {
        let [a, b, d] = t.map(|i| mesh.vertices[i].position);
        let n = (b - a).cross(&(d - a));
        let to = |p: Point3| Vec3::new(p.x, p.y, p.z);
        volume += to(a).dot(&to(b).cross(&to(d))) / 6.0;
        let g = Vec3::new(
            (a.x + b.x + d.x) / 3.0 - c.x,
            (a.y + b.y + d.y) / 3.0 - c.y,
            (a.z + b.z + d.z) / 3.0 - c.z,
        );
        let on_ball = body.faces[*face]
            .surface
            .as_ref()
            .is_some_and(|s| matches!(s.geometry, SurfaceGeometry::Nurbs(_)));
        if on_ball && (n.dot(&g) > 0.0) == cavity {
            turned += 0.5 * n.length();
        }
    }
    assert!(
        turned == 0.0,
        "{name}: facet area {turned} turned the wrong way"
    );
    volume
}

#[test]
fn a_nurbs_ball_meshes_on_the_atlas() {
    let [x, y, z] = CENTRE;
    let [a, b, c] = SCALE;
    let ball = body(&format!(
        "ball(m, [{x}, {y}, {z}], {RADIUS}, [{a}, {b}, {c}]);"
    ));
    let v = assert_meshes("ball", &ball, false);
    // The measure integrates the atlas's triangles on the exact surface
    // (its seam triangles anchored on u = 0 put their far corners past a
    // clamped domain and fell back to facets: 1.3e-5 low).
    let q = popkern_measure::body_volume(&ball, &Tolerance::default()).expect("volume");
    let want = ellipsoid_volume();
    assert!((q - want).abs() < 1e-9 * want, "volume {q} want {want}");
    // The chord deficit: under by at most twice the chord height over
    // the smallest semi-axis, relative.
    let chord = TessellationOptions::default().chord_height;
    let rel = (v - want) / want;
    assert!(
        (-3.0 * 2.0 * chord / (0.95 * RADIUS)..1e-9).contains(&rel),
        "mesh volume {v} for {want}"
    );
}

#[test]
fn a_nurbs_ball_cavity_meshes_on_the_atlas() {
    let [x, y, z] = CENTRE;
    let [a, b, c] = SCALE;
    let cavity = body(&format!(
        "let box = m.box(4.0, 4.0, 4.0, #{{ at: [-2.0, -2.0, -2.0] }});\n\
         box.subtract(m.tool(|m| ball(m, [{x}, {y}, {z}], {RADIUS}, [{a}, {b}, {c}])));"
    ));
    assert_meshes("cavity", &cavity, true);
    let v = popkern_measure::body_volume(&cavity, &Tolerance::default()).expect("volume");
    let want = 64.0 - ellipsoid_volume();
    assert!((v - want).abs() < 1e-6 * want, "volume {v} want {want}");
}
