//! Operand kinds (design §3): one per surface kind, three parameter
//! variants each, built once per process and moved rigidly per scene.

use std::collections::HashMap;
use std::sync::{Arc, Mutex, OnceLock};

use popkern_core::geometry::Point3;
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use serde::{Deserialize, Serialize};

/// The operand kinds (design §3).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Kind {
    /// Plane: a box.
    P,
    /// Cylinder.
    Y,
    /// Cone with its apex.
    C,
    /// Cone frustum (no apex), revolved.
    F,
    /// Sphere (chart poles).
    S,
    /// Ring torus.
    T,
    /// NURBS of revolution: a spline profile revolved.
    R,
    /// Closed NURBS tube with flat discs (`sheet_nurbs` closed round u).
    Q,
    /// NURBS ball (both v rows collapsed).
    B,
    /// Freeform NURBS: a box split by a bicubic dome.
    W,
    /// Analytic blends: a box with every edge filleted.
    L,
    /// G2 blends: a box with G2 rounds (Gregory corners, NURBS strips).
    G,
}

/// Every solid kind, in table order.
pub const KINDS: [Kind; 12] = [
    Kind::P,
    Kind::Y,
    Kind::C,
    Kind::F,
    Kind::S,
    Kind::T,
    Kind::R,
    Kind::Q,
    Kind::B,
    Kind::W,
    Kind::L,
    Kind::G,
];

/// A built operand with the data the generator reads.
pub struct Operand {
    pub kind: Kind,
    pub variant: u8,
    pub brep: BRep,
    /// An outward-wound tessellation (signed volume positive).
    pub mesh: popkern_mesh::TriMesh,
    /// Triangles of the faces whose surface is the kind's.
    pub relevant: Vec<usize>,
    /// Cumulative areas of `relevant` (area-weighted sampling: a cone's
    /// apex fan of small triangles must not be over-sampled).
    pub relevant_area: Vec<f64>,
    /// B-rep quadrature volume.
    pub volume: f64,
    /// Tessellation divergence volume.
    pub mesh_volume: f64,
    /// Surface area.
    pub area: f64,
    /// Bounding-box diagonal.
    pub size: f64,
    /// The mesh's chord bound.
    pub chord: f64,
}

/// Whether a face surface belongs to the kind's relevant set.
pub fn relevant_surface(kind: Kind, s: &SurfaceGeometry) -> bool {
    use SurfaceGeometry as G;
    match kind {
        Kind::P => matches!(s, G::Plane { .. }),
        Kind::Y => matches!(s, G::Cylinder { .. }),
        Kind::C | Kind::F => matches!(s, G::Cone { .. }),
        Kind::S => matches!(s, G::Sphere { .. }),
        Kind::T => matches!(s, G::Torus { .. }),
        Kind::R | Kind::Q | Kind::B | Kind::W => matches!(s, G::Nurbs(_)),
        Kind::L => !matches!(s, G::Plane { .. }),
        Kind::G => matches!(s, G::Nurbs(_) | G::Gregory(_)),
    }
}

fn script(source: &str, label: &str) -> BRep {
    let mut host = popkern_rhai::ScriptHost::with_options(Default::default());
    let run = host
        .run(source, label)
        .unwrap_or_else(|e| panic!("operand {label}: {e}"));
    run.model.expect("model").to_parts_all().expect("parts")[0]
        .brep()
        .clone()
}

/// The rational biquadratic net of a ball of radius `r` scaled by `s`
/// (Rhai), as the sheet_nurbs tests build it.
const BALL_NET: &str = r#"fn ball_net(r, s) {
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
            row.push([s[0] * mer[j][0] * ci[i][0], s[1] * mer[j][0] * ci[i][1], s[2] * mer[j][1]]);
            wr.push(cw[i] * mw[j]);
        }
        net.push(row);
        wts.push(wr);
    }
    #{ net: net, weights: wts }
}
"#;

/// A closed tube net: the rational circle (u, 9 rows) at heights with
/// radii `radii` (v), degree 2 along both.
fn tube_source(radii: &[f64], height: f64) -> String {
    let n = radii.len();
    let rs: Vec<String> = radii.iter().map(|r| format!("{r:.6}")).collect();
    let mut vknots = vec![0.0, 0.0, 0.0];
    for k in 1..n - 2 {
        vknots.push(k as f64 / (n - 2) as f64);
    }
    vknots.extend([1.0, 1.0, 1.0]);
    let vk: Vec<String> = vknots.iter().map(|k| format!("{k:.6}")).collect();
    format!(
        r#"let m = model();
let w = sqrt(0.5);
let ci = [[1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [-1.0, 1.0], [-1.0, 0.0], [-1.0, -1.0], [0.0, -1.0], [1.0, -1.0], [1.0, 0.0]];
let cw = [1.0, w, 1.0, w, 1.0, w, 1.0, w, 1.0];
let radii = [{rs}];
let net = [];
let wts = [];
for i in 0..9 {{
    let row = [];
    let wr = [];
    for j in 0..{n} {{
        row.push([radii[j] * ci[i][0], radii[j] * ci[i][1], {height} * j / {last}.0]);
        wr.push(cw[i]);
    }}
    net.push(row);
    wts.push(wr);
}}
m.sheet_nurbs(net, #{{ degree: [2, 2], weights: wts,
    knots: [[0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0], [{vk}]] }});
m"#,
        rs = rs.join(", "),
        n = n,
        height = height,
        last = n - 1,
        vk = vk.join(", ")
    )
}

fn build(kind: Kind, variant: u8, tol: Tolerance) -> BRep {
    let k = variant as f64;
    let moved = |b: BRep, z: f64| {
        b.transform(&popkern_core::geometry::Transform3::translation(
            0.0, 0.0, z,
        ))
        .expect("move")
    };
    match kind {
        Kind::P => {
            let (x, y, z) = (2.0 + 0.5 * k, 1.6 + 0.3 * k, 1.2 + 0.4 * k);
            let rect = |h: f64| {
                vec![
                    Point3::new(-x / 2.0, -y / 2.0, h),
                    Point3::new(x / 2.0, -y / 2.0, h),
                    Point3::new(x / 2.0, y / 2.0, h),
                    Point3::new(-x / 2.0, y / 2.0, h),
                ]
            };
            BRep::polygon_frustum(&rect(-z / 2.0), &rect(z / 2.0), tol).expect("box")
        }
        Kind::Y => {
            let (r, h) = (0.8 + 0.2 * k, 2.0 + 0.5 * k);
            moved(BRep::cylinder(r, h, tol).expect("cylinder"), -h / 2.0)
        }
        Kind::C => {
            let (r, h) = (0.9 + 0.2 * k, 1.8 + 0.4 * k);
            moved(BRep::cone(r, h, tol).expect("cone"), -h / 2.0)
        }
        Kind::S => BRep::sphere(0.9 + 0.25 * k, tol).expect("sphere"),
        Kind::T => BRep::torus(1.2 + 0.2 * k, 0.4 + 0.1 * k, tol).expect("torus"),
        Kind::F => {
            let (r1, r2, h) = (1.0 + 0.2 * k, 0.5 + 0.1 * k, 1.6 + 0.3 * k);
            let b = script(
                &format!(
                    "let m = model();\nm.profile(xz()).path(0.0, 0.0).line_to({r1}, 0.0).line_to({r2}, {h})\
                     .line_to(0.0, {h}).close().sketch(\"frustum\")\
                     .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());\nm"
                ),
                "frustum",
            );
            moved(b, -h / 2.0)
        }
        Kind::R => {
            let h = 2.0 + 0.4 * k;
            let b = script(
                &format!(
                    "let m = model();\nm.profile(xz()).path(0.0, 0.0).line_to(0.9, 0.0)\
                     .spline_to([[{a}, {h1}], [{b}, {h2}], [0.6, {h}]], #{{}})\
                     .line_to(0.0, {h}).close().sketch(\"vase\")\
                     .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());\nm",
                    a = 1.25 + 0.1 * k,
                    b = 0.7 + 0.05 * k,
                    h1 = 0.35 * h,
                    h2 = 0.7 * h,
                ),
                "vase",
            );
            moved(b, -h / 2.0)
        }
        Kind::Q => {
            let h = 2.0 + 0.4 * k;
            let radii = [0.8, 1.05 + 0.1 * k, 0.75, 0.95, 0.85];
            moved(script(&tube_source(&radii, h), "tube"), -h / 2.0)
        }
        Kind::B => {
            let r = 0.9 + 0.2 * k;
            let s = [1.0, 1.0 + 0.15 * k, 0.85 + 0.1 * k];
            script(
                &format!(
                    "{BALL_NET}\nlet m = model();\nlet b = ball_net({r}, [{}, {}, {}]);\n\
                     m.sheet_nurbs(b.net, #{{ degree: [2, 2], weights: b.weights,\n\
                     knots: [[0.0, 0.0, 0.0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0], [0.0, 0.0, 0.0, 0.5, 0.5, 1.0, 1.0, 1.0]] }});\nm",
                    s[0], s[1], s[2]
                ),
                "nurbs_ball",
            )
        }
        Kind::W => {
            let a = 2.0 + 0.3 * k;
            let b = script(
                &format!(
                    r#"let m = model();
let z = [[1.2, 1.5, 1.6, 1.3], [1.4, 2.2, 2.0, 1.5], [1.3, 2.1, 2.4, 1.6], [1.1, 1.4, 1.6, 1.4]];
let net = [];
for i in 0..4 {{ let row = []; for j in 0..4 {{ row.push([{a} * i / 3.0 - {a} / 2.0, {a} * j / 3.0 - {a} / 2.0, z[i][j] * {s}]); }} net.push(row); }}
let dome = m.sheet_nurbs(net, #{{ degree: 3 }});
let b = m.box({a}, {a}, 4.0, #{{ at: [-{a} / 2.0, -{a} / 2.0, 0.0] }}).split([dome], containing([0.0, 0.0, 0.2]))[0];
m"#,
                    s = 1.0 + 0.1 * k
                ),
                "dome_block",
            );
            moved(b, -1.0)
        }
        Kind::L | Kind::G => {
            let (x, y, z, r) = (2.0 + 0.3 * k, 1.6 + 0.2 * k, 1.3 + 0.2 * k, 0.25 + 0.05 * k);
            let g2 = if kind == Kind::G {
                ", #{ continuity: \"g2\" }"
            } else {
                ""
            };
            script(
                &format!(
                    "let m = model();\nlet b = m.box({x}, {y}, {z}, #{{ at: [-{x} / 2.0, -{y} / 2.0, -{z} / 2.0] }});\n\
                     b.fillet(edges(), {r}{g2});\nm"
                ),
                if kind == Kind::G {
                    "g2_block"
                } else {
                    "round_block"
                },
            )
        }
    }
}

fn outward_mesh(brep: &BRep, chord: f64) -> popkern_mesh::TriMesh {
    let options = popkern_mesh::TessellationOptions {
        chord_height: chord,
        ..Default::default()
    };
    let mut mesh = popkern_mesh::tessellate(brep, &options).expect("tessellate operand");
    if super::geometry::mesh_signed_volume(&mesh) < 0.0 {
        for t in &mut mesh.triangles {
            t.swap(1, 2);
        }
    }
    mesh
}

fn assemble(kind: Kind, variant: u8) -> Operand {
    let tol = Tolerance::default();
    let brep = build(kind, variant, tol);
    let errors = brep.validate_all();
    assert!(
        errors.is_empty(),
        "operand {kind:?}/{variant}: {:?}",
        errors.first()
    );
    // The size from a coarse tessellation: `body_bounding_box` reads
    // vertices only (PK-1258) and misses every curved face's extent.
    let rough = popkern_measure::body_bounding_box(&brep)
        .map(|b| b.min.distance_to(&b.max))
        .unwrap_or(1.0)
        .max(1e-9);
    let coarse = outward_mesh(&brep, 1e-2 * rough);
    let (lo, hi) = super::distance::mesh_box(&coarse);
    let size = lo.distance_to(&hi);
    let chord = 2e-4 * size;
    let mesh = outward_mesh(&brep, chord);
    let relevant = mesh
        .face_map
        .iter()
        .enumerate()
        .filter(|(_, f)| {
            brep.faces
                .get(**f)
                .and_then(|face| face.surface.as_ref())
                .is_some_and(|s| relevant_surface(kind, &s.geometry))
        })
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    assert!(
        !relevant.is_empty(),
        "operand {kind:?}/{variant}: no relevant face"
    );
    let mut acc = 0.0;
    let relevant_area: Vec<f64> = relevant
        .iter()
        .map(|&t| {
            let tri = mesh.triangles[t];
            let p = |k: usize| mesh.vertices[tri[k]].position;
            acc += 0.5 * (p(1) - p(0)).cross(&(p(2) - p(0))).length();
            acc
        })
        .collect();
    let volume = popkern_measure::body_volume(&brep, &tol).expect("volume");
    let area = popkern_measure::body_surface_area(&brep, &tol).expect("area");
    let mesh_volume = super::geometry::mesh_signed_volume(&mesh);
    Operand {
        kind,
        variant,
        brep,
        mesh,
        relevant,
        relevant_area,
        volume,
        mesh_volume,
        area,
        size,
        chord,
    }
}

type Cache = Mutex<HashMap<(Kind, u8), Arc<Operand>>>;

/// The operand of `kind` and `variant` (0..3), built once per process.
pub fn operand(kind: Kind, variant: u8) -> Arc<Operand> {
    static CACHE: OnceLock<Cache> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    if let Some(op) = cache.lock().expect("cache").get(&(kind, variant)) {
        return op.clone();
    }
    let built = Arc::new(assemble(kind, variant % 3));
    cache
        .lock()
        .expect("cache")
        .entry((kind, variant))
        .or_insert(built)
        .clone()
}
