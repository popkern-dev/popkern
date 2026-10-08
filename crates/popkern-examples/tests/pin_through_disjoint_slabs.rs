//! PK-705: a cylinder unioned through several disjoint slabs of the
//! other operand — a comb (teeth joined by a back plate) with a pin
//! running along the teeth, through the gaps between them. The
//! print-in-place hinge's integral pin is this boolean: one pin
//! through every knuckle of leaf A, exposed in leaf B's zones.

use popkern_bool::merge_fallback::take_merge_fallback_events;
use popkern_bool::{BoolOperation, boolean};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_sketch::{AnalyticProfile, ProfileShape, SketchPlane};

const PIN_RADIUS: f64 = 2.0;
/// Tooth spans along z.
type Teeth = [(f64, f64)];
/// Teeth span x in [-TOOTH, TOOTH], y in [-TOOTH, TOOTH]; the plate
/// runs back to x = -PLATE.
const TOOTH: f64 = 3.0;
const PLATE: f64 = 15.0;

fn tol() -> Tolerance {
    Tolerance::default()
}

fn extrude(shape: ProfileShape<f64>, plane: SketchPlane, depth: f64) -> BRep {
    let sheet = AnalyticProfile {
        shape,
        center: [0.0, 0.0],
    }
    .to_brep(plane, tol(), None)
    .expect("sheet");
    let face = sheet.regions[0].face;
    let mut b = sheet.brep;
    let normal = plane.normal;
    b.extrude(face, normal, depth, None).expect("extrude");
    b
}

/// The comb: teeth at the given z spans joined by a plate behind them,
/// sketched in the xz plane and extruded across y in [-TOOTH, TOOTH].
fn comb(teeth: &[(f64, f64)]) -> BRep {
    let first = teeth[0].0;
    let last = teeth[teeth.len() - 1].1;
    let mut points = vec![[-PLATE, first], [TOOTH, first]];
    for (j, &(_, b)) in teeth.iter().enumerate() {
        points.push([TOOTH, b]);
        if let Some(&(next, _)) = teeth.get(j + 1) {
            points.push([-TOOTH, b]);
            points.push([-TOOTH, next]);
            points.push([TOOTH, next]);
        }
    }
    points.push([-PLATE, last]);
    // u = x, v = z: the normal x × z is -y, so start at y = +TOOTH.
    let plane = SketchPlane {
        origin: Point3::new(0.0, TOOTH, 0.0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 0.0, 1.0),
        normal: Vec3::new(0.0, -1.0, 0.0),
    };
    extrude(ProfileShape::Polygon { points }, plane, 2.0 * TOOTH)
}

fn pin(z0: f64, z1: f64) -> BRep {
    let plane = SketchPlane {
        origin: Point3::new(0.0, 0.0, z0),
        u_axis: Vec3::new(1.0, 0.0, 0.0),
        v_axis: Vec3::new(0.0, 1.0, 0.0),
        normal: Vec3::new(0.0, 0.0, 1.0),
    };
    extrude(ProfileShape::Circle { radius: PIN_RADIUS }, plane, z1 - z0)
}

/// Divergence-theorem volume of the display mesh, and its welded open
/// edge count.
fn mesh_volume(body: &BRep) -> (f64, usize) {
    let mesh = popkern_mesh::tessellate(body, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    let mut volume = 0.0;
    for t in &mesh.triangles {
        let [a, b, c] = t.map(|i| mesh.vertices[i].position);
        volume += (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
            + a.z * (b.x * c.y - b.y * c.x))
            / 6.0;
    }
    (volume, popkern_mesh::welded_open_edges(&mesh))
}

/// Comb volume plus the pin's length outside every tooth.
fn analytic_volume(teeth: &[(f64, f64)], z0: f64, z1: f64) -> f64 {
    let first = teeth[0].0;
    let last = teeth[teeth.len() - 1].1;
    let tooth_len: f64 = teeth.iter().map(|(a, b)| b - a).sum();
    let profile = (PLATE - TOOTH) * (last - first) + 2.0 * TOOTH * tooth_len;
    let inside: f64 = teeth
        .iter()
        .map(|&(a, b)| (b.min(z1) - a.max(z0)).max(0.0))
        .sum();
    let exposed = (z1 - z0) - inside;
    profile * 2.0 * TOOTH + std::f64::consts::PI * PIN_RADIUS * PIN_RADIUS * exposed
}

#[test]
fn a_pin_through_several_disjoint_teeth_unions_watertight() {
    // The hinge repro: five zones on L = 60, leaf A owns teeth 0, 2, 4;
    // the pin runs from the middle of the first tooth to the middle of
    // the last — two exposed gaps, three buried stretches.
    let hinge = [(0.0, 11.88), (24.06, 35.94), (48.12, 60.0)];
    let two = [(0.0, 10.0), (20.0, 30.0)];
    let four = [(0.0, 8.0), (14.0, 22.0), (28.0, 36.0), (42.0, 50.0)];
    let cases: [(&Teeth, f64, f64); 6] = [
        (&hinge, 5.94, 54.06),
        (&two, 5.0, 25.0),
        // Exposed end: the pin leaves the last tooth.
        (&two, 5.0, 36.0),
        (&hinge, 5.94, 66.0),
        (&four, 4.0, 46.0),
        (&four, -6.0, 46.0),
    ];
    for (teeth, z0, z1) in cases {
        let label = format!("{} teeth, pin z [{z0}, {z1}]", teeth.len());
        let _ = take_merge_fallback_events();
        let result = boolean(BoolOperation::Union, comb(teeth), pin(z0, z1), tol(), None)
            .unwrap_or_else(|e| panic!("{label}: {e}"));
        result
            .validate()
            .unwrap_or_else(|e| panic!("{label}: {e:?}"));
        let (volume, open) = mesh_volume(&result);
        assert_eq!(open, 0, "{label}: welded display mesh has open edges");
        let expected = analytic_volume(teeth, z0, z1);
        assert!(
            (volume - expected).abs() < 2e-3 * expected,
            "{label}: volume {volume:.3} vs analytic {expected:.3}"
        );
    }
}
