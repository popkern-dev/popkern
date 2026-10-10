//! PK-295 regression pin: the served rim-filleted-cylinder minus
//! strip-crossing-cuboid subtraction must TESSELLATE WATERTIGHT.
//! Promoted from the PK-295 diagnostic probe — the historical result
//! shipped 5 faces and 20 open mesh edges with a clean validate_all,
//! and mid-campaign states shipped 57 open edges along the seam-side
//! bulge until the slit-seam dissolution pass landed.

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
use popkern_core::geometry::{Point3, Transform3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use std::collections::HashMap;

fn rim_filleted_cylinder(tol: Tolerance) -> BRep {
    let brep = BRep::cylinder(1.0, 2.0, tol).expect("cylinder");
    let mut rim = Vec::new();
    for (e, d) in brep.edges.iter() {
        let he = &brep.half_edges[d.half_edge];
        let a = brep.vertices[he.origin].position;
        let b = brep.vertices[brep.half_edges[he.twin].origin].position;
        let f1 = brep.loops[he.loop_].face;
        let f2 = brep.loops[brep.half_edges[he.twin].loop_].face;
        if f1 != f2 && (a.z - 2.0).abs() < 1e-9 && (b.z - 2.0).abs() < 1e-9 {
            rim.push(e);
        }
    }
    let specs: Vec<FilletSpec> = rim
        .iter()
        .map(|&edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(0.2),
        })
        .collect();
    fillet_edges(&brep, &specs, &FilletOptions::default(), tol, None)
        .expect("rim fillet must complete")
        .brep
}

fn cube_cutter(tol: Tolerance) -> BRep {
    BRep::unit_cube(tol)
        .expect("cube")
        .transform(&Transform3::scale_uniform(1.8))
        .expect("scale")
        .transform(&Transform3::translation(-0.9, -0.9, 0.6))
        .expect("translate")
}

/// Open mesh edges of `result`'s tessellation (welded at 1e-6).
fn open_mesh_edges(result: &BRep) -> usize {
    let mesh = popkern_mesh::tessellate(result, &popkern_mesh::TessellationOptions::default())
        .expect("tessellation");
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let mut canon: Vec<Point3> = Vec::new();
    let id = |p: Point3, key: &mut HashMap<(i64, i64, i64), usize>, canon: &mut Vec<Point3>| {
        let k = (
            (p.x * 1.0e6).round() as i64,
            (p.y * 1.0e6).round() as i64,
            (p.z * 1.0e6).round() as i64,
        );
        *key.entry(k).or_insert_with(|| {
            canon.push(p);
            canon.len() - 1
        })
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position, &mut key, &mut canon),
            id(mesh.vertices[t[1]].position, &mut key, &mut canon),
            id(mesh.vertices[t[2]].position, &mut key, &mut canon),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    counts.values().filter(|&&c| c == 1).count()
}

#[test]
fn served_subtraction_watertightness() {
    let tol = Tolerance::default();
    let result = popkern_bool::subtraction(rim_filleted_cylinder(tol), cube_cutter(tol), tol, None)
        .expect("difference");
    let open = open_mesh_edges(&result);
    assert_eq!(
        open, 0,
        "served subtraction must tessellate watertight (PK-295), {open} open mesh edges"
    );
}

/// The fillet torus band's boundary on the cylinder lies EXACTLY on the
/// torus's v = 0 tube seam (PK-774). A rim circle carrying centre
/// roundoff (a fitted 1.7999999999999987 for 1.8) put the boundary a
/// hair below the seam, the record points projected to v = 2π − ε, the
/// walk split the band's records across the seam and the subtraction
/// left a loop open by 0.87. The band's embeddings now fold into its own
/// circular v window, so every nudge serves the same watertight body.
#[test]
fn a_band_bounded_on_its_v_seam_serves_whatever_side_roundoff_puts_it() {
    let tol = Tolerance::default();
    let mut faces = Vec::new();
    for dz in [0.0, -1.3e-15, 1.0e-14, -1.0e-14, 3.0e-14] {
        let mut part = rim_filleted_cylinder(tol);
        let mut nudged = 0;
        for d in part.edges.values_mut() {
            if let Some(c) = d.curve.as_mut()
                && let popkern_core::geometry::curve::CurveGeometry::Arc { center, radius, .. } =
                    &mut c.geometry
                && (center.z - 1.8).abs() < 1e-9
                && (*radius - 1.0).abs() < 1e-9
            {
                center.z += dz;
                nudged += 1;
            }
        }
        assert_eq!(nudged, 1, "the band's cylinder-side rim circle");
        let result = popkern_bool::subtraction(part, cube_cutter(tol), tol, None)
            .unwrap_or_else(|e| panic!("dz {dz:e}: {e:?}"));
        assert_eq!(open_mesh_edges(&result), 0, "dz {dz:e}: open mesh edges");
        faces.push(result.faces.len());
    }
    assert!(
        faces.windows(2).all(|w| w[0] == w[1]),
        "face counts {faces:?}"
    );
}

/// Volume of a closed tessellation by the divergence theorem.
fn mesh_volume(result: &BRep) -> f64 {
    let mesh = popkern_mesh::tessellate(result, &popkern_mesh::TessellationOptions::default())
        .expect("tessellation");
    mesh.triangles
        .iter()
        .map(|t| {
            let [a, b, c] = t.map(|i| mesh.vertices[i].position);
            let (a, b, c) = (
                popkern_core::geometry::Vec3::new(a.x, a.y, a.z),
                popkern_core::geometry::Vec3::new(b.x, b.y, b.z),
                popkern_core::geometry::Vec3::new(c.x, c.y, c.z),
            );
            a.dot(&b.cross(&c)) / 6.0
        })
        .sum()
}

/// The same scene rigidly rotated must serve the same body (PK-775).
/// A rotation re-anchors a rotational chart: the cylinder's and the
/// torus band's u-seams move to the shift angle (u ∈ [−0.002, 6.281]
/// at 7°), and the torus×plane marcher may author a record on the band's
/// 2π v-branch. Each angle has failed differently: a record straddling
/// the band's moved seam met its rim only on another lift; a closed
/// record whose closure point lies on the rim saw one crossing; a 2π-
/// branch record found none. Each case dropped as Outside, the cellular
/// assembly went non-manifold and graph-merge left the 0.87 gap or
/// curveless edges.
#[test]
fn the_rigidly_rotated_scene_serves_the_same_body_at_every_angle() {
    let tol = Tolerance::default();
    let base = popkern_bool::subtraction(rim_filleted_cylinder(tol), cube_cutter(tol), tol, None)
        .expect("0°");
    let (faces0, volume0) = (base.faces.len(), mesh_volume(&base));
    for deg in [7.0_f64, 13.0, 30.0, 41.0, 77.0] {
        let t = Transform3::rotation_around_axis(
            Point3::new(0.0, 0.0, 0.0),
            popkern_core::geometry::Vec3::new(1.0, 0.3, 0.0),
            deg.to_radians(),
        )
        .expect("rotation");
        let part = rim_filleted_cylinder(tol).transform(&t).expect("part");
        let cutter = cube_cutter(tol).transform(&t).expect("cutter");
        let result = popkern_bool::subtraction(part, cutter, tol, None)
            .unwrap_or_else(|e| panic!("{deg}°: {e:?}"));
        assert_eq!(open_mesh_edges(&result), 0, "{deg}°: open mesh edges");
        assert_eq!(result.faces.len(), faces0, "{deg}°: face count");
        let v = mesh_volume(&result);
        assert!(
            (v - volume0).abs() / volume0 < 1e-3,
            "{deg}°: volume {v} vs {volume0} at 0°"
        );
    }
}
