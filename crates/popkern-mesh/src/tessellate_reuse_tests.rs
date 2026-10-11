//! The final face pass reuses the encroachment rounds' boundary-first
//! meshes (PK-794): same mesh, fewer face builds.

use popkern_core::progress::{ProgressEvent, ProgressReporter, StageTracker};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::{TessellationOptions, TriMesh};

#[derive(Default)]
struct Stages(Vec<ProgressEvent>);

impl ProgressReporter for Stages {
    fn on_stage(&mut self, event: ProgressEvent) {
        self.0.push(event);
    }
}

/// The PK-295 served subtraction: a rim-filleted cylinder minus a
/// strip-crossing cuboid, which drives two encroachment rounds.
fn encroaching_subtraction() -> BRep {
    use popkern_blend::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};
    let tol = Tolerance::default();
    let cyl = BRep::cylinder(1.0, 2.0, tol).unwrap();
    let rim: Vec<_> = cyl
        .edges
        .iter()
        .filter(|(_, d)| {
            let he = &cyl.half_edges[d.half_edge];
            let a = cyl.vertices[he.origin].position;
            let b = cyl.vertices[cyl.half_edges[he.twin].origin].position;
            let f1 = cyl.loops[he.loop_].face;
            let f2 = cyl.loops[cyl.half_edges[he.twin].loop_].face;
            f1 != f2 && (a.z - 2.0).abs() < 1e-9 && (b.z - 2.0).abs() < 1e-9
        })
        .map(|(e, _)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(0.2),
        })
        .collect();
    let part = fillet_edges(&cyl, &rim, &FilletOptions::default(), tol, None)
        .unwrap()
        .brep;
    let cutter = BRep::unit_cube(tol)
        .unwrap()
        .transform(&popkern_core::geometry::Transform3::scale_uniform(1.8))
        .unwrap()
        .transform(&popkern_core::geometry::Transform3::translation(
            -0.9, -0.9, 0.6,
        ))
        .unwrap();
    popkern_bool::subtraction(part, cutter, tol, None).unwrap()
}

fn mesh(brep: &BRep, reuse: bool) -> (TriMesh, Vec<ProgressEvent>) {
    let mut stages = Stages::default();
    let mesh = {
        let mut tracker = StageTracker::new(None, Some(&mut stages), 64);
        super::tessellate_brep_reusing(brep, &TessellationOptions::default(), &mut tracker, reuse)
            .unwrap()
    };
    (mesh, stages.0)
}

fn count(events: &[ProgressEvent], stage: &str) -> Vec<usize> {
    events
        .iter()
        .filter(|e| e.stage == stage)
        .map(|e| e.entity_count)
        .collect()
}

#[test]
fn reused_round_meshes_match_a_full_rebuild_bit_for_bit() {
    let brep = encroaching_subtraction();
    let faces = brep.faces.values().filter(|f| f.surface.is_some()).count();
    let (reused, reused_stages) = mesh(&brep, true);
    let (rebuilt, rebuilt_stages) = mesh(&brep, false);

    // The fixture splits edges, so faces were checked more than once
    // and some checks went stale.
    let splits = count(&reused_stages, "tessellate.encroachment_round");
    assert!(
        splits.iter().any(|&n| n > 0),
        "no encroachment splits: {splits:?}"
    );
    // Every face's last check ran on the final samples, so nothing is
    // rebuilt; the reference rebuilds every face.
    assert_eq!(count(&reused_stages, "tessellate.face_pass"), vec![0]);
    assert_eq!(count(&rebuilt_stages, "tessellate.face_pass"), vec![faces]);

    assert_eq!(reused.triangles, rebuilt.triangles);
    assert_eq!(reused.face_map, rebuilt.face_map);
    assert_eq!(reused.vertices.len(), rebuilt.vertices.len());
    for (a, b) in reused.vertices.iter().zip(&rebuilt.vertices) {
        let bits =
            |p: &popkern_core::geometry::Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
        assert_eq!(bits(&a.position), bits(&b.position));
        let nbits =
            |n: &popkern_core::geometry::Vec3| (n.x.to_bits(), n.y.to_bits(), n.z.to_bits());
        assert_eq!(nbits(&a.normal), nbits(&b.normal));
    }
}
