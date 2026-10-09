//! Instanced rendering (PK-759): a body placed many times is prepared
//! once and placed per occurrence; the picture equals rendering the
//! placed copies one by one.

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_viz::cpu_render::{Instance, prepare_instances};

fn translate(x: f64, y: f64, z: f64) -> [[f64; 4]; 4] {
    [
        [1.0, 0.0, 0.0, x],
        [0.0, 1.0, 0.0, y],
        [0.0, 0.0, 1.0, z],
        [0.0, 0.0, 0.0, 1.0],
    ]
}

#[test]
fn instances_of_one_body_are_its_mesh_placed_per_occurrence() {
    let cube = BRep::unit_cube(Tolerance::default()).expect("cube");
    let mesh = popkern_mesh::tessellate(&cube, &popkern_mesh::TessellationOptions::display())
        .expect("mesh");
    let one = prepare_instances(
        &[Instance {
            brep: &cube,
            mesh: &mesh,
            placement: translate(0.0, 0.0, 0.0),
        }],
        true,
    );
    // A quarter turn about z, then 5 along x.
    let turn = [
        [0.0, -1.0, 0.0, 5.0],
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ];
    let two = prepare_instances(
        &[
            Instance {
                brep: &cube,
                mesh: &mesh,
                placement: translate(0.0, 0.0, 0.0),
            },
            Instance {
                brep: &cube,
                mesh: &mesh,
                placement: turn,
            },
        ],
        true,
    );
    assert_eq!(two.0.len(), 2 * one.0.len());
    assert_eq!(two.1.len(), 2 * one.1.len());
    assert_eq!(two.2.len(), 2 * one.2.len());
    let n = one.0.len();
    for (k, v) in one.0.iter().enumerate() {
        let w = &two.0[n + k];
        let [x, y, z] = v.position;
        let want = [5.0 - y, x, z];
        for (got, want_c) in w.position.iter().zip(want) {
            assert!(
                (got - want_c).abs() < 1e-5,
                "vertex {k}: {:?} vs {want:?}",
                w.position
            );
        }
        // Normals turn with the body, unscaled.
        let [nx, ny, nz] = v.normal;
        let wn = [-ny, nx, nz];
        for (got, want_c) in w.normal.iter().zip(wn) {
            assert!((got - want_c).abs() < 1e-5);
        }
        assert_eq!(w.color, v.color);
    }
    // The second instance's triangles index its own vertices.
    assert!(two.1[one.1.len()..].iter().all(|&i| i as usize >= n));
}
