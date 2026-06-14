use std::collections::HashMap;

use popkern_core::geometry::gregory::RectangularGregoryBiNonic;
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::super::BoundarySide;
use super::super::cell_reauthor::AuthoredSubChain;
use super::super::chain_discovery::ChainTrim;
use super::super::chain_fit::ChainRowKey;
use super::{joint_endpoint_leg_equations, joint_endpoint_plane_constraints, leg_vector_targets};

fn base_patch() -> RectangularGregoryBiNonic {
    let mut patch = RectangularGregoryBiNonic::zeros();
    for i in 0..10 {
        for j in 0..10 {
            patch.grid[i][j] = Point3::new(i as f64 / 9.0, j as f64 / 9.0, 0.0);
        }
    }
    patch
}

fn set_boundary_cross(patch: &mut RectangularGregoryBiNonic, side: BoundarySide, cross: Vec3) {
    let offset = cross * (1.0 / 9.0);
    match side {
        BoundarySide::UMin => {
            for k in 0..10 {
                patch.grid[1][k] = patch.grid[0][k] + offset;
            }
        }
        BoundarySide::UMax => {
            for k in 0..10 {
                patch.grid[8][k] = patch.grid[9][k] + offset;
            }
        }
        BoundarySide::VMin => {
            for k in 0..10 {
                patch.grid[k][1] = patch.grid[k][0] + offset;
            }
        }
        BoundarySide::VMax => {
            for k in 0..10 {
                patch.grid[k][8] = patch.grid[k][9] + offset;
            }
        }
    }
}

fn chain_with_endpoint_derivatives(start: Vec3, end: Vec3) -> [Point3; 10] {
    let mut cps = [Point3::new(0.0, 0.0, 0.0); 10];
    cps[1] = cps[0] + start * (1.0 / 9.0);
    cps[9] = Point3::new(1.0, 0.0, 0.0);
    cps[8] = cps[9] - end * (1.0 / 9.0);
    cps
}

#[test]
fn trim_endpoint_targets_chain_owned_cp1_to_interior_peer_vector() {
    let face: FaceHandle = FaceHandle::default();
    let start_vertex: VertexHandle = VertexHandle::default();
    let end_vertex: VertexHandle = VertexHandle::default();
    let row_a = ChainRowKey {
        face,
        side: BoundarySide::UMin,
    };
    let row_b = ChainRowKey {
        face,
        side: BoundarySide::UMax,
    };
    let owner_a = ChainRowKey {
        face,
        side: BoundarySide::VMax,
    };
    let mut patches = HashMap::new();
    let mut patch = base_patch();
    set_boundary_cross(&mut patch, BoundarySide::UMin, Vec3::new(0.0, 0.5, 0.0));
    set_boundary_cross(&mut patch, BoundarySide::UMax, Vec3::new(0.0, 2.0, 0.0));
    patches.insert(face, patch);
    let authored = vec![AuthoredSubChain {
        chain_cps: [Point3::new(0.0, 0.0, 0.0); 10],
        rows: vec![(owner_a, (0.0, 1.0))],
        nodes: Vec::new(),
    }];
    let trims = vec![ChainTrim {
        row_a,
        interval_a: (0.0, 1.0),
        row_b,
        // End node pairs row_a endpoint 1.0 with row_b interior
        // 0.5, the T-junction case that cannot be fixed by
        // copying row-1 interiors after the fact.
        interval_b: (0.5, 0.0),
        node_start: start_vertex,
        node_end: end_vertex,
    }];
    let (targets, count, max_delta) = leg_vector_targets(&authored, &trims, &patches);
    assert_eq!(count, 1, "only the chain-owned endpoint should be targeted");
    assert!(
        (max_delta - 1.5).abs() < 1.0e-12,
        "target delta should measure endpoint leg to interior peer"
    );
    let rows = targets
        .get(&0)
        .expect("owner chain should receive a target");
    assert_eq!(rows.len(), 1);
    let (station, target, weight) = rows[0];
    assert!((station - 0.0).abs() < 1.0e-12);
    assert!((weight - 5.0).abs() < 1.0e-12);
    assert!(
        (target - Vec3::new(0.0, 2.0, 0.0))
            .dot(&(target - Vec3::new(0.0, 2.0, 0.0)))
            .sqrt()
            < 1.0e-12,
        "endpoint owner should be retargeted to the interior side's cross derivative"
    );
}

#[test]
fn trim_endpoint_planes_chain_owned_leg_against_interior_peer_plane() {
    let face: FaceHandle = FaceHandle::default();
    let row_a = ChainRowKey {
        face,
        side: BoundarySide::UMin,
    };
    let row_b = ChainRowKey {
        face,
        side: BoundarySide::UMax,
    };
    let owner_a = ChainRowKey {
        face,
        side: BoundarySide::VMax,
    };
    let mut patch = base_patch();
    set_boundary_cross(&mut patch, BoundarySide::UMax, Vec3::new(1.0, 0.0, 0.0));
    let mut patches = HashMap::new();
    patches.insert(face, patch);
    let authored = vec![AuthoredSubChain {
        chain_cps: chain_with_endpoint_derivatives(
            Vec3::new(0.0, 0.0, 3.0),
            Vec3::new(0.0, 0.0, 0.0),
        ),
        rows: vec![(owner_a, (0.0, 1.0))],
        nodes: Vec::new(),
    }];
    let trims = vec![ChainTrim {
        row_a,
        interval_a: (0.0, 1.0),
        row_b,
        interval_b: (0.5, 0.0),
        node_start: VertexHandle::default(),
        node_end: VertexHandle::default(),
    }];
    let constraints = joint_endpoint_plane_constraints(&authored, &trims, &patches);
    assert_eq!(constraints.len(), 1);
    let constraint = constraints[0];
    assert_eq!(constraint.endpoint.chain_index, 0);
    assert!((constraint.endpoint.station - 0.0).abs() < 1.0e-12);
    assert!(
        (constraint.residual_before - 3.0).abs() < 1.0e-12,
        "endpoint leg's out-of-peer-plane component should be measured"
    );
}

#[test]
fn trim_endpoint_pairs_chain_owned_legs_for_joint_solve() {
    let face: FaceHandle = FaceHandle::default();
    let row_a = ChainRowKey {
        face,
        side: BoundarySide::UMin,
    };
    let row_b = ChainRowKey {
        face,
        side: BoundarySide::UMax,
    };
    let owner_a = ChainRowKey {
        face,
        side: BoundarySide::VMin,
    };
    let owner_b = ChainRowKey {
        face,
        side: BoundarySide::VMax,
    };
    let authored = vec![
        AuthoredSubChain {
            chain_cps: chain_with_endpoint_derivatives(
                Vec3::new(1.0, 0.0, 0.0),
                Vec3::new(0.0, 0.0, 0.0),
            ),
            rows: vec![(owner_a, (0.0, 1.0))],
            nodes: Vec::new(),
        },
        AuthoredSubChain {
            chain_cps: chain_with_endpoint_derivatives(
                Vec3::new(0.0, 0.0, 0.0),
                Vec3::new(-2.0, 0.0, 0.0),
            ),
            rows: vec![(owner_b, (0.0, 1.0))],
            nodes: Vec::new(),
        },
    ];
    let trims = vec![ChainTrim {
        row_a,
        interval_a: (0.0, 1.0),
        row_b,
        interval_b: (0.0, 1.0),
        node_start: VertexHandle::default(),
        node_end: VertexHandle::default(),
    }];
    let equations = joint_endpoint_leg_equations(&authored, &trims);
    assert_eq!(equations.len(), 2, "both trim endpoints are chain-owned");
    let first = equations
        .iter()
        .find(|equation| equation.first.station == 0.0)
        .expect("start endpoint equation should be present");
    assert_eq!(first.first.chain_index, 0);
    assert_eq!(first.second.chain_index, 1);
    assert!((first.first.derivative_scale - 1.0).abs() < 1.0e-12);
    assert!((first.second.derivative_scale + 1.0).abs() < 1.0e-12);
    assert!((first.residual_before - 1.0).abs() < 1.0e-12);
}
