//! Stage-D repeated-apex topology rewrite acceptance (PK-092, re-homed
//! from PK-005 item 8): under the strict V-free contract no production
//! face retains the old apex vertex in its outer loop, cap-to-cylinder
//! seams are true twin pairs (one cap side, one cylinder side), and the
//! production seam-ownership inventory is zero (no production side may
//! twin residual scaffolding or dangle unpaired).

use crate::direct_apex_validate_helpers::{
    EnvGuard, build_under_direct_apex, direct_apex_fallback_reasons,
};
use popkern_core::topology::TopologyRole;

#[test]
#[ignore = "diagnostic gate: uses process-wide strict direct-apex env and must run in isolation"]
fn direct_apex_stage_d_v_free_topology_acceptance() {
    let _guard = EnvGuard::set("POPKERN_DIRECT_APEX_REQUIRE_V_FREE_REPEATED_APEX", "1");
    let result = build_under_direct_apex(true, &[0.04, 0.06, 0.05, 0.07]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "Stage-D acceptance requires the native strict path; got fallbacks:\n{}",
        reasons.join("\n"),
    );
    let brep = &result.brep;

    // 1. No production face's outer loop visits the old apex (0,0,0.8).
    let mut apex_holders = Vec::new();
    for (face, face_data) in &brep.faces {
        if face_data.role == TopologyRole::ConstructionResidual {
            continue;
        }
        let Some(loop_data) = brep.loops.get(face_data.outer_loop) else {
            continue;
        };
        let start = loop_data.half_edge;
        let mut cursor = start;
        for _ in 0..=brep.half_edges.len() {
            let Some(half_edge) = brep.half_edges.get(cursor) else {
                break;
            };
            if let Some(vertex) = brep.vertices.get(half_edge.origin) {
                let p = vertex.position;
                if (p.x.powi(2) + p.y.powi(2)).sqrt() < 1.0e-9 && (p.z - 0.8).abs() < 1.0e-9 {
                    apex_holders.push(face);
                }
            }
            cursor = half_edge.next;
            if cursor == start {
                break;
            }
        }
    }
    assert!(
        apex_holders.is_empty(),
        "production faces retain the old apex vertex: {apex_holders:?}"
    );

    // 2 + 3. Production seam-ownership inventory: every half-edge owned
    // by a production face must twin a live half-edge that is itself
    // owned by a production face. Violations are the handover-completion
    // backlog (weld or retire); Stage-D demands zero.
    let mut violations = Vec::new();
    for (heh, half_edge) in &brep.half_edges {
        let Some(loop_data) = brep.loops.get(half_edge.loop_) else {
            continue;
        };
        let owner_is_production = brep
            .faces
            .get(loop_data.face)
            .is_some_and(|face| face.role != TopologyRole::ConstructionResidual);
        if !owner_is_production {
            continue;
        }
        let twin_owner = brep
            .half_edges
            .get(half_edge.twin)
            .and_then(|twin| brep.loops.get(twin.loop_))
            .map(|twin_loop| twin_loop.face);
        let twin_is_production = twin_owner.is_some_and(|face| {
            brep.faces
                .get(face)
                .is_some_and(|face_data| face_data.role != TopologyRole::ConstructionResidual)
        });
        if !twin_is_production {
            violations.push((heh, loop_data.face, twin_owner));
        }
    }
    assert!(
        violations.is_empty(),
        "production seam-ownership inventory must be zero; {} production side(s) twin residual \
         scaffolding or dangle unpaired:\n{}",
        violations.len(),
        violations
            .iter()
            .map(|(heh, face, twin_owner)| format!("  {heh:?} in {face:?} -> {twin_owner:?}"))
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
#[ignore = "diagnostic gate: uses process-wide strict direct-apex env and must run in isolation"]
fn direct_apex_strict_v_free_repeated_apex_succeeds_without_fallback() {
    let _guard = EnvGuard::set("POPKERN_DIRECT_APEX_REQUIRE_V_FREE_REPEATED_APEX", "1");
    let result = build_under_direct_apex(true, &[0.04, 0.06, 0.05, 0.07]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "strict V-free repeated-apex mode must complete without fallback; got:\n{}",
        reasons.join("\n"),
    );
    assert!(
        !reasons
            .iter()
            .any(|reason| reason.contains("repeated-apex construction retained apex vertex")),
        "strict V-free repeated-apex mode must not retain production apex topology; got:\n{}",
        reasons.join("\n"),
    );
}

/// PK-126 leg 3: pyramid3_uniform under strict V-free completes
/// natively since the residual-strip duplication falls back to the
/// boundary outcome's LIVE seam when the plan-recorded cylinder
/// candidate was consumed by the N=3 trim-rewrite ordering
/// (measured: plan error -> legacy g1 ~1.0 before; native
/// g0 2.2e-4 / g1 3.9e-4 after).
#[test]
#[ignore = "diagnostic gate: uses process-wide strict direct-apex env and must run in isolation"]
fn direct_apex_strict_v_free_pyramid3_uniform_succeeds_without_fallback() {
    let _guard = EnvGuard::set("POPKERN_DIRECT_APEX_REQUIRE_V_FREE_REPEATED_APEX", "1");
    let result = build_under_direct_apex(true, &[0.05, 0.05, 0.05]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "strict V-free pyramid3_uniform must complete without fallback; got:\n{}",
        reasons.join("\n"),
    );
}

/// PK-129 acceptance: cube_mixed under strict V-free completes
/// natively — the pocket junction is settled before boundary
/// materialization, so cap/strip/residual assemblies see no split
/// junctions.
#[test]
#[ignore = "diagnostic gate: uses process-wide strict direct-apex env and must run in isolation"]
fn direct_apex_strict_v_free_cube_mixed_succeeds_without_fallback() {
    let _guard = EnvGuard::set("POPKERN_DIRECT_APEX_REQUIRE_V_FREE_REPEATED_APEX", "1");
    let result = build_under_direct_apex(false, &[0.08, 0.10, 0.12]);
    let reasons = direct_apex_fallback_reasons(&result);
    assert!(
        reasons.is_empty(),
        "strict V-free cube_mixed must complete without fallback; got:\n{}",
        reasons.join("\n"),
    );
}
