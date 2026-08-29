use super::*;

impl AdjacentLoopRewriteKind {
    pub(super) fn diagnostic_label(self) -> &'static str {
        match self {
            AdjacentLoopRewriteKind::ReverseAdjacentTwinCycle => "ReverseAdjacentTwinCycle",
            AdjacentLoopRewriteKind::BoundaryPreservingAdjacentLoopRepairRequired => {
                "BoundaryBearingAdjacentLoopRepairRequired"
            }
            AdjacentLoopRewriteKind::ConstructionResidualAdjacentLoopRepairRequired => {
                "ConstructionResidualAdjacentLoopRepairRequired"
            }
            AdjacentLoopRewriteKind::AttachTwinToLaterStripLoop => "AttachTwinToLaterStripLoop",
        }
    }

    pub(super) fn reason(self) -> &'static str {
        match self {
            AdjacentLoopRewriteKind::ReverseAdjacentTwinCycle => {
                "live adjacent loop must be reauthored because replacing a stolen half-edge with \
                 its twin reverses local loop direction"
            }
            AdjacentLoopRewriteKind::BoundaryPreservingAdjacentLoopRepairRequired => {
                "live adjacent loop has boundary sides; whole-loop twin reversal would use null \
                 twins, so the adjacent loop must receive a same-direction boundary replacement \
                 span before the cap steals the live half-edge"
            }
            AdjacentLoopRewriteKind::ConstructionResidualAdjacentLoopRepairRequired => {
                "adjacent loop belongs to construction-residual retirement topology; it must \
                 receive a same-direction boundary replacement span before the cap steals the \
                 live half-edge"
            }
            AdjacentLoopRewriteKind::AttachTwinToLaterStripLoop => {
                "new FilletEnd twin is not yet assigned to a live adjacent loop; caller must \
                 attach its twin to the later fillet strip loop"
            }
        }
    }
}

pub(super) fn author_boundary_preserving_adjacent_loop_repair(
    brep: &mut BRep,
    rewrite: &AdjacentLoopRewriteRequirement,
    loop_handle: LoopHandle,
    face: FaceHandle,
) -> Result<(), ApexBuildError> {
    let stolen_index = rewrite
        .current_neighbor_cycle
        .iter()
        .position(|half_edge| *half_edge == rewrite.stolen_cap_side)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot repair boundary-bearing adjacent loop \
                 {loop_handle:?}: stolen cap side {:?} is absent from current cycle {:?}",
                rewrite.stolen_cap_side, rewrite.current_neighbor_cycle
            ),
        })?;
    let next =
        rewrite.current_neighbor_cycle[(stolen_index + 1) % rewrite.current_neighbor_cycle.len()];
    let start = brep
        .half_edges
        .get(rewrite.stolen_cap_side)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot repair boundary-bearing adjacent loop \
                 {loop_handle:?}: stale stolen side {:?}",
                rewrite.stolen_cap_side
            ),
        })?;
    let expected_end = brep
        .half_edges
        .get(next)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot repair boundary-bearing adjacent loop \
                 {loop_handle:?}: stale successor {next:?}"
            ),
        })?;
    let curve = brep
        .half_edges
        .get(rewrite.stolen_cap_side)
        .and_then(|half_edge| brep.edges.get(half_edge.edge))
        .and_then(|edge| edge.curve.clone())
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot repair boundary-bearing adjacent loop \
                 {loop_handle:?}: stolen side {:?} has no curve to clone for the replacement span",
                rewrite.stolen_cap_side
            ),
        })?;
    let (edge, replacement) = insert_boundary_half_edge(brep, start);
    brep.edges[edge].curve = Some(curve);
    brep.edges[edge].role = TopologyRole::Production;
    if support_split_half_edge_end_vertex(brep, replacement) != Some(expected_end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot repair boundary-bearing adjacent loop \
                 {loop_handle:?}: authored replacement {replacement:?} does not span \
                 {start:?}->{expected_end:?}"
            ),
        });
    }
    let mut repaired_cycle = rewrite.current_neighbor_cycle.clone();
    repaired_cycle[stolen_index] = replacement;
    validate_half_edge_cycle(
        brep,
        &repaired_cycle,
        "support boundary-preserving adjacent-loop repair",
        face,
    )?;
    wire_existing_loop(brep, loop_handle, face, &repaired_cycle)?;
    validate_boundary_preserving_adjacent_loop_repair(
        brep,
        rewrite,
        loop_handle,
        replacement,
        &repaired_cycle,
    )
}

pub(super) fn clone_residual_owned_sides_for_adjacent_rewrite(
    brep: &mut BRep,
    cap_loop: LoopHandle,
    target_loop: LoopHandle,
    reauthored_cycle: &[HalfEdgeHandle],
) -> Result<Vec<HalfEdgeHandle>, ApexBuildError> {
    let mut protected_cycle = Vec::with_capacity(reauthored_cycle.len());
    for (side_index, source) in reauthored_cycle.iter().copied().enumerate() {
        if !half_edge_belongs_to_third_live_loop(brep, source, cap_loop, target_loop) {
            protected_cycle.push(source);
            continue;
        }
        let next = reauthored_cycle[(side_index + 1) % reauthored_cycle.len()];
        protected_cycle.push(clone_same_direction_boundary_side(
            brep, source, next, side_index,
        )?);
    }
    Ok(protected_cycle)
}

pub(super) fn protect_source_loops_before_new_loop_claim(
    brep: &mut BRep,
    claiming_cycle: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    let mut claims_by_loop: std::collections::BTreeMap<LoopHandle, Vec<HalfEdgeHandle>> =
        std::collections::BTreeMap::new();
    for half_edge in claiming_cycle {
        let Some(loop_handle) = brep.half_edges.get(*half_edge).map(|data| data.loop_) else {
            continue;
        };
        if loop_handle != LoopHandle::default() && brep.loops.contains_key(loop_handle) {
            claims_by_loop
                .entry(loop_handle)
                .or_default()
                .push(*half_edge);
        }
    }
    for (loop_handle, claimed_sides) in claims_by_loop {
        let face = brep
            .loops
            .get(loop_handle)
            .map(|loop_data| loop_data.face)
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot protect source loop {loop_handle:?}: \
                     loop disappeared before strip claim"
                ),
            })?;
        let mut repaired_cycle = collect_loop_half_edges(brep, loop_handle)?;
        for claimed_side in claimed_sides {
            let Some(side_index) = repaired_cycle
                .iter()
                .position(|half_edge| *half_edge == claimed_side)
            else {
                continue;
            };
            let next = repaired_cycle[(side_index + 1) % repaired_cycle.len()];
            repaired_cycle[side_index] =
                clone_same_direction_boundary_side(brep, claimed_side, next, side_index)?;
        }
        validate_half_edge_cycle(
            brep,
            &repaired_cycle,
            "support source-loop protection before strip claim",
            face,
        )?;
        wire_existing_loop(brep, loop_handle, face, &repaired_cycle)?;
        if brep
            .faces
            .get(face)
            .is_some_and(|face_data| face_data.role == TopologyRole::ConstructionResidual)
        {
            for half_edge in &repaired_cycle {
                if let Some(edge) = brep.half_edges.get(*half_edge).map(|data| data.edge)
                    && let Some(edge_data) = brep.edges.get_mut(edge)
                {
                    edge_data.role = TopologyRole::ConstructionResidual;
                }
            }
        }
    }
    Ok(())
}

fn clone_same_direction_boundary_side(
    brep: &mut BRep,
    source: HalfEdgeHandle,
    next: HalfEdgeHandle,
    side_index: usize,
) -> Result<HalfEdgeHandle, ApexBuildError> {
    let source_data =
        brep.half_edges
            .get(source)
            .cloned()
            .ok_or_else(|| ApexBuildError::PlanMismatch {
                reason: format!(
                    "planned repeated-apex partition cannot clone residual-owned adjacent rewrite \
                 side[{side_index}]: stale source {source:?}"
                ),
            })?;
    let expected_end = brep
        .half_edges
        .get(next)
        .map(|half_edge| half_edge.origin)
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot clone residual-owned adjacent rewrite \
                 side[{side_index}]: stale successor {next:?}"
            ),
        })?;
    let curve = brep
        .edges
        .get(source_data.edge)
        .and_then(|edge| edge.curve.clone())
        .or_else(|| {
            let start = brep.vertices.get(source_data.origin)?.position;
            let end = brep.vertices.get(expected_end)?.position;
            CurveOnEdge::line_segment(start, end, brep.tolerance.linear())
        })
        .ok_or_else(|| ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cannot clone residual-owned adjacent rewrite \
                 side[{side_index}] {source:?}: source edge has no curve and endpoints cannot \
                 author a bounded line"
            ),
        })?;
    let (edge, replacement) = insert_boundary_half_edge(brep, source_data.origin);
    brep.edges[edge].curve = Some(curve);
    brep.edges[edge].role = TopologyRole::Production;
    if support_split_half_edge_end_vertex(brep, replacement) != Some(expected_end) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition cloned residual-owned adjacent rewrite \
                 side[{side_index}] {replacement:?} does not end at {expected_end:?}"
            ),
        });
    }
    Ok(replacement)
}

fn half_edge_belongs_to_third_live_loop(
    brep: &BRep,
    half_edge: HalfEdgeHandle,
    cap_loop: LoopHandle,
    target_loop: LoopHandle,
) -> bool {
    let Some(half_edge_data) = brep.half_edges.get(half_edge) else {
        return false;
    };
    half_edge_data.loop_ != LoopHandle::default()
        && half_edge_data.loop_ != cap_loop
        && half_edge_data.loop_ != target_loop
        && brep.loops.contains_key(half_edge_data.loop_)
}

fn validate_boundary_preserving_adjacent_loop_repair(
    brep: &BRep,
    rewrite: &AdjacentLoopRewriteRequirement,
    loop_handle: LoopHandle,
    replacement: HalfEdgeHandle,
    repaired_cycle: &[HalfEdgeHandle],
) -> Result<(), ApexBuildError> {
    if repaired_cycle.contains(&rewrite.stolen_cap_side) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition boundary-bearing adjacent loop repair kept stolen \
                 side {:?} in loop {loop_handle:?}",
                rewrite.stolen_cap_side
            ),
        });
    }
    let replacement_data = brep.half_edges.get(replacement).ok_or_else(|| {
        ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition boundary-bearing adjacent loop repair produced \
                 stale replacement {replacement:?}"
            ),
        }
    })?;
    if replacement_data.loop_ != loop_handle {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition boundary-bearing adjacent loop repair assigned \
                 replacement {replacement:?} to {:?}, expected {loop_handle:?}",
                replacement_data.loop_
            ),
        });
    }
    if !brep.boundary_edges.contains(&replacement_data.edge) {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition boundary-bearing adjacent loop repair replacement \
                 {replacement:?} is not a boundary edge"
            ),
        });
    }
    let actual_cycle = collect_loop_half_edges(brep, loop_handle)?;
    if actual_cycle != repaired_cycle {
        return Err(ApexBuildError::PlanMismatch {
            reason: format!(
                "planned repeated-apex partition boundary-bearing adjacent loop repair did not \
                 commit the expected loop cycle for {loop_handle:?}: got {actual_cycle:?}, \
                 expected {repaired_cycle:?}"
            ),
        });
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use popkern_core::geometry::Point3;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::{Face, Loop, Shell, Solid, Vertex};

    use super::*;

    #[test]
    fn boundary_bearing_adjacent_loop_repairs_with_same_direction_boundary_span() {
        let mut brep = BRep::new(Tolerance::default());
        let solid = brep.solids.insert(Solid {
            shell: Default::default(),
        });
        let shell = brep.shells.insert(Shell {
            face: Default::default(),
            solid,
        });
        brep.solids[solid].shell = shell;
        let vertices = [
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
        ]
        .map(|position| {
            brep.vertices.insert(Vertex {
                position,
                outgoing: Default::default(),
                name: None,
                tolerance: None,
            })
        });
        let cap_face = insert_empty_face(&mut brep, shell);
        let adjacent_face = insert_empty_face(&mut brep, shell);
        let adjacent_loop = brep.faces[adjacent_face].outer_loop;
        let (_shared_edge, stolen_cap_side, _cap_side_twin) =
            insert_unwired_edge_pair(&mut brep, vertices[0], vertices[1]);
        brep.edges[brep.half_edges[stolen_cap_side].edge].curve = Some(
            CurveOnEdge::line_segment(
                brep.vertices[vertices[0]].position,
                brep.vertices[vertices[1]].position,
                brep.tolerance.linear(),
            )
            .expect("shared edge has nondegenerate curve"),
        );
        let boundary_abutting = insert_boundary_span(&mut brep, vertices[1], vertices[2]);
        let boundary_closing = insert_boundary_span(&mut brep, vertices[2], vertices[0]);
        wire_existing_loop(
            &mut brep,
            adjacent_loop,
            adjacent_face,
            &[stolen_cap_side, boundary_abutting, boundary_closing],
        )
        .expect("synthetic adjacent loop wires");

        let rewrites = plan_adjacent_loop_rewrites_for_reauthored_cap_cycle(
            &brep,
            cap_face,
            &[stolen_cap_side],
        )
        .expect("adjacent loop rewrite plans");
        assert_eq!(rewrites.len(), 1);
        assert_eq!(
            rewrites[0].rewrite_kind,
            AdjacentLoopRewriteKind::BoundaryPreservingAdjacentLoopRepairRequired
        );
        assert!(
            rewrites[0]
                .diagnostic_label()
                .contains("BoundaryBearingAdjacentLoopRepairRequired")
        );
        let rewrite_plan = CapBoundaryRewritePlan {
            cap_face,
            cap_loop: brep.faces[cap_face].outer_loop,
            required: true,
            fillet_end_count: 1,
            face_edge_count: 0,
            current_cap_cycle: vec![stolen_cap_side],
            reauthored_cap_cycle: vec![stolen_cap_side],
            adjacent_loop_rewrites: rewrites,
        };

        rewrite_plan
            .validate_adjacent_loop_rewrite_readiness()
            .expect("boundary-preserving repair is a supported rewrite kind");
        rewrite_plan
            .apply_adjacent_loop_rewrites(&mut brep)
            .expect("boundary-preserving repair applies");

        let repaired_cycle =
            collect_loop_half_edges(&brep, adjacent_loop).expect("repaired adjacent loop");
        assert_eq!(repaired_cycle.len(), 3);
        assert!(!repaired_cycle.contains(&stolen_cap_side));
        assert_eq!(repaired_cycle[1], boundary_abutting);
        assert_eq!(repaired_cycle[2], boundary_closing);
        assert_eq!(
            brep.half_edges[repaired_cycle[0]].origin,
            brep.half_edges[stolen_cap_side].origin
        );
        assert_eq!(
            support_split_half_edge_end_vertex(&brep, repaired_cycle[0]),
            Some(vertices[1])
        );
        assert!(
            brep.boundary_edges
                .contains(&brep.half_edges[repaired_cycle[0]].edge)
        );
    }

    fn insert_empty_face(
        brep: &mut BRep,
        shell: popkern_core::topology::ShellHandle,
    ) -> FaceHandle {
        let face = brep.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: vec![],
            shell,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let loop_ = brep.loops.insert(Loop {
            half_edge: Default::default(),
            face,
        });
        brep.faces[face].outer_loop = loop_;
        face
    }

    fn insert_boundary_span(
        brep: &mut BRep,
        start: VertexHandle,
        end: VertexHandle,
    ) -> HalfEdgeHandle {
        let (edge, half_edge) = insert_boundary_half_edge(brep, start);
        brep.edges[edge].curve = Some(
            CurveOnEdge::line_segment(
                brep.vertices[start].position,
                brep.vertices[end].position,
                brep.tolerance.linear(),
            )
            .expect("boundary span has nondegenerate curve"),
        );
        half_edge
    }
}
