//! PK-278 setback cap presolve — the parallel solve phase of the
//! two-phase corner construction.
//!
//! The expensive numeric work of a setback corner (KP20 border
//! extraction + §4.4 cap assembly) is a pure function of the corner's
//! fillet surfaces and apex position — it never reads the BRep. This
//! module precomputes that solve for every prospective setback corner
//! in parallel, before the serial vertex-blend loop runs.
//!
//! # Correctness contract
//!
//! A precomputed entry is consumed only when the consumer's inputs
//! match the entry's recorded inputs **exactly** (structural `==` on
//! the surfaces, trim curves, and apex position — bit-level for every
//! float). On any mismatch the consumer solves inline, byte-identical
//! to the serial pipeline. The equality gate is deliberately exact,
//! not tolerance-based: a tolerance would admit near-miss inputs and
//! break bit-identical output.
//!
//! Determinism: the solve is single-threaded pure math (all faer
//! calls in the KP20 path run `Parallelism::None`), so the same input
//! bits produce the same output bits on any thread. The apply phase
//! consumes entries in the unchanged serial corner order, so handle
//! allocation is unchanged.

use std::collections::BTreeMap;

use popkern_core::geometry::Point3;
use popkern_core::geometry::gregory::{Kp20PaperCapReport, Kp20PaperGuideBorder};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, VertexHandle};

use crate::construct::EdgeGeometry;
use crate::construct::blend_topo::build_fillet_surface_list;
use crate::construct::kp20_paper_border_extract::{
    Kp20PaperBorderExtractReport, SetbackJunctionPolicy, kp20_paper_border_from_fillet_setback_with,
};
use crate::fillet_surface::FilletSurfaceResult;

/// The pure solve output for one setback corner.
#[derive(Debug, Clone)]
pub(crate) struct SetbackCapSolve {
    /// The extracted guide border (consumed by the cap fan builder).
    pub border: Kp20PaperGuideBorder,
    /// Extraction evidence, including the setback arcs the surgery
    /// inserts into the topology.
    pub extract: Kp20PaperBorderExtractReport,
    /// The §4.4 cap assembly, or its error formatted exactly as the
    /// inline fan-builder call site formats it.
    pub assembly: Result<Kp20PaperCapReport, String>,
}

/// One presolved corner: the exact inputs used plus the solve result.
#[derive(Debug)]
pub(crate) struct SetbackCapPresolveEntry {
    fillet_surfaces: Vec<FilletSurfaceResult>,
    vertex_position: Point3,
    solve: Result<SetbackCapSolve, String>,
}

impl SetbackCapPresolveEntry {
    /// Return the presolved result only when `fillet_surfaces` and
    /// `vertex_position` match the recorded inputs exactly.
    pub(crate) fn matching(
        &self,
        fillet_surfaces: &[FilletSurfaceResult],
        vertex_position: Point3,
    ) -> Option<&Result<SetbackCapSolve, String>> {
        (self.vertex_position == vertex_position && self.fillet_surfaces == fillet_surfaces)
            .then_some(&self.solve)
    }
}

/// The presolved corners of one `apply_blend_geometry` run.
pub(crate) type SetbackCapPresolve = BTreeMap<VertexHandle, SetbackCapPresolveEntry>;

/// The star point of the KP20 cap assembly: 25% of the way from the
/// border centroid toward the apex. Shared between the presolve and
/// the inline fan builder so both compute identical bits.
pub(crate) fn kp20_cap_star(border: &Kp20PaperGuideBorder, vertex_position: Point3) -> Point3 {
    let mut acc = [0.0_f64; 3];
    let mut count = 0.0;
    for sector in &border.sectors {
        for piece in sector {
            for cp in &piece.rows[0] {
                acc[0] += cp.x;
                acc[1] += cp.y;
                acc[2] += cp.z;
                count += 1.0;
            }
        }
    }
    let centroid = Point3::new(acc[0] / count, acc[1] / count, acc[2] / count);
    Point3::new(
        centroid.x + (vertex_position.x - centroid.x) * 0.25,
        centroid.y + (vertex_position.y - centroid.y) * 0.25,
        centroid.z + (vertex_position.z - centroid.z) * 0.25,
    )
}

/// The setback surgery's border: reflex face corners of a saddle take
/// bridge sides (PK-527). A convex corner has none and is the
/// corners-only border bit for bit. Shared by the presolve and both
/// surgery sites, so a presolve match stays bitwise.
pub(crate) fn setback_surgery_border(
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
) -> Result<
    (
        popkern_core::geometry::gregory::Kp20PaperGuideBorder,
        Kp20PaperBorderExtractReport,
    ),
    String,
> {
    kp20_paper_border_from_fillet_setback_with(
        fillet_surfaces,
        vertex_position,
        SetbackJunctionPolicy::BridgeReflexFaces,
    )
}

/// The full pure solve for one corner: border extraction plus the
/// §4.4 cap assembly. Error strings match the inline call sites
/// exactly (`surgery_inner` propagates the extractor error raw; the
/// fan builder wraps the assembly error as `"cap assembly: {e}"`).
pub(crate) fn solve_setback_cap(
    fillet_surfaces: &[FilletSurfaceResult],
    vertex_position: Point3,
) -> Result<SetbackCapSolve, String> {
    let (border, extract) = setback_surgery_border(fillet_surfaces, vertex_position)?;
    let star = kp20_cap_star(&border, vertex_position);
    let assembly = popkern_core::geometry::gregory::assemble_kp20_paper_cap_conditioned(
        &border,
        star,
        popkern_core::geometry::gregory::Kp20StripEnforcement::HardG2,
        0.0,
    )
    .map_err(|e| format!("cap assembly: {e}"));
    Ok(SetbackCapSolve {
        border,
        extract,
        assembly,
    })
}

/// Precompute the setback cap solve for every prospective setback
/// corner, in parallel when the `parallel` feature is enabled.
///
/// Eligibility mirrors the serial dispatch (`insert_vertex_blend`'s
/// setback rung) but is only a work filter: presolving a corner the
/// rung never reaches wastes a bounded amount of parallel CPU, and
/// missing one only means that corner solves inline. Correctness is
/// carried entirely by [`SetbackCapPresolveEntry::matching`].
///
/// `cancelled` skips the presolve when a §17 cancel is already
/// pending — the presolve is one uninterruptible chunk bounded by
/// the existing one-corner cancellation latency envelope.
pub(crate) fn presolve_setback_caps(
    brep: &BRep,
    vertex_edges: &BTreeMap<VertexHandle, Vec<usize>>,
    edge_geometries: &[EdgeGeometry],
    tolerance: &Tolerance,
    cancelled: bool,
) -> SetbackCapPresolve {
    if cancelled || !crate::construct::setback_vertex_topology::setback_topology_commit_enabled() {
        return BTreeMap::new();
    }
    let tasks: Vec<(VertexHandle, Vec<FilletSurfaceResult>, Point3)> = vertex_edges
        .iter()
        .filter_map(|(vertex, edge_indices)| {
            let position = brep.vertices.get(*vertex).map(|v| v.position)?;
            let fillet_surfaces = build_fillet_surface_list(edge_indices, edge_geometries);
            let min_n = crate::construct::blend_insert::setback_min_n_for(
                edge_indices,
                edge_geometries,
                &fillet_surfaces,
                tolerance,
            );
            (edge_indices.len() >= min_n && fillet_surfaces.len() >= 3).then_some((
                *vertex,
                fillet_surfaces,
                position,
            ))
        })
        .collect();

    #[cfg(feature = "parallel")]
    let solved: Vec<(VertexHandle, SetbackCapPresolveEntry)> = {
        use rayon::prelude::*;
        tasks
            .into_par_iter()
            .map(|(vertex, fillet_surfaces, position)| {
                let solve = solve_setback_cap(&fillet_surfaces, position);
                (
                    vertex,
                    SetbackCapPresolveEntry {
                        fillet_surfaces,
                        vertex_position: position,
                        solve,
                    },
                )
            })
            .collect()
    };
    #[cfg(not(feature = "parallel"))]
    let solved: Vec<(VertexHandle, SetbackCapPresolveEntry)> = tasks
        .into_iter()
        .map(|(vertex, fillet_surfaces, position)| {
            let solve = solve_setback_cap(&fillet_surfaces, position);
            (
                vertex,
                SetbackCapPresolveEntry {
                    fillet_surfaces,
                    vertex_position: position,
                    solve,
                },
            )
        })
        .collect();

    solved.into_iter().collect()
}
