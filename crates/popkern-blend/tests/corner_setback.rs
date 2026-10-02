//! PK-623: user-set corner setbacks. A corner with a user setback is
//! served at exactly those stations, or the fillet fails typed; it never
//! falls back to an automatic setback.

use std::collections::BTreeSet;

use popkern_blend::{
    BlendError, BlendFailure, BlendResult, ContinuityLevel, CornerSetback, FilletOptions,
    FilletRadius, FilletSpec, VertexMode, fillet_edges,
};
use popkern_core::geometry::{Point3, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, EdgeHandle, VertexHandle};

/// The steep variable trihedral apex (PK-616 / PK-619 fixture).
const RADII: [f64; 3] = [0.02, 0.04, 0.06];

fn apex() -> (BRep, VertexHandle, Vec<EdgeHandle>) {
    let body = BRep::n_gonal_pyramid(3, 1.0, 1.2, Tolerance::default()).expect("pyramid");
    let apex = body
        .vertices
        .iter()
        .max_by(|a, b| a.1.position.z.total_cmp(&b.1.position.z))
        .map(|(v, _)| v)
        .expect("apex");
    let mut edges: Vec<_> = body
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect();
    edges.sort();
    edges.dedup();
    (body, apex, edges)
}

fn fillet(setback: CornerSetback) -> Result<BlendResult, BlendError> {
    let (body, _, edges) = apex();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(RADII)
        .map(|(&edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let options = FilletOptions {
        continuity: ContinuityLevel::G1,
        corner_setback: setback,
        ..Default::default()
    };
    fillet_edges(&body, &specs, &options, Tolerance::default(), None)
}

/// The out-of-range failure a request produced, if it did.
fn out_of_range(result: &Result<BlendResult, BlendError>) -> Option<(EdgeHandle, f64, f64, f64)> {
    let Err(BlendError::PartialFailure { failures, .. }) = result else {
        return None;
    };
    failures.iter().find_map(|(_, f)| match f {
        BlendFailure::SetbackOutOfRange {
            edge,
            requested,
            minimum,
            maximum,
            ..
        } => Some((*edge, *requested, *minimum, *maximum)),
        _ => None,
    })
}

/// Every cap vertex on a round's cylinder, as (axis direction, vertex
/// position): the round's end section at its station.
fn cap_feet(result: &BlendResult) -> Vec<(popkern_core::geometry::Vec3, Point3)> {
    let brep = &result.brep;
    let caps: BTreeSet<_> = result.corner_patches.iter().map(|p| p.face).collect();
    let mut feet = Vec::new();
    for (_, he) in brep.half_edges.iter() {
        if !caps.contains(&brep.loops[he.loop_].face) {
            continue;
        }
        let twin = &brep.half_edges[he.twin];
        let across = brep.loops[twin.loop_].face;
        if let Some(SurfaceGeometry::Cylinder { axis, .. }) =
            brep.faces[across].surface.as_ref().map(|s| &s.geometry)
        {
            feet.push((*axis, brep.vertices[he.origin].position));
        }
    }
    feet
}

/// The setback of a cap vertex on a round: its distance along the
/// round's axis from the apex.
fn station(apex: Point3, axis: popkern_core::geometry::Vec3, p: Point3) -> f64 {
    (p - apex).dot(&axis).abs()
}

#[test]
fn a_uniform_setback_places_every_round_end_at_the_requested_station() {
    let (body, vertex, _) = apex();
    let apex_at = body.vertices[vertex].position;
    // The corner's own refusal names each round's admissible range.
    let (_, _, minimum, _) = out_of_range(&fillet(CornerSetback::uniform(1e-3)))
        .expect("a 1e-3 setback lies before the rail crossings");
    assert!(minimum > 1e-3);
    let requested = 1.5 * minimum;
    let result = match fillet(CornerSetback::uniform(requested)) {
        Ok(result) => result,
        Err(BlendError::PartialFailure { failures, .. }) => panic!("{failures:#?}"),
        Err(e) => panic!("{e}"),
    };
    assert!(result.brep.validate().is_ok());
    assert!(result.corner_patches.iter().all(|p| p.setback_cap));
    let feet = cap_feet(&result);
    assert!(feet.len() >= 6, "{} feet", feet.len());
    for (axis, p) in feet {
        let s = station(apex_at, axis, p);
        assert!(
            (s - requested).abs() < 1e-6,
            "foot at {s:.6}, requested {requested:.6}"
        );
    }
}

#[test]
fn an_end_override_moves_only_its_round() {
    let (body, vertex, edges) = apex();
    let apex_at = body.vertices[vertex].position;
    let (_, _, minimum, _) = out_of_range(&fillet(CornerSetback::uniform(1e-3)))
        .expect("a 1e-3 setback lies before the rail crossings");
    let (uniform, moved) = (1.5 * minimum, 2.0 * minimum);
    let result = fillet(CornerSetback::uniform(uniform).with_end(vertex, edges[0], moved))
        .unwrap_or_else(|e| panic!("{e}"));
    // The overridden edge's cylinder carries the moved station; the others the
    // uniform one.
    let on_overridden_edge = |axis: popkern_core::geometry::Vec3| {
        let (a, b) = (body.vertices[vertex].position, {
            let he = body
                .half_edges
                .values()
                .find(|h| h.edge == edges[0] && h.origin == vertex)
                .expect("edge 0 leaves the apex");
            body.vertices[body.half_edges[he.twin].origin].position
        });
        let dir = (b - a).normalized().expect("edge direction");
        dir.dot(&axis).abs() > 1.0 - 1e-9
    };
    let mut seen = [false, false];
    for (axis, p) in cap_feet(&result) {
        let s = station(apex_at, axis, p);
        let expected = if on_overridden_edge(axis) {
            moved
        } else {
            uniform
        };
        seen[usize::from(on_overridden_edge(axis))] = true;
        assert!(
            (s - expected).abs() < 1e-6,
            "foot at {s:.6}, expected {expected:.6}"
        );
    }
    assert_eq!(seen, [true, true]);
}

#[test]
fn a_setback_outside_the_corner_range_fails_typed() {
    let (_, _, edges) = apex();
    let short = fillet(CornerSetback::uniform(1e-3));
    let (edge, requested, minimum, maximum) = out_of_range(&short).expect("short setback");
    assert!(edges.contains(&edge));
    assert!(
        requested < minimum && minimum < maximum,
        "{requested} {minimum} {maximum}"
    );
    let long = fillet(CornerSetback::uniform(50.0));
    let (_, requested, _, maximum) = out_of_range(&long).expect("long setback");
    assert!(requested > maximum, "{requested} {maximum}");
}

#[test]
fn invalid_setbacks_are_refused_before_any_geometry() {
    let (body, vertex, edges) = apex();
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(RADII)
        .map(|(&edge, r)| FilletSpec {
            edge,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let run = |setback: CornerSetback, vertex_mode: VertexMode| {
        let options = FilletOptions {
            corner_setback: setback,
            vertex_mode,
            ..Default::default()
        };
        fillet_edges(&body, &specs, &options, Tolerance::default(), None)
    };
    let invalid = |r: Result<BlendResult, BlendError>| matches!(r, Err(BlendError::InvalidParameter { ref parameter, .. }) if parameter == "options.corner_setback");
    assert!(invalid(run(
        CornerSetback::uniform(-0.1),
        VertexMode::Smooth
    )));
    assert!(invalid(run(
        CornerSetback::uniform(f64::NAN),
        VertexMode::Smooth
    )));
    assert!(invalid(run(
        CornerSetback::uniform(0.1),
        VertexMode::Mitred
    )));
    // An override on an edge the fillet does not round.
    let other = body
        .edges
        .keys()
        .find(|e| !edges.contains(e))
        .expect("a base edge");
    assert!(invalid(run(
        CornerSetback::default().with_end(vertex, other, 0.1),
        VertexMode::Smooth
    )));
}
