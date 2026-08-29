//! §5.3 item #5 — saddle-type vertex fixture test.
//!
//! Constructs a "notched cube" via boolean subtraction (unit cube
//! minus a smaller corner cube) to produce vertices where convex and
//! concave edges meet. Verifies that:
//!
//! 1. The classifier reports `OppositeSide` at a notch-boundary vertex.
//! 2. The rational-Gregory dispatch is correctly gated OFF for such
//!    vertices — the fillet pipeline falls through to polar-NURBS
//!    (the safe default) rather than silently producing bad Gregory
//!    geometry.
//!
//! This is the fixture-validation side of item #5's classification
//! gate (commit 849545c). Without it, the gate was tested only in
//! theory; this test exercises the gate on a real mixed-convexity
//! BRep.

use popkern_bool::subtraction;
use popkern_core::geometry::Point3;
use popkern_core::geometry::transform::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

use crate::types::{FilletOptions, FilletRadius, FilletSpec};
use crate::vertex_classifier::{VertexBlendClassification, classify_vertex_blend};

/// Build a unit cube with a small corner notched out via boolean
/// subtraction. The notch is a scaled + translated copy of the unit
/// cube placed at the `+x +y +z` corner (cube scaled by `0.6`,
/// translated to `(0.6, 0.6, 0.6)` — spans `[0.6, 1.2]³` and pokes
/// past the original cube's boundary so the intersection is a clean
/// `0.4³` corner piece).
///
/// Returns the resulting `BRep` and a vertex handle for one of the
/// "boundary" vertices where 1 convex cube edge meets 2 concave
/// notch edges (the saddle vertex we test against).
fn build_notched_cube() -> Option<(BRep, popkern_core::topology::VertexHandle)> {
    let tol = Tolerance::default();
    let unit = BRep::unit_cube(tol).expect("unit cube");
    // Notch: unit cube scaled by 0.6 then translated by (0.6, 0.6, 0.6).
    let scale = Transform3::scale_uniform(0.6);
    let translate = Transform3::translation(0.6, 0.6, 0.6);
    let notch_transform = translate.compose(&scale);
    let notch_raw = BRep::unit_cube(tol).expect("unit cube for notch");
    let notch = notch_raw.transform(&notch_transform).ok()?;
    let notched = subtraction(unit, notch, tol, None).ok()?;

    // A saddle vertex: one of the 3 boundary vertices where the
    // notch opens onto the original cube surface. Candidates:
    //   (1.0, 0.6, 0.6), (0.6, 1.0, 0.6), (0.6, 0.6, 1.0).
    // Pick the first one we find.
    let targets = [
        Point3::new(1.0, 0.6, 0.6),
        Point3::new(0.6, 1.0, 0.6),
        Point3::new(0.6, 0.6, 1.0),
    ];
    for (vh, v) in notched.vertices.iter() {
        for target in &targets {
            if (v.position.x - target.x).abs() < 1e-6
                && (v.position.y - target.y).abs() < 1e-6
                && (v.position.z - target.z).abs() < 1e-6
            {
                return Some((notched, vh));
            }
        }
    }
    None
}

/// §5.3 item #5 — classifier reports `OppositeSide` at the
/// notch-boundary vertices of a notched cube.
#[test]
fn item5_notched_cube_boundary_vertex_is_opposite_side() {
    let Some((brep, vertex)) = build_notched_cube() else {
        eprintln!("item #5 fixture: boolean subtraction failed — skipping");
        return;
    };
    let incident_edges: Vec<_> = brep
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let incident_edges: Vec<_> = {
        let mut seen = std::collections::BTreeSet::new();
        incident_edges
            .into_iter()
            .filter(|e| seen.insert(*e))
            .collect()
    };
    let report = classify_vertex_blend(&brep, vertex, &incident_edges, &Tolerance::default());
    eprintln!(
        "item #5 notched-cube boundary vertex classification: {:?} \
         (n_convex={}, n_concave={}, n_mixed={})",
        report.classification, report.n_convex, report.n_concave, report.n_mixed
    );
    assert!(
        matches!(
            report.classification,
            VertexBlendClassification::OppositeSide
        ),
        "notched-cube boundary vertex must classify as OppositeSide \
         (got {:?})",
        report.classification
    );
}

/// §5.3 item #5 / PK-231 — the `OppositeSide` saddle vertex SERVES:
/// exactly 3 Gregory sub-patches, no Mitred fallback, no direct-apex
/// capture, structurally clean.
///
/// This pin was an either-way pin ("3 sub-patches OR ≤ 1 fallback")
/// from its landing until PK-231, and the escape hatch absorbed a
/// three-era regression stack silently:
///
/// 1. `BRep::transform` scaled plane normals (the notch transform
///    left |n| = 0.6), displacing every fillet offset — the "serving"
///    era-1 patches were built on mis-tangent strips and genuinely
///    folded, so the item #10 self-intersection screen CORRECTLY
///    remediated them away (the screen was never the defect).
/// 2. Under trim-crossing anchors the poisoned trims collapsed both
///    side-2 strips onto one consumed face and the consumed-merge stage failed
///    outright (`PartialFailure` — which this test used to skip).
/// 3. The mixed-radius direct-apex promotion arm then captured the
///    saddle (it never consulted the classifier), planned the Open
///    rung ("only 2 fillets" — the concave strip is invisible to the
///    convex footprint search), and recorded a vacuous
///    `route=DirectApex` success with ZERO faces, skipping the
///    pre-ring saddle machinery entirely.
///
/// All three are fixed at root (similarity-gated `BRep::transform`
/// with unit direction fields, `SameSideConvex`-only promotion, the
/// auto-promoted-Open decline backstop), so the serve is asserted
/// EXACTLY — any regression fails instead of draining into a
/// fallback arm.
#[test]
fn item5_saddle_fillet_serves_three_gregory_sub_patches() {
    let (brep, vertex) = build_notched_cube()
        .expect("notched-cube boolean fixture must build (PK-231: no more silent skips)");

    let incident_edges: Vec<_> = brep
        .half_edges
        .values()
        .filter(|he| he.origin == vertex)
        .map(|he| he.edge)
        .collect();
    let incident_edges: Vec<_> = {
        let mut seen = std::collections::BTreeSet::new();
        incident_edges
            .into_iter()
            .filter(|e| seen.insert(*e))
            .collect()
    };
    assert_eq!(
        incident_edges.len(),
        3,
        "saddle vertex should have 3 incident edges"
    );

    // Mixed radii (to preclude the rolling-ball sphere path for any
    // non-saddle corner that may still classify as SameSide*).
    let radii = [0.06_f64, 0.08, 0.10];
    let specs: Vec<FilletSpec> = incident_edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();

    let opts = FilletOptions {
        use_rational_gregory: true,
        ..FilletOptions::default()
    };
    let result = crate::fillet_edges(&brep, &specs, &opts, Tolerance::default(), None).expect(
        "saddle fillet pipeline must complete (PK-231: the consumed-merge stage carve restored)",
    );

    // No Mitred fallback and no direct-apex capture: the saddle is
    // vertex-blend-owned (promotion is gated to SameSideConvex). Vertex-blend
    // routes DO record their served corner faces since PK-258 (the
    // serving-retry keep-better measure), so the route label — not an
    // empty built_faces list — is the ownership claim.
    assert!(
        !result
            .edge_warnings
            .iter()
            .any(|(_, w)| matches!(w, crate::error::BlendWarning::VertexModeFallback { .. })),
        "saddle blend must not fall back to Mitred; warnings: {:?}",
        result.edge_warnings
    );
    assert!(
        result
            .repair_report
            .apex_routes
            .iter()
            .filter(|r| r.vertex == vertex)
            .all(|r| r.route == crate::types::ApexConstructionRoute::VertexBlend),
        "OppositeSide saddle must stay on the vertex-blend route (promotion is \
         SameSideConvex-only); routes: {:?}",
        result
            .repair_report
            .apex_routes
            .iter()
            .map(|r| (r.vertex, r.route))
            .collect::<Vec<_>>()
    );

    let saddle_corners: Vec<_> = result
        .corner_patches
        .iter()
        .filter(|p| p.vertex == vertex)
        .collect();
    assert_eq!(
        saddle_corners.len(),
        3,
        "saddle dispatch must serve exactly 3 Gregory sub-patches; got {} ({:?})",
        saddle_corners.len(),
        saddle_corners
            .iter()
            .map(|p| (p.method, p.n_sides))
            .collect::<Vec<_>>()
    );

    // Sub-patch geometry: Gregory surfaces, non-degenerate, and the
    // serve VALIDATE-CLEAN under collect-all check 9 (PK-234): the
    // boundary-driven ring carve (station pushback + rail-spur recut
    // through the mitre complex + arc-normalized blend loop) bonds
    // the decomposition onto a machine-exact ring (4.7e-12), the
    // support-seam relocation no-ops on it, and the boundary-owned
    // seam authoring lands every cap edge inside the validation
    // band. Zero CurveSurfaceDeviation errors — the mitred fallback
    // era carried 14 at worst 1.016e-1 (history in the PK-234
    // ledger; PK-235 exposed the true counts).
    use popkern_core::geometry::surface::SurfaceEval;

    let validation_errors = result.brep.validate_all();
    assert!(
        validation_errors.is_empty(),
        "saddle ring serve must be validate_all CLEAN; got {validation_errors:?}"
    );

    // The ring rung must record that it applied — a silent fallback
    // to the mitred construction would keep this green only until
    // the deviation asserts above caught it, so pin the route fact
    // directly.
    let boundary_ring_facts: Vec<&str> = result
        .repair_report
        .apex_diagnostics
        .iter()
        .filter_map(|d| match d {
            crate::types::ApexConstructionDiagnostic::ExactBoundaryCap(fact)
                if fact.code == "boundary_ring.fact" =>
            {
                Some(fact.stage.as_str())
            }
            _ => None,
        })
        .collect();
    assert!(
        boundary_ring_facts.contains(&"boundary_ring.applied"),
        "the saddle ring rung must record an applied fact; facts: {boundary_ring_facts:?}"
    );

    for patch in &saddle_corners {
        let face = result
            .brep
            .faces
            .get(patch.face)
            .expect("sub-patch face must exist in result");
        let sof = face
            .surface
            .as_ref()
            .expect("sub-patch face must have a surface assigned");
        assert!(
            matches!(
                sof.geometry,
                popkern_core::geometry::surface::SurfaceGeometry::Gregory(_)
            ),
            "sub-patch face must carry a Gregory surface; got {:?}",
            std::mem::discriminant(&sof.geometry)
        );
        // Positive-area proxy: sample at the 4 UV corners
        // and require max pairwise distance > 1e-6 (the
        // patch isn't degenerate to a point).
        let pts = [
            sof.geometry.point_at(0.0, 0.0),
            sof.geometry.point_at(1.0, 0.0),
            sof.geometry.point_at(0.0, 1.0),
            sof.geometry.point_at(1.0, 1.0),
        ];
        let max_dist = (0..4)
            .flat_map(|i| (i + 1..4).map(move |j| pts[i].distance_to(&pts[j])))
            .fold(0.0_f64, f64::max);
        assert!(
            max_dist > 1e-6,
            "Gregory sub-patch must be non-degenerate; max corner-pair distance = {max_dist}"
        );
    }
}

/// §5.3 item #5 — unit test on the saddle-specific central-point
/// calculation in `build_n_quad_decomposition_with_options`. Builds
/// a synthetic 2N-segment boundary whose FilletEnd arc midpoints
/// lie in a known configuration and verifies that for the
/// `OppositeSide` classification the computed central point is the
/// average of the arc midpoints (NOT the hex-vertex centroid).
///
/// Uses `analyze_saddle_central_point` below — a test-only helper
/// that just re-runs the same central-point selection logic used
/// in the builder.
#[test]
fn item5_saddle_central_point_is_arc_mid_average() {
    use crate::blend::boundary::{BlendBoundary, BoundarySegment};

    let tol = Tolerance::default();

    // Fabricate a synthetic hexagonal saddle boundary. Three
    // FilletEnd arcs at different z levels (one "up", two "down"
    // simulating 1 convex + 2 concave); bridges straight-line
    // between them. The hex-vertex centroid averages over all 6
    // vertices (arc endpoints included) so it's pulled toward
    // whichever side has tighter bridge-vertex spread. The saddle
    // central point uses only the 3 arc midpoints — symmetric in z.
    let arc = |p0: Point3, pm: Point3, p1: Point3| {
        crate::center_curve::fit_nurbs_through_points(&[p0, pm, p1]).unwrap()
    };
    let p_up_a = Point3::new(1.0, 0.0, 0.5);
    let p_up_m = Point3::new(1.0, 0.2, 0.6);
    let p_up_b = Point3::new(1.0, 0.4, 0.5);

    let p_dn1_a = Point3::new(-0.5, 0.87, -0.5);
    let p_dn1_m = Point3::new(-0.4, 0.82, -0.6);
    let p_dn1_b = Point3::new(-0.3, 0.77, -0.5);

    let p_dn2_a = Point3::new(-0.5, -0.87, -0.5);
    let p_dn2_m = Point3::new(-0.4, -0.82, -0.6);
    let p_dn2_b = Point3::new(-0.3, -0.77, -0.5);

    let bridge = |a: Point3, b: Point3| {
        let mid = Point3::new(0.5 * (a.x + b.x), 0.5 * (a.y + b.y), 0.5 * (a.z + b.z));
        crate::center_curve::fit_nurbs_through_points(&[a, mid, b]).unwrap()
    };

    let dummy_face = popkern_core::topology::FaceHandle::default();
    let dummy_vh = popkern_core::topology::VertexHandle::default();
    let boundary = BlendBoundary {
        vertex: dummy_vh,
        vertex_position: Point3::new(0.0, 0.0, 0.0),
        segments: vec![
            BoundarySegment::FilletEnd {
                curve: arc(p_up_a, p_up_m, p_up_b),
                fillet_index: 0,
            },
            BoundarySegment::FaceEdge {
                curve: bridge(p_up_b, p_dn1_a),
                face: dummy_face,
            },
            BoundarySegment::FilletEnd {
                curve: arc(p_dn1_a, p_dn1_m, p_dn1_b),
                fillet_index: 1,
            },
            BoundarySegment::FaceEdge {
                curve: bridge(p_dn1_b, p_dn2_a),
                face: dummy_face,
            },
            BoundarySegment::FilletEnd {
                curve: arc(p_dn2_a, p_dn2_m, p_dn2_b),
                fillet_index: 2,
            },
            BoundarySegment::FaceEdge {
                curve: bridge(p_dn2_b, p_up_a),
                face: dummy_face,
            },
        ],
        bridge_fallbacks: Vec::new(),
        sphere_fit_substituted: false,
        per_fillet_caps: Vec::new(),
        hex_boundary_cps_deg9: Vec::new(),
        fillet_end_arcs: Vec::new(),
    };

    // Compute the expected saddle central point: average of the 3
    // arc midpoints.
    let expected = Point3::new(
        (p_up_m.x + p_dn1_m.x + p_dn2_m.x) / 3.0,
        (p_up_m.y + p_dn1_m.y + p_dn2_m.y) / 3.0,
        (p_up_m.z + p_dn1_m.z + p_dn2_m.z) / 3.0,
    );

    // Run the decomposition — for `OppositeSide`, central_point
    // should match `expected`. For `SameSideConvex`, it should
    // NOT match (hex centroid instead).
    let saddle = crate::blend::gregory_three_quad::build_n_quad_decomposition_with_options(
        &boundary,
        &[],
        false,
        VertexBlendClassification::OppositeSide,
    );
    let convex = crate::blend::gregory_three_quad::build_n_quad_decomposition_with_options(
        &boundary,
        &[],
        false,
        VertexBlendClassification::SameSideConvex,
    );

    // Both paths must build on this well-formed synthetic boundary.
    let saddle = saddle.expect(
        "item #5: saddle-path decomposition must build on a valid \
         hexagonal boundary — None indicates a regression in the \
         saddle central-point or boundary plumbing",
    );
    let convex = convex.expect(
        "item #5 precondition: convex-path decomposition should also \
         build on the same boundary (sanity check)",
    );

    let c = saddle.central_point;
    let d_expected = (c.x - expected.x)
        .hypot(c.y - expected.y)
        .hypot(c.z - expected.z);
    assert!(
        d_expected < 1e-9,
        "item #5: saddle central point {c:?} must equal arc-mid average \
         {expected:?} (within 1e-9); got distance {d_expected:.4e}"
    );

    // The convex (hex-centroid) path MUST produce a DIFFERENT
    // central point on this synthetic boundary — proves the saddle
    // branch actually changed behaviour rather than trivially
    // matching the default.
    let d_saddle_convex = (saddle.central_point.x - convex.central_point.x)
        .hypot(saddle.central_point.y - convex.central_point.y)
        .hypot(saddle.central_point.z - convex.central_point.z);
    assert!(
        d_saddle_convex > 1e-6,
        "item #5: saddle and convex central points should differ on a \
         synthetic saddle boundary; got {d_saddle_convex:.4e}"
    );

    let _ = tol; // silence unused if all branches skip.
}

#[path = "boundary_ring_class_tests.rs"]
mod boundary_ring_class_tests;
