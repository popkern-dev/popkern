//! PK-276 tessellation interactive control: pre-signalled functional
//! test, never-fired bitwise identity, and progress-stream shape.
//! (Tessellation of the small primitive fixtures is fast — the < 1 s
//! latency claim is structural via the shared detach layer, gated
//! mid-flight on the heavier boolean/fillet/smooth-mesh ops.)

use popkern_core::cancellation::CancellationToken;
use popkern_core::progress::{CaptureReporter, OperationControl};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_mesh::{MeshError, TessellationOptions, tessellate, tessellate_with_control};

#[test]
fn pre_signalled_token_cancels_tessellation_before_any_work() {
    let brep = BRep::cylinder(0.35, 1.0, Tolerance::default()).expect("cylinder builds");
    let token = CancellationToken::new();
    token.cancel();
    let result = tessellate_with_control(
        &brep,
        &TessellationOptions::default(),
        OperationControl::with_cancel(&token),
    );
    assert!(
        matches!(result, Err(MeshError::Cancelled)),
        "pre-signalled token must surface Cancelled, got {:?}",
        result.err()
    );
    assert!(brep.validate().is_ok(), "input untouched by cancelled run");
}

/// Never-fired token: the detached tessellation is bitwise identical
/// to the plain call (positions compared at the bit level — the
/// PK-277 determinism bar).
#[test]
fn never_fired_token_matches_plain_tessellation_bitwise() {
    let brep = BRep::cylinder(0.35, 1.0, Tolerance::default()).expect("cylinder builds");
    let options = TessellationOptions::default();
    let plain = tessellate(&brep, &options).expect("plain tessellation serves");

    let token = CancellationToken::new();
    let mut reporter = CaptureReporter::default();
    let controlled = tessellate_with_control(
        &brep,
        &options,
        OperationControl::new(&token, &mut reporter),
    )
    .expect("controlled tessellation serves");

    assert_eq!(plain.triangles, controlled.triangles, "triangles differ");
    assert_eq!(plain.vertices.len(), controlled.vertices.len());
    for (i, (l, r)) in plain
        .vertices
        .iter()
        .zip(controlled.vertices.iter())
        .enumerate()
    {
        let bits =
            |p: &popkern_core::geometry::Point3| (p.x.to_bits(), p.y.to_bits(), p.z.to_bits());
        assert_eq!(
            bits(&l.position),
            bits(&r.position),
            "vertex {i} position differs"
        );
    }

    // Progress shape: monotone indices, honest ceiling, known labels,
    // edge sampling first and assembly last.
    assert!(
        reporter.events.len() >= 5,
        "expected at least sampling + round + face pass + sew + orient + assembly, got {:?}",
        reporter.events.iter().map(|e| e.stage).collect::<Vec<_>>()
    );
    for (i, event) in reporter.events.iter().enumerate() {
        assert_eq!(event.stage_index, i);
        assert!(event.stage.starts_with("tessellate."));
        assert!(event.total_stages >= reporter.events.len());
    }
    assert_eq!(reporter.events[0].stage, "tessellate.edge_samples");
    assert_eq!(
        reporter.events.last().map(|e| e.stage),
        Some("tessellate.assemble_zipper")
    );
}
