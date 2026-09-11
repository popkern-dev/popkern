//! PK-478 stage-2 sizing: canonical-forward audit across a
//! representative producer battery. The summary counts are printed;
//! run by hand (ignored) while the producer-fix campaign is open:
//! `cargo nextest run -p popkern-examples canonical_forward_census \
//!  --run-ignored all --no-capture`

use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_core::validation::CanonicalForwardState;
use popkern_part::prelude::*;

fn audit(label: &str, brep: &BRep) -> (usize, usize, usize) {
    let (mut ok, mut violated, mut indet) = (0, 0, 0);
    for v in brep.audit_canonical_forward() {
        match v.state {
            CanonicalForwardState::Ok | CanonicalForwardState::OkByChain { .. } => ok += 1,
            CanonicalForwardState::Violated | CanonicalForwardState::ViolatedByChain { .. } => {
                violated += 1;
                eprintln!("[cf-census] {label}: VIOLATED {:?} ({:?})", v.edge, v.state);
            }
            CanonicalForwardState::Indeterminate(reason) => {
                indet += 1;
                eprintln!("[cf-census] {label}: indeterminate {:?} ({reason})", v.edge);
            }
        }
    }
    eprintln!("[cf-census] {label}: ok {ok} violated {violated} indeterminate {indet}");
    (ok, violated, indet)
}

#[test]
#[ignore = "PK-478 stage-2 sizing census — run by hand"]
fn canonical_forward_producer_census() {
    let tol = Tolerance::default();
    let mut total_violated = 0usize;

    // Constructors.
    for (label, body) in [
        ("cylinder", BRep::cylinder(2.0, 5.0, tol).expect("cylinder")),
        ("sphere", BRep::sphere(2.0, tol).expect("sphere")),
        ("torus", BRep::torus(3.0, 1.0, tol).expect("torus")),
    ] {
        total_violated += audit(label, &body).1;
    }

    // Part-level producers.
    let disk = Profile::on(Plane::XY)
        .circle(2.0)
        .extrude(5.0)
        .expect("extrude");
    total_violated += audit("extruded-cylinder", disk.brep()).1;

    let plate = Profile::on(Plane::XY)
        .rect(10.0, 10.0)
        .extrude(4.0)
        .expect("plate");
    let bored = plate
        .subtract(
            &Profile::on(Plane::XY)
                .circle(1.5)
                .extrude(4.0)
                .expect("drill"),
        )
        .expect("bored");
    total_violated += audit("bored-plate", bored.brep()).1;

    let tilted = Profile::on(Plane::XY)
        .rect(10.0, 10.0)
        .extrude(4.0)
        .expect("plate")
        .subtract(
            &Profile::on(Plane::XY)
                .circle(1.5)
                .extrude(9.0)
                .expect("drill")
                .rotate_x(0.4_f64)
                .expect("tilt")
                .translate(0.0, 0.0, -2.0)
                .expect("place"),
        )
        .expect("tilted bore");
    total_violated += audit("tilted-bored-plate", tilted.brep()).1;

    let filleted = Profile::on(Plane::XY)
        .circle(2.0)
        .extrude(4.0)
        .expect("puck")
        .fillet(0.5, |e| e.circular().at_height(4.0))
        .expect("filleted puck rim");
    total_violated += audit("filleted-puck", filleted.brep()).1;

    let chamfered = Profile::on(Plane::XY)
        .circle(2.0)
        .extrude(4.0)
        .expect("puck")
        .chamfer(0.4, |e| e.circular().at_height(4.0))
        .expect("chamfered puck rim");
    total_violated += audit("chamfered-puck", chamfered.brep()).1;

    // Revolve (sweep-level): a rectangle meridian into a full ring.
    let ring_wire = popkern_core::geometry::NurbsCurve {
        degree: 1,
        knots: vec![0.0, 0.0, 0.25, 0.5, 0.75, 1.0, 1.0],
        control_points: vec![
            popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
            popkern_core::geometry::Point3::new(2.0, 0.0, 0.0),
            popkern_core::geometry::Point3::new(2.0, 0.0, 3.0),
            popkern_core::geometry::Point3::new(1.0, 0.0, 3.0),
            popkern_core::geometry::Point3::new(1.0, 0.0, 0.0),
        ],
        weights: vec![1.0; 5],
    };
    let revolved = popkern_sweep::revolve(
        &popkern_sweep::SweepProfile::Wire(ring_wire),
        &popkern_core::geometry::Point3::new(0.0, 0.0, 0.0),
        &popkern_core::geometry::Vec3::new(0.0, 0.0, 1.0),
        std::f64::consts::TAU,
        tol,
        None,
    )
    .expect("revolved ring");
    total_violated += audit("revolved-ring", &revolved.brep).1;

    eprintln!("[cf-census] TOTAL violated: {total_violated}");
}

/// Transient localization probe for the "pcurve traversal
/// unmeasurable" indeterminate class (PK-478 stage 3 mechanism b):
/// which ingredient is missing per side on the boolean-output and
/// blend-output closed edges.
#[test]
#[ignore = "PK-478 stage-3 localization — run by hand"]
fn unmeasurable_side_localization() {
    let tilted = Profile::on(Plane::XY)
        .rect(10.0, 10.0)
        .extrude(4.0)
        .expect("plate")
        .subtract(
            &Profile::on(Plane::XY)
                .circle(1.5)
                .extrude(9.0)
                .expect("drill")
                .rotate_x(0.4_f64)
                .expect("tilt")
                .translate(0.0, 0.0, -2.0)
                .expect("place"),
        )
        .expect("tilted bore");
    let chamfered = Profile::on(Plane::XY)
        .circle(2.0)
        .extrude(4.0)
        .expect("puck")
        .chamfer(0.4, |e| e.circular().at_height(4.0))
        .expect("chamfered puck rim");
    for (label, brep) in [
        ("tilted-bored", tilted.brep()),
        ("chamfered-puck", chamfered.brep()),
    ] {
        for v in brep.audit_canonical_forward() {
            let CanonicalForwardState::Indeterminate(reason) = v.state else {
                continue;
            };
            let edge = &brep.edges[v.edge];
            let he = &brep.half_edges[edge.half_edge];
            let tw = &brep.half_edges[he.twin];
            let side = |h: &popkern_core::topology::HalfEdge| {
                let face = brep.loops[h.loop_].face;
                let surf = brep.faces[face]
                    .surface
                    .as_ref()
                    .map(|s| format!("{:?}", std::mem::discriminant(&s.geometry)));
                let pc = h.pcurve.as_ref().map(|p| {
                    let (t0, t1) = p.parameter_range();
                    format!("range ({t0:.3},{t1:.3})")
                });
                format!("face {face:?} surf {surf:?} pcurve {pc:?}")
            };
            eprintln!(
                "[loc3] {label} {:?} ({reason}) curve {} | canonical: {} | twin: {}",
                v.edge,
                edge.curve.is_some(),
                side(he),
                side(tw)
            );
            // Step-replicate the traversal claim on each side.
            for (tag, h) in [("canonical", he), ("twin", tw)] {
                use popkern_core::geometry::SurfaceEval as _;
                let Some(pc) = h.pcurve.as_ref() else {
                    continue;
                };
                let face = brep.loops[h.loop_].face;
                let Some(sof) = brep.faces[face].surface.as_ref() else {
                    continue;
                };
                let Some(curve) = edge.curve.as_ref() else {
                    continue;
                };
                let (t0, t1) = pc.parameter_range();
                let tm = t0 + (t1 - t0) * 0.5;
                let te = t0 + (t1 - t0) * 0.505;
                let a = pc.point_at(tm);
                let b = pc.point_at(te);
                let p1 = sof.geometry.point_at(a.x, a.y);
                let p2 = sof.geometry.point_at(b.x, b.y);
                let step = (p2 - p1).length();
                let inv1 = popkern_core::geometry::curve_inversion::invert_curve_at_point(
                    curve,
                    p1,
                    &brep.tolerance,
                );
                let inv2 = popkern_core::geometry::curve_inversion::invert_curve_at_point(
                    curve,
                    p2,
                    &brep.tolerance,
                );
                eprintln!(
                    "[loc3]   {tag}: uv ({:.4},{:.4})->({:.4},{:.4}) step3d {step:.3e} inv {:?} {:?} curve-range ({:.3},{:.3})",
                    a.x, a.y, b.x, b.y, inv1, inv2, curve.t_start, curve.t_end
                );
            }
        }
    }
}

/// Transient localization probe for the blend strip rim lie.
#[test]
#[ignore = "PK-478 stage-2 localization — run by hand"]
fn blend_rim_lie_localization() {
    let chamfered = Profile::on(Plane::XY)
        .circle(2.0)
        .extrude(4.0)
        .expect("puck")
        .chamfer(0.4, |e| e.circular().at_height(4.0))
        .expect("chamfered puck rim");
    let brep = chamfered.brep();
    for (fh, f) in &brep.faces {
        let kind = f.surface.as_ref().map(|s| {
            format!(
                "{:?} flag {}",
                std::mem::discriminant(&s.geometry),
                s.orientation
            )
        });
        eprintln!("[loc] face {fh:?} {kind:?}");
    }
    for v in brep.audit_canonical_forward() {
        eprintln!("[loc] verdict {:?} {:?}", v.edge, v.state);
        if let popkern_core::validation::CanonicalForwardState::ViolatedByChain {
            lying_half_edge,
        } = v.state
        {
            let he = &brep.half_edges[lying_half_edge];
            let face = brep.loops[he.loop_].face;
            let edge = &brep.edges[he.edge];
            let canonical_is_liar = edge.half_edge == lying_half_edge;
            eprintln!(
                "[loc] lying he {lying_half_edge:?} on face {face:?}; canonical-is-liar {canonical_is_liar}; has pcurve {}",
                he.pcurve.is_some()
            );
            if let Some(pc) = &he.pcurve {
                let (t0, t1) = pc.parameter_range();
                let a = pc.point_at(t0);
                let b = pc.point_at(t1);
                eprintln!(
                    "[loc] lying pcurve ({:.3},{:.3})->({:.3},{:.3})",
                    a.x, a.y, b.x, b.y
                );
            }
            let tw = &brep.half_edges[he.twin];
            let twf = brep.loops[tw.loop_].face;
            if let Some(pc) = &tw.pcurve {
                let (t0, t1) = pc.parameter_range();
                let a = pc.point_at(t0);
                let b = pc.point_at(t1);
                eprintln!(
                    "[loc] twin he {:?} on face {twf:?} pcurve ({:.3},{:.3})->({:.3},{:.3})",
                    he.twin, a.x, a.y, b.x, b.y
                );
            }
        }
    }
}
