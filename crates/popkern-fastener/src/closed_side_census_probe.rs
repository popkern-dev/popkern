//! PK-478 E1 measured-fact pins + the closed-edge witness census
//! dump (ignored; re-run by hand when revising the sense-aware-merge
//! design). The census itself is the promoted structured diagnostic
//! `popkern_bool::closed_side_census`.

#[cfg(test)]
mod tests {
    use popkern_bool::closed_side_census;
    use popkern_core::tolerance::Tolerance;
    use popkern_core::topology::BRep;
    use popkern_core::validation::MaterialSense;

    fn extruded_shank(tol: Tolerance) -> BRep {
        let spec = crate::spec::ThreadSpec::parse("M6").expect("M6");
        crate::engine::extrude_shape(
            popkern_sketch::ProfileShape::Circle {
                radius: spec.major_diameter() / 2.0,
            },
            0.0,
            4.0,
            "probe",
            tol,
        )
        .expect("shank")
    }

    /// PK-478 E1 measured fact (2026-09-11): per-face material
    /// folding is HETEROGENEOUS even on a fresh, winding-honest,
    /// majority-outward extrusion — the bottom cap's flag-folded
    /// normal points INTO the material (flag true, normal_at +z,
    /// material above z = 0) while wall and top cap fold outward.
    /// `component_material_sense` is a majority statistic, not a
    /// per-face invariant; any witness or walk convention keyed on
    /// the flag-folded normal being outward is falsified by this
    /// body. If this pin ever breaks, the PK-478 design assumptions
    /// (docs/architecture/pk478-sense-aware-merge.md) must be
    /// revisited before touching the closed-edge witness.
    #[test]
    fn extrusion_folding_is_per_face_heterogeneous() {
        use popkern_core::geometry::SurfaceEval as _;
        let tol = Tolerance::default();
        let shank = extruded_shank(tol);
        let seed = shank.faces.iter().next().expect("faces").0;
        assert_eq!(
            shank.component_material_sense(seed),
            Some(MaterialSense::FoldedOutward),
            "majority sense"
        );
        let mut bottom_folded_z = None;
        for (_, f) in &shank.faces {
            let Some(sof) = &f.surface else { continue };
            if let popkern_core::geometry::SurfaceGeometry::Plane { origin, .. } = &sof.geometry
                && origin.z.abs() < 1e-9
            {
                let n = sof.geometry.normal_at(0.0, 0.0);
                let folded = if sof.orientation { n } else { n * -1.0 };
                bottom_folded_z = Some(folded.z);
            }
        }
        let z = bottom_folded_z.expect("bottom cap found");
        assert!(
            z > 0.0,
            "bottom cap folded normal now points {z:+.2} — the measured \
             heterogeneity fact changed; revisit the PK-478 design"
        );
    }

    /// Exploration dump for design-doc revisions: witness candidates
    /// per closed edge across the E1 populations. Run by hand:
    /// `cargo nextest run -p popkern-fastener e1_census_dump \
    ///  --run-ignored all --no-capture`.
    #[test]
    #[ignore = "PK-478 E1 exploration dump — run by hand when revising the sense-aware-merge design"]
    fn e1_census_dump() {
        let tol = Tolerance::default();
        let spec = crate::spec::ThreadSpec::parse("M6").expect("M6");
        let shank = extruded_shank(tol);
        let mut reversed = shank.clone();
        reversed.reverse_body_orientation();
        let chamfered = crate::builders::tip_chamfer(&shank, &spec, 4.0, tol).expect("chamfer");
        let mut rechamfered = chamfered.clone();
        rechamfered.reverse_body_orientation();
        for (label, body) in [
            ("outward-cylinder", &shank),
            ("reversed-cylinder", &reversed),
            ("inward-chamfered", &chamfered),
            ("reversed-chamfered", &rechamfered),
        ] {
            let seed = body.faces.iter().next().expect("faces").0;
            eprintln!(
                "[e1] === {label}: sense {:?} ===",
                body.component_material_sense(seed)
            );
            for row in closed_side_census(body) {
                eprintln!("[e1] {label} {row:?}");
            }
        }
    }
    /// PK-478 stage-1 lie census: audit the canonical-forward
    /// convention across the E1 populations + a drilled plate (the
    /// PK-465 bore-wall class). Prints verdicts; promoted assertions
    /// follow once the population is understood. Run by hand.
    #[test]
    #[ignore = "PK-478 stage-1 lie census — run by hand"]
    fn canonical_forward_lie_census() {
        let tol = Tolerance::default();
        let spec = crate::spec::ThreadSpec::parse("M6").expect("M6");
        let shank = extruded_shank(tol);
        let mut reversed = shank.clone();
        reversed.reverse_body_orientation();
        let chamfered = crate::builders::tip_chamfer(&shank, &spec, 4.0, tol).expect("chamfer");
        let plate = crate::engine::extrude_shape(
            popkern_sketch::ProfileShape::Polygon {
                points: vec![[-10.0, -10.0], [10.0, -10.0], [10.0, 10.0], [-10.0, 10.0]],
            },
            0.0,
            8.0,
            "probe",
            tol,
        )
        .expect("plate");
        let drill = crate::engine::extrude_shape(
            popkern_sketch::ProfileShape::Circle { radius: 2.0 },
            -1.0,
            10.0,
            "probe",
            tol,
        )
        .expect("drill");
        let bored = crate::engine::subtract(&plate, &drill, "probe", tol).expect("bored plate");
        for (label, body) in [
            ("outward-cylinder", &shank),
            ("reversed-cylinder", &reversed),
            ("inward-chamfered", &chamfered),
            ("bored-plate", &bored),
        ] {
            eprintln!("[cf] === {label} ===");
            for v in body.audit_canonical_forward() {
                eprintln!("[cf] {label} {:?} {:?}", v.edge, v.state);
            }
        }
    }

    /// PK-478 stage-1 lie specimen (measured 2026-09-11): the blend
    /// tip chamfer emits its tip rim with the canonical-forward
    /// convention VIOLATED (the lying half-edge convicted by UV
    /// chain continuity). This pin asserts the audit DETECTS the
    /// producer lie; when PK-478 stage 2 fixes the blend producer,
    /// flip this pin to assert no violations instead.
    #[test]
    fn audit_detects_blend_tip_rim_canonical_lie() {
        use popkern_core::validation::CanonicalForwardState;
        let tol = Tolerance::default();
        let spec = crate::spec::ThreadSpec::parse("M6").expect("M6");
        let shank = extruded_shank(tol);
        let chamfered = crate::builders::tip_chamfer(&shank, &spec, 4.0, tol).expect("chamfer");
        let violations = chamfered
            .audit_canonical_forward()
            .into_iter()
            .filter(|v| {
                matches!(
                    v.state,
                    CanonicalForwardState::Violated | CanonicalForwardState::ViolatedByChain { .. }
                )
            })
            .count();
        assert_eq!(
            violations, 1,
            "expected exactly the tip-rim producer lie (stage 2 flips this pin to 0)"
        );
    }
}
