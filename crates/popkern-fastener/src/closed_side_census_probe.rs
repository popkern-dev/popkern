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

    /// PK-478 E1 fact pin, revised by PK-565 (2026-09-28). E1
    /// measured the extruded bottom cap folding INTO the material
    /// (flag true, normal_at +z, material above z = 0). That was a
    /// sweep producer defect: the self-loop start cap read its
    /// half-edge as running "forward" along the profile circle. Caps
    /// now follow the PK-303 canonical-half-edge convention, and every
    /// face of a fresh extrusion folds OUTWARD. The PK-478 design
    /// still keys nothing on flag-folded normals, because foreign
    /// bodies can be heterogeneous (see the E1 revision note in
    /// docs/architecture/pk478-sense-aware-merge.md).
    #[test]
    fn extrusion_folding_is_per_face_outward() {
        use popkern_core::geometry::SurfaceEval as _;
        let tol = Tolerance::default();
        let shank = extruded_shank(tol);
        let seed = shank.faces.iter().next().expect("faces").0;
        assert_eq!(
            shank.component_material_sense(seed),
            Some(MaterialSense::FoldedOutward),
            "majority sense"
        );
        let mut cap_folded_z = Vec::new();
        for (_, f) in &shank.faces {
            let Some(sof) = &f.surface else { continue };
            if let popkern_core::geometry::SurfaceGeometry::Plane { origin, .. } = &sof.geometry {
                let n = sof.geometry.normal_at(0.0, 0.0);
                let folded = if sof.orientation { n } else { n * -1.0 };
                cap_folded_z.push((origin.z, folded.z));
            }
        }
        assert_eq!(cap_folded_z.len(), 2, "two planar caps");
        for (z, folded_z) in cap_folded_z {
            let outward = if z.abs() < 1e-9 { -1.0 } else { 1.0 };
            assert!(
                folded_z * outward > 0.99,
                "cap at z={z} folds {folded_z:+.2} — not outward; PK-565 regressed"
            );
        }
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

    /// PK-478 stage 2: the blend's serving boundary repairs the
    /// closed-rim canonical-forward lies (stage 1 measured the tip
    /// rim ViolatedByChain: the cap-side pcurve wound against its
    /// loop and the canonical bit pointed backward). Blend output
    /// must now audit clean.
    #[test]
    fn blend_tip_rim_audits_canonical_clean() {
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
            violations, 0,
            "blend output must audit canonical-forward clean (stage 2 repair at the serving boundary)"
        );
    }
}
