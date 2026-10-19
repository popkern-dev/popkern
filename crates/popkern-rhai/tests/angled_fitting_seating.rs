//! The angled fitting's cap screw seats straight on the bracket's
//! inclined face (PK-849, PK-897): the tap's entry cone clears the
//! screw's head-side thread runout, so the head bears on the face with
//! no washer and the threads mate. Assembly interference reads the pair
//! touching; `seat` supplies the tap's insertion frame. Negative
//! control: the same screw turned half a turn about
//! its axis is half a lead out of phase, and its threads overlap the
//! tap's.

fn interference(source: &str) -> popkern_assembly::InterferenceReport {
    let run = popkern_rhai::ScriptHost::new()
        .run(source, "angled_fitting")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    run.assembly
        .expect("assembly")
        .interference(&Default::default())
        .expect("interference")
}

#[test]
fn the_cap_screw_seats_on_the_inclined_face() {
    let source = std::fs::read_to_string(format!(
        "{}/examples/angled_fitting.rhai",
        env!("CARGO_MANIFEST_DIR")
    ))
    .expect("script");
    let screw = |report: &popkern_assembly::InterferenceReport| {
        report
            .pairs
            .iter()
            .find(|p| p.b == "angled_screw" || p.a == "angled_screw")
            .map(|p| p.status.clone())
    };
    let seated = interference(&source);
    assert!(
        matches!(
            screw(&seated),
            Some(popkern_assembly::PairStatus::Touching { .. })
        ),
        "{seated:?}"
    );
    // `seat` frames the screw with x = +y (the tap's seed); the same
    // frame with x = -y is the half turn.
    let placement = "a.place(cap_screw, \"angled_screw\").seat(tap);";
    assert!(
        source.contains(placement),
        "the script's screw placement moved"
    );
    let half_turn = source.replace(
        placement,
        "a.place(cap_screw, \"angled_screw\").frame([12.0, -12.0, 9.5], \
         [-9.0 / n_len, 0.0, -24.0 / n_len], [0.0, -1.0, 0.0]);",
    );
    let off_phase = interference(&half_turn);
    assert!(
        !matches!(
            screw(&off_phase),
            Some(popkern_assembly::PairStatus::Touching { .. })
                | Some(popkern_assembly::PairStatus::Clear { .. })
                | None
        ),
        "a half-turn phase error must overlap: {off_phase:?}"
    );
}
