// Assembly guide for the angled fitting: each fastener slides in along
// its own axis and turns only where its thread engages, one lead a
// turn — the M5 flat-heads (0.8) for the last 5 mm, into the bench
// under the flange; the M6 cap screw (1.0) through its 8 mm in the
// angled tap. `insert` and `screw_in` read each screw's axis (its
// insertion frame on its hole) and lead from the built assembly. A film
// reads the built assembly and never changes it; the parts it names are
// the assembly's occurrences.
let film = storyboard(#{ look: "toon", camera: #{ azimuth: -55.0, elevation: 28.0, zoom: 1.05 } });

film.title("Angled fitting", "A wedge bracket that holds a part at an angle", 2.5, |s| {
    s.orbit(30.0);
    for part in ["angled_screw", "mount_screw_a", "mount_screw_b"] {
        s.show(part, true, #{ at: 2.5 });
    }
});

film.step("1. Screw the bracket down", "Two M5 flat-heads, flush in their countersinks", 3.5, |s| {
    // 23 mm straight down the clearance hole, then the last 5 mm turning
    // as the thread bites.
    for screw in [["mount_screw_a", 0.0], ["mount_screw_b", 0.8]] {
        let at = screw[1];
        s.insert(screw[0], #{ distance: 23.0, at: at, duration: 0.8 });
        s.screw_in(screw[0], #{ distance: 5.0, at: at + 0.8, duration: 1.4, ease: "linear" });
    }
    s.camera(#{ azimuth: -20.0, elevation: 35.0, zoom: 1.5, target: "bracket_1" });
});

film.step("2. Thread the cap screw into the angled tap", "M6 x 8, on the tap's own axis", 3.5, |s| {
    // 18 mm straight in, then its 8 mm thread turning into the tap.
    s.insert("angled_screw", #{ distance: 18.0, duration: 1.0 });
    s.screw_in("angled_screw", #{ distance: 8.0, at: 1.0, duration: 2.3, ease: "linear" });
    s.camera(#{ azimuth: -105.0, elevation: 24.0, zoom: 1.35 });
});

film.step("Ready", "Four parts, every screw on its hole's axis", 3.0, |s| {
    s.camera(#{ azimuth: 40.0, elevation: 25.0, zoom: 1.0, target: "all" });
});
film.hold(0.5);
film
