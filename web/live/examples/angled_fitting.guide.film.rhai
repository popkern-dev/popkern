// Assembly guide for the angled fitting: each fastener slides in along
// its own axis and turns only where its thread engages, one lead a
// turn — the M5 flat-heads (0.8) for the last 5 mm, into the bench
// under the flange; the M6 cap screw (1.0) through its 8 mm in the
// angled tap. A film reads the built assembly and never changes it; the parts
// it names are the assembly's occurrences.
let doc = document();
let film = storyboard(#{ look: "toon", camera: #{ azimuth: -55.0, elevation: 28.0, zoom: 1.05 } });

// The tap's axis (into the bracket) and the way out of it.
let n = sqrt(9.0 * 9.0 + 24.0 * 24.0);
let tap = [-9.0 / n, 0.0, -24.0 / n];
let out = [9.0 / n * 26.0, 0.0, 24.0 / n * 26.0];
// The cap screw's 26 mm out: 18 straight, then its 8 mm thread.
let reach = [9.0 / n * 18.0, 0.0, 24.0 / n * 18.0];
let engaged = [9.0 / n * 8.0, 0.0, 24.0 / n * 8.0];

film.title("Angled fitting", "A wedge bracket that holds a part at an angle", 2.5, |s| {
    s.orbit(30.0);
    for part in ["washer_1", "angled_screw", "mount_screw_a", "mount_screw_b"] {
        s.show(part, true, #{ at: 2.5 });
    }
});

film.step("1. Screw the bracket down", "Two M5 flat-heads, flush in their countersinks", 3.5, |s| {
    // 23 mm straight down the clearance hole, then the last 5 mm turning
    // 5 / 0.8 = 6.25 times (2250°) as the thread bites.
    for screw in [["mount_screw_a", -5.5, 0.0], ["mount_screw_b", -18.5, 0.8]] {
        let at = screw[2];
        s.move(screw[0], #{ from: [0.0, 0.0, 23.0], at: at, duration: 0.8 });
        s.move(screw[0], #{ from: [0.0, 0.0, 5.0], at: at + 0.8, duration: 1.4, ease: "linear" });
        s.turn(screw[0], axis_at([34.0, screw[1], 5.0], [0.0, 0.0, 1.0]), 2250.0, 0.0, #{ at: at + 0.8, duration: 1.4, ease: "linear" });
    }
    s.camera(#{ azimuth: -20.0, elevation: 35.0, zoom: 1.5, target: "bracket_1" });
});

film.step("2. A washer on the inclined face", "ISO 7089 M6, flat on the 20.6° face", 2.5, |s| {
    s.move("washer_1", #{ from: out, duration: 1.6 });
    s.highlight("washer_1", #{ at: 1.2, duration: 1.3 });
    s.camera(#{ azimuth: -80.0, elevation: 32.0, zoom: 1.5, target: "bracket_1" });
});

film.step("3. Thread the cap screw into the angled tap", "M6 x 8, on the tap's own axis", 3.5, |s| {
    s.move("angled_screw", #{ from: reach, duration: 1.0 });
    s.move("angled_screw", #{ from: engaged, at: 1.0, duration: 2.3, ease: "linear" });
    s.turn("angled_screw", axis_at([12.0, -12.0, 9.5], tap), -2880.0, 0.0, #{ at: 1.0, duration: 2.3, ease: "linear" });
    s.camera(#{ azimuth: -105.0, elevation: 24.0, zoom: 1.35 });
});

film.step("Ready", "Five parts, every screw on its hole's axis", 3.0, |s| {
    s.camera(#{ azimuth: 40.0, elevation: 25.0, zoom: 1.0, target: "all" });
});
film.hold(0.5);
film
