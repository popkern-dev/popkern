// The arm's joints, one after another: a film reads the built assembly
// and sweeps its joint parameters (placement-only: nothing rebuilds),
// starting from the pose the document is in.
let doc = document();
let film = storyboard(#{ look: "studio", camera: #{ azimuth: -60.0, elevation: 18.0, zoom: 1.0 } });
let q = [];
for k in 1..7 { q.push(doc.value("j" + k)); }

film.title("Six-axis arm", "Six cycloidal joints, every part modelled", 3.0, |s| {
    s.orbit(30.0);
});
film.step("Base and shoulder", "j1 turns the column, j2 lifts the upper arm", 3.5, |s| {
    s.joint("j1", q[0], q[0] + 50.0);
    s.joint("j2", q[1], q[1] - 25.0, #{ at: 0.6 });
});
film.step("Elbow", "j3", 2.5, |s| {
    s.joint("j3", q[2], q[2] + 40.0);
});
film.step("Wrist", "j4, j5 and j6", 3.5, |s| {
    s.joint("j4", q[3], q[3] + 90.0);
    s.joint("j5", q[4], q[4] - 45.0, #{ at: 0.6 });
    s.joint("j6", q[5], q[5] + 180.0, #{ at: 1.2 });
});
film.step("Grip", "The fingers close", 2.0, |s| {
    s.joint("grip", doc.value("grip"), 8.0);
});
film.step("Back to the start", "", 3.5, |s| {
    s.joint("j1", q[0] + 50.0, q[0]);
    s.joint("j2", q[1] - 25.0, q[1]);
    s.joint("j3", q[2] + 40.0, q[2]);
    s.joint("j4", q[3] + 90.0, q[3]);
    s.joint("j5", q[4] - 45.0, q[4]);
    s.joint("j6", q[5] + 180.0, q[5]);
    s.joint("grip", 8.0, doc.value("grip"));
    s.orbit(-30.0);
});
film
