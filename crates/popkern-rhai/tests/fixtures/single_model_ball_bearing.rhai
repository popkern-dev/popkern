// A deep-groove ball bearing: two races revolved from sections whose
// grooves are exact arcs, and a ring of balls placed by trigonometry
// in a loop — ten bodies from one script, every surface an exact
// plane, cylinder, sphere or torus.
let m = model();

// Groove geometry: balls of radius 0.15 on a pitch circle of radius
// 0.7, grooves cut 0.01 deeper. The groove arc crosses each race
// wall at z = 0.2 +/- 0.1058.
let pitch = 0.7;
let ball_r = 0.15;
let groove_r = 0.16;
let lip = sqrt(groove_r * groove_r - 0.12 * 0.12);

// ---- Inner race: bore 0.35, wall at 0.58, groove to 0.54 ----
let inner = m.profile(xz())
    .path(0.35, 0.0).line_to(0.58, 0.0).line_to(0.58, 0.2 - lip)
    .arc_to([pitch - groove_r, 0.2], [0.58, 0.2 + lip])
    .line_to(0.58, 0.4).line_to(0.35, 0.4).close()
    .sketch("inner_section")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
inner = inner.chamfer(edges().circular().at_height(0.0).or(edges().circular().at_height(0.4)), 0.02,
    #{ label: "inner_lead_ins" });

// ---- Outer race: wall at 0.82, groove to 0.86, rim 1.05 ----
let outer = m.profile(xz())
    .path(0.82, 0.2 - lip).line_to(0.82, 0.0).line_to(1.05, 0.0)
    .line_to(1.05, 0.4).line_to(0.82, 0.4).line_to(0.82, 0.2 + lip)
    .arc_to([pitch + groove_r, 0.2], [0.82, 0.2 - lip]).close()
    .sketch("outer_section")
    .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
outer = outer.chamfer(edges().circular().at_height(0.0).or(edges().circular().at_height(0.4)), 0.02,
    #{ label: "outer_lead_ins" });

// ---- Eight balls, placed by the loop, not by hand ----
let pi = 3.141592653589793;
let count = 8;
for i in 0..count {
    let angle = i * 2.0 * pi / count;
    let ball = m.profile(xz())
        .path(0.0, 0.2 - ball_r).arc_to([ball_r, 0.2], [0.0, 0.2 + ball_r])
        .line_to(0.0, 0.2 - ball_r).close()
        .sketch("ball_" + i)
        .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
    ball.translate(pitch * cos(angle), pitch * sin(angle), 0.0);
}

m.assert_watertight();
m
