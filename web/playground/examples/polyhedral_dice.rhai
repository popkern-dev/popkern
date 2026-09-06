// The five Platonic solids as a dice set — d4, d6, d8, d12 and d20
// in a row, every face an exact plane. The face normals are computed
// in this script (the golden ratio builds the d12 and d20), and each
// die is carved from a block by intersecting one slab per face pair.
// The d6 is filleted and dimpled with 21 subtracted spherical pips,
// opposite faces summing to seven like a real die.
let m = model();

let phi = (1.0 + sqrt(5.0)) / 2.0;

fn norm(v) {
    let l = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    [v[0] / l, v[1] / l, v[2] / l]
}

// Carve a die centred at the origin: keep the material within
// inradius r of each (unit) face normal. Paired normals cut both the
// face and its antipode with one centred slab; the tetrahedron's
// faces have no antipodes, so its slabs reach from far below.
fn carve(m, tag, normals, r, paired, big) {
    let die = m.profile(xy_at(-big)).rect(2.0 * big, 2.0 * big)
        .sketch(tag + "_block").extrude(blind(2.0 * big));
    let i = 0;
    for raw in normals {
        let n = norm(raw);
        let slab = if paired {
            m.tool(|t| t.profile(plane_at([-r * n[0], -r * n[1], -r * n[2]], n))
                .rect(4.0 * big, 4.0 * big).sketch(tag + "_slab_" + i)
                .extrude(blind(2.0 * r)))
        } else {
            let deep = 3.0 * big;
            m.tool(|t| t.profile(plane_at([-deep * n[0], -deep * n[1], -deep * n[2]], n))
                .rect(6.0 * big, 6.0 * big).sketch(tag + "_slab_" + i)
                .extrude(blind(deep + r)))
        };
        die = die.intersect(slab, #{ label: tag + "_cut_" + i });
        i += 1;
    }
    die
}

// Rest a carved die on the table: rotate the chosen face normal to
// straight down, then set it on z = 0 at its row position.
fn rest(m, die, tag, down, r, x, y) {
    let pi = 3.141592653589793;
    let d = norm(down);
    let c = -d[2];
    if c >= 1.0 { c = 1.0 }
    if c <= -1.0 { c = -1.0 }
    let ang = acos(c) * 180.0 / pi;
    // axis = d x (0,0,-1)
    let axis = [-d[1], d[0], 0.0];
    let al = sqrt(axis[0] * axis[0] + axis[1] * axis[1]);
    if al < 1.0e-12 { axis = [1.0, 0.0, 0.0]; }
    if ang > 1.0e-9 {
        die = die.rotate(axis_at([0.0, 0.0, 0.0], axis), ang, #{ label: tag + "_rest" });
    }
    die.translate(x, y, r, #{ label: tag + "_place" })
}

// ---- d4: tetrahedron, vertices at four corners of a cube ----
let r4 = 0.18;
let d4 = carve(m, "d4",
    [[-1.0, -1.0, -1.0], [-1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [1.0, 1.0, -1.0]],
    r4, false, 3.2 * r4);
d4 = rest(m, d4, "d4", [-1.0, -1.0, -1.0], r4, -0.8, -0.1);

// ---- d8: octahedron, faces toward the cube corners ----
let r8 = 0.26;
let d8 = carve(m, "d8",
    [[1.0, 1.0, 1.0], [1.0, -1.0, 1.0], [-1.0, 1.0, 1.0], [-1.0, -1.0, 1.0]],
    r8, true, 1.9 * r8);
d8 = rest(m, d8, "d8", [-1.0, -1.0, -1.0], r8, 0.62, -0.55);

// ---- d12: dodecahedron, faces toward the icosahedron vertices ----
let r12 = 0.28;
let d12 = carve(m, "d12",
    [[0.0, 1.0, phi], [0.0, -1.0, phi], [1.0, phi, 0.0],
     [1.0, -phi, 0.0], [phi, 0.0, 1.0], [phi, 0.0, -1.0]],
    r12, true, 1.5 * r12);
d12 = rest(m, d12, "d12", [0.0, -1.0, -phi], r12, -0.15, 0.5);

// ---- d20: icosahedron, faces toward the dodecahedron vertices ----
let r20 = 0.30;
let d20 = carve(m, "d20",
    [[1.0, 1.0, 1.0], [1.0, 1.0, -1.0], [1.0, -1.0, 1.0], [1.0, -1.0, -1.0],
     [0.0, 1.0 / phi, phi], [0.0, -1.0 / phi, phi], [1.0 / phi, phi, 0.0],
     [-1.0 / phi, phi, 0.0], [phi, 0.0, 1.0 / phi], [phi, 0.0, -1.0 / phi]],
    r20, true, 1.5 * r20);
d20 = rest(m, d20, "d20", [-1.0, -1.0, -1.0], r20, 0.95, 0.5);

// ---- d6: the classic pip die, built resting and filleted ----
let s6 = 0.55;
let h6 = s6 / 2.0;
let d6 = m.profile(xy_at(-h6)).rect(s6, s6).sketch("d6_block").extrude(blind(s6));
d6 = d6.fillet(edges().straight(), 0.05, #{ label: "d6_rounds" });

// One pip: a sphere (a semicircle revolved about z, as the bearing
// builds its balls), moved to the pip centre just outside the face so
// its cap dips in as a dimple. A sphere is rotation-invariant, so one
// revolve axis serves every face.
fn pip(m, tag, center, rp) {
    m.tool(|t| {
        let ball = t.profile(xz())
            .path(0.0, -rp).arc_to([rp, 0.0], [0.0, rp]).line_to(0.0, -rp)
            .close().sketch(tag)
            .revolve(axis_at([0.0, 0.0, 0.0], [0.0, 0.0, 1.0]), full_circle());
        ball.translate(center[0], center[1], center[2])
    })
}

// Pip layouts on a face's (u, v) grid; opposite faces sum to seven.
fn pips_for(value, g) {
    switch value {
        1 => [[0.0, 0.0]],
        2 => [[-1.0 * g, -1.0 * g], [g, g]],
        3 => [[-1.0 * g, -1.0 * g], [0.0, 0.0], [g, g]],
        4 => [[-1.0 * g, -1.0 * g], [-1.0 * g, g], [g, -1.0 * g], [g, g]],
        5 => [[-1.0 * g, -1.0 * g], [-1.0 * g, g], [0.0, 0.0], [g, -1.0 * g], [g, g]],
        _ => [[-1.0 * g, -1.0 * g], [-1.0 * g, 0.0], [-1.0 * g, g],
              [g, -1.0 * g], [g, 0.0], [g, g]],
    }
}

let rp = 0.072;          // pip sphere radius
let sink = 0.052;        // sphere centre this far outside the face
let g = 0.15;            // pip grid pitch
// Each face: [normal, u axis, v axis, value].
let d6_faces = [
    [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1],
    [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, -1.0, 0.0], 6],
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], 3],
    [[-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], 4],
    [[0.0, 1.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 2],
    [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0], 5],
];
for face in d6_faces {
    let n = face[0];
    let u = face[1];
    let v = face[2];
    let value = face[3];
    let k = 0;
    for uv in pips_for(value, g) {
        let c = [
            n[0] * (h6 + sink) + u[0] * uv[0] + v[0] * uv[1],
            n[1] * (h6 + sink) + u[1] * uv[0] + v[1] * uv[1],
            n[2] * (h6 + sink) + u[2] * uv[0] + v[2] * uv[1],
        ];
        let tag = "d6_pip_" + value + "_" + k;
        d6 = d6.subtract(pip(m, tag, c, rp), #{ label: tag + "_cut" });
        k += 1;
    }
}
d6 = d6.translate(-0.1, -0.42, h6, #{ label: "d6_place" });

m.assert_watertight();
m
