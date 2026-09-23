// A pentagonal pyramid carved from a blank by five exact tilted
// planes — no pyramid primitive needed. Its five slant edges get
// constant rounds, and where the five rounds meet at the apex the
// kernel builds the setback Gregory vertex cap (twenty patch faces)
// in place of the sharp point. Requesting g2 here measures honestly
// as a downgrade today: the cap interior is curvature-continuous but
// the strip-to-cap rim rows are an open frontier (PK-515).
let m = model();

let pi = 3.141592653589793;
let sides = 5;
let base_r = 1.0;          // base circumradius
let height = 1.2;
let apothem = base_r * cos(pi / sides);

// The blank: a block taller and wider than the finished pyramid.
let pyramid = m.profile(xy()).rect(3.0, 3.0).sketch("blank").extrude(blind(1.6));

// Each lateral face is one slab intersection. The face plane holds
// the apex (0, 0, height) and one base edge; its outward normal is
// height * u + apothem * z for the face's outward direction u.
let deep = 4.0;
for i in 0..sides {
    let ang = 2.0 * pi * i / sides;
    let u = [cos(ang), sin(ang)];
    let len = sqrt(height * height + apothem * apothem);
    let n = [height * u[0] / len, height * u[1] / len, apothem / len];
    let d = apothem * height / len;   // plane offset: n . apex
    let p0 = [(d - deep) * n[0], (d - deep) * n[1], (d - deep) * n[2]];
    let slab = m.tool(|t| t.profile(plane_at(p0, n))
        .rect(8.0, 8.0).sketch("face_slab_" + i)
        .extrude(blind(deep)));
    pyramid = pyramid.intersect(slab, #{ label: "face_cut_" + i });
}

// The five slant edges are the straight edges of slant length; the
// base ring (a shorter chord) stays sharp.
let slant = sqrt(base_r * base_r + height * height);
pyramid.fillet(
    edges().straight().chord_length_between(slant - 0.05, slant + 0.05),
    0.16, #{ continuity: "g1", label: "apex_rounds" });

m.assert_watertight();
m
