fn ball(m, r) {
    m.sphere(r)
}

// Two balls of radius 5, on the origin and on (9.999999, 0, 0), overlapping 1e-6 deep: their
// common part.
let m = model();
let lens = ball(m, 5.0);
let other = m.tool(|m| ball(m, 5.0).translate(9.999999, 0.0, 0.0));
lens = lens.intersect(other);
m
