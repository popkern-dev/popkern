//! Small vector helpers for selector matching.

/// a − b.
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// a · b.
pub fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// a × b.
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean length (sqrt is correctly rounded; no portable wrapper).
pub fn norm(a: [f64; 3]) -> f64 {
    dot(a, a).sqrt()
}

/// `a` scaled to unit length, or `None` for a zero or non-finite vector.
pub fn unit(a: [f64; 3]) -> Option<[f64; 3]> {
    let n = norm(a);
    (n.is_finite() && n > 0.0).then(|| a.map(|c| c / n))
}

/// Distance from `p` to the line through `origin` along `axis`.
pub fn point_line_distance(p: [f64; 3], origin: [f64; 3], axis: [f64; 3]) -> f64 {
    match unit(axis) {
        Some(u) => norm(cross(sub(p, origin), u)),
        None => norm(sub(p, origin)),
    }
}

fn add_scaled(a: [f64; 3], b: [f64; 3], s: f64) -> [f64; 3] {
    [a[0] + s * b[0], a[1] + s * b[1], a[2] + s * b[2]]
}

/// Distance from `p` to triangle `abc` (closest point by Voronoi
/// regions; Ericson, *Real-Time Collision Detection*, §5.1.5).
pub fn point_triangle_distance(p: [f64; 3], a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    let ab = sub(b, a);
    let ac = sub(c, a);
    let ap = sub(p, a);
    let d1 = dot(ab, ap);
    let d2 = dot(ac, ap);
    if d1 <= 0.0 && d2 <= 0.0 {
        return norm(ap);
    }
    let bp = sub(p, b);
    let d3 = dot(ab, bp);
    let d4 = dot(ac, bp);
    if d3 >= 0.0 && d4 <= d3 {
        return norm(bp);
    }
    let vc = d1 * d4 - d3 * d2;
    if vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0 {
        let v = d1 / (d1 - d3);
        return norm(sub(p, add_scaled(a, ab, v)));
    }
    let cp = sub(p, c);
    let d5 = dot(ab, cp);
    let d6 = dot(ac, cp);
    if d6 >= 0.0 && d5 <= d6 {
        return norm(cp);
    }
    let vb = d5 * d2 - d1 * d6;
    if vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0 {
        let w = d2 / (d2 - d6);
        return norm(sub(p, add_scaled(a, ac, w)));
    }
    let va = d3 * d6 - d5 * d4;
    if va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0 {
        let w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return norm(sub(p, add_scaled(b, sub(c, b), w)));
    }
    let denom = va + vb + vc;
    if denom == 0.0 || !denom.is_finite() {
        // Degenerate triangle: fall back to the nearest corner.
        return norm(ap).min(norm(bp)).min(norm(cp));
    }
    let v = vb / denom;
    let w = vc / denom;
    norm(sub(p, add_scaled(add_scaled(a, ab, v), ac, w)))
}
