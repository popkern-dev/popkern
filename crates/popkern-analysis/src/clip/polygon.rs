//! Convex polygons clipped against axis-aligned planes.

/// Most vertices a triangle can have after clipping by six planes.
pub const MAX_VERTICES: usize = 9;

/// A convex planar polygon with at most [`MAX_VERTICES`] vertices.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Polygon {
    len: usize,
    vertices: [[f64; 3]; MAX_VERTICES],
}

impl Polygon {
    /// The polygon with no vertices.
    pub const fn empty() -> Polygon {
        Polygon {
            len: 0,
            vertices: [[0.0; 3]; MAX_VERTICES],
        }
    }

    /// The triangle `abc`.
    pub fn triangle(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Polygon {
        let mut vertices = [[0.0; 3]; MAX_VERTICES];
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
        Polygon { len: 3, vertices }
    }

    /// The vertices, in order.
    pub fn vertices(&self) -> &[[f64; 3]] {
        &self.vertices[..self.len]
    }

    /// Number of vertices.
    pub fn len(&self) -> usize {
        self.len
    }

    /// True when clipping removed every vertex.
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }

    fn push(&mut self, p: [f64; 3]) {
        // A convex polygon cut by one plane gains at most one vertex,
        // so six cuts of a triangle stay within MAX_VERTICES.
        if self.len < MAX_VERTICES {
            self.vertices[self.len] = p;
            self.len += 1;
        }
    }

    /// Area by the fan from the first vertex.
    pub fn area(&self) -> f64 {
        let v = self.vertices();
        let mut sum = [0.0; 3];
        for i in 1..v.len().saturating_sub(1) {
            let c = cross(sub(v[i], v[0]), sub(v[i + 1], v[0]));
            for k in 0..3 {
                sum[k] += c[k];
            }
        }
        0.5 * (sum[0] * sum[0] + sum[1] * sum[1] + sum[2] * sum[2]).sqrt()
    }

    /// True when every vertex has coordinate exactly `value` on `axis`.
    pub fn lies_on(&self, axis: usize, value: f64) -> bool {
        !self.is_empty() && self.vertices().iter().all(|v| v[axis] == value)
    }

    /// The part of the polygon on the kept side of the plane
    /// `p[axis] = value`: `p[axis] >= value` when `keep_above`, else
    /// `p[axis] <= value`. Points on the plane are kept.
    pub fn clip(&self, axis: usize, value: f64, keep_above: bool) -> Polygon {
        let inside = |p: &[f64; 3]| {
            if keep_above {
                p[axis] >= value
            } else {
                p[axis] <= value
            }
        };
        let mut out = Polygon::empty();
        let v = self.vertices();
        for i in 0..v.len() {
            let current = v[i];
            let next = v[(i + 1) % v.len()];
            let (ci, ni) = (inside(&current), inside(&next));
            if ci {
                out.push(current);
            }
            if ci != ni {
                out.push(crossing(current, next, axis, value));
            }
        }
        out
    }
}

/// Where segment pq crosses the plane `axis = value`, computed from the
/// endpoints in lexicographic order so both cells sharing the plane get
/// the same point, and snapped onto the plane.
fn crossing(p: [f64; 3], q: [f64; 3], axis: usize, value: f64) -> [f64; 3] {
    let (a, b) = if p <= q { (p, q) } else { (q, p) };
    let t = (value - a[axis]) / (b[axis] - a[axis]);
    let mut x = [0.0; 3];
    for k in 0..3 {
        x[k] = a[k] + t * (b[k] - a[k]);
    }
    x[axis] = value;
    x
}

/// a − b.
pub fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

/// a × b.
pub fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
