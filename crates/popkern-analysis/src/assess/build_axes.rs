//! A stress read against the layers of a print (notes 11,
//! orientation-aware rule): the normal stress across the layers, the
//! shear on them, and the in-layer block. Without a build direction
//! the worst case over every orientation is read from the principal
//! stresses instead.

use crate::recover::tensor::{Voigt, principal, von_mises};

/// The unit vector along `d`, `None` for a zero or non-finite one.
pub fn unit(d: [f64; 3]) -> Option<[f64; 3]> {
    let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
    (l.is_finite() && l > 0.0).then(|| d.map(|c| c / l))
}

/// σ·b for Voigt `s` (xx, yy, zz, yz, xz, xy; tensor shears).
fn apply(s: &Voigt, b: [f64; 3]) -> [f64; 3] {
    [
        s[0] * b[0] + s[5] * b[1] + s[4] * b[2],
        s[5] * b[0] + s[1] * b[1] + s[3] * b[2],
        s[4] * b[0] + s[3] * b[1] + s[2] * b[2],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// What the strength checks read at one stress sample, MPa.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Measures {
    /// Von Mises (of the in-layer block with a build direction).
    pub von_mises: f64,
    /// Largest principal stress, at least 0 (in-layer block likewise).
    pub tension: f64,
    /// Smallest principal stress, at most 0 (in-layer block likewise).
    pub compression: f64,
    /// The first invariant tr σ (of the in-layer block likewise).
    pub trace: f64,
    /// Tension across the layers, at least 0: (d·σ·d)⁺, or the largest
    /// principal at the worst-case orientation.
    pub across: f64,
    /// Shear on the layer planes: |σ·d − (d·σ·d) d|, or the largest
    /// shear on any plane, (σ₁ − σ₃)/2, at the worst-case orientation.
    pub layer_shear: f64,
}

/// The measures of `s` for unit build direction `d` (`None`: the
/// worst-case orientation).
pub fn measures(s: &Voigt, d: Option<[f64; 3]>) -> Measures {
    let Some(d) = d else {
        let p = principal(s);
        return Measures {
            von_mises: von_mises(s),
            tension: p[0].max(0.0),
            compression: p[2].min(0.0),
            trace: s[0] + s[1] + s[2],
            across: p[0].max(0.0),
            layer_shear: 0.5 * (p[0] - p[2]),
        };
    };
    let t = apply(s, d);
    let normal = dot(d, t);
    let shear = [0, 1, 2].map(|k| t[k] - normal * d[k]);
    // Any axis not parallel to d seeds the in-layer frame.
    let seed = if d[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let e1 = unit([0, 1, 2].map(|k| seed[k] - dot(seed, d) * d[k])).unwrap_or([1.0, 0.0, 0.0]);
    let e2 = [
        d[1] * e1[2] - d[2] * e1[1],
        d[2] * e1[0] - d[0] * e1[2],
        d[0] * e1[1] - d[1] * e1[0],
    ];
    let (a, b, c) = (
        dot(e1, apply(s, e1)),
        dot(e2, apply(s, e2)),
        dot(e1, apply(s, e2)),
    );
    let half = 0.5 * (a - b);
    let radius = (half * half + c * c).sqrt();
    Measures {
        von_mises: (a * a - a * b + b * b + 3.0 * c * c).max(0.0).sqrt(),
        tension: (0.5 * (a + b) + radius).max(0.0),
        compression: (0.5 * (a + b) - radius).min(0.0),
        trace: a + b,
        across: normal.max(0.0),
        layer_shear: dot(shear, shear).sqrt(),
    }
}
