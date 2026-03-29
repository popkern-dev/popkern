//! Geometry helper functions for pattern operations.
//!
//! Rotation, reflection, and synchronisation of cloned arc/conic
//! parameters. Split from `sketch_pattern_geo.rs` to stay within the
//! 500-line file limit.

use crate::elements::{ArcHandle, ElementHandle, PointHandle};
use crate::Sketch;

/// Rotate a point around a centre.
pub(crate) fn rotate_point(
    px: f64,
    py: f64,
    cx: f64,
    cy: f64,
    cos_a: f64,
    sin_a: f64,
) -> (f64, f64) {
    let dx = px - cx;
    let dy = py - cy;
    (cx + dx * cos_a - dy * sin_a, cy + dx * sin_a + dy * cos_a)
}

/// Reflect a point across a line defined by two points.
pub(crate) fn reflect_point(
    px: f64,
    py: f64,
    lx1: f64,
    ly1: f64,
    lx2: f64,
    ly2: f64,
) -> (f64, f64) {
    let dx = lx2 - lx1;
    let dy = ly2 - ly1;
    let len2 = dx * dx + dy * dy;
    if len2 < 1e-30 {
        return (px, py);
    }
    let vx = px - lx1;
    let vy = py - ly1;
    let t = (vx * dx + vy * dy) / len2;
    let fx = lx1 + t * dx;
    let fy = ly1 + t * dy;
    (2.0 * fx - px, 2.0 * fy - py)
}

/// Set a point's parameter values directly.
pub(crate) fn set_point_params(sketch: &mut Sketch, ph: PointHandle, x: f64, y: f64) {
    if let Some(p) = sketch.points.get(ph) {
        let [ix, iy] = p.param_indices();
        sketch.params[ix] = x;
        sketch.params[iy] = y;
    }
}

/// Rotate the centres of cloned arcs/conics to match rotated positions.
#[allow(clippy::too_many_arguments)]
pub(crate) fn rotate_cloned_centres(
    sketch: &mut Sketch,
    seed: &[ElementHandle],
    remap: &crate::pattern::HandleRemap,
    cx: f64,
    cy: f64,
    cos_a: f64,
    sin_a: f64,
) {
    for elem in seed {
        match elem {
            ElementHandle::Arc(ah) => {
                if let Some(&new_ah) = remap.arcs.get(ah) {
                    rotate_centre_2(sketch, *ah, new_ah, cx, cy, cos_a, sin_a, true);
                }
            }
            ElementHandle::Ellipse(eh) => {
                if let Some(&new_eh) = remap.ellipses.get(eh) {
                    let Some(oe) = sketch.ellipses.get(*eh) else {
                        continue;
                    };
                    let os = oe.param_start;
                    let (rx, ry) = rotate_point(
                        sketch.params[os],
                        sketch.params[os + 1],
                        cx,
                        cy,
                        cos_a,
                        sin_a,
                    );
                    let Some(ce) = sketch.ellipses.get(new_eh) else {
                        continue;
                    };
                    let cs = ce.param_start;
                    sketch.params[cs] = rx;
                    sketch.params[cs + 1] = ry;
                }
            }
            ElementHandle::Parabola(ph) => {
                if let Some(&new_ph) = remap.parabolas.get(ph) {
                    let Some(op) = sketch.parabolas.get(*ph) else {
                        continue;
                    };
                    let os = op.param_start;
                    let (rx, ry) = rotate_point(
                        sketch.params[os],
                        sketch.params[os + 1],
                        cx,
                        cy,
                        cos_a,
                        sin_a,
                    );
                    let Some(cp) = sketch.parabolas.get(new_ph) else {
                        continue;
                    };
                    let cs = cp.param_start;
                    sketch.params[cs] = rx;
                    sketch.params[cs + 1] = ry;
                }
            }
            ElementHandle::Hyperbola(hh) => {
                if let Some(&new_hh) = remap.hyperbolas.get(hh) {
                    let Some(oh) = sketch.hyperbolas.get(*hh) else {
                        continue;
                    };
                    let os = oh.param_start;
                    let (rx, ry) = rotate_point(
                        sketch.params[os],
                        sketch.params[os + 1],
                        cx,
                        cy,
                        cos_a,
                        sin_a,
                    );
                    let Some(ch) = sketch.hyperbolas.get(new_hh) else {
                        continue;
                    };
                    let cs = ch.param_start;
                    sketch.params[cs] = rx;
                    sketch.params[cs + 1] = ry;
                }
            }
            _ => {}
        }
    }
}

/// Rotate an arc centre.
#[allow(clippy::too_many_arguments)]
fn rotate_centre_2(
    sketch: &mut Sketch,
    orig: ArcHandle,
    cloned: ArcHandle,
    cx: f64,
    cy: f64,
    cos_a: f64,
    sin_a: f64,
    _is_arc: bool,
) {
    let Some(oa) = sketch.arcs.get(orig) else {
        return;
    };
    let os = oa.param_start;
    let (rx, ry) = rotate_point(
        sketch.params[os],
        sketch.params[os + 1],
        cx,
        cy,
        cos_a,
        sin_a,
    );
    let Some(ca) = sketch.arcs.get(cloned) else {
        return;
    };
    let cs = ca.param_start;
    sketch.params[cs] = rx;
    sketch.params[cs + 1] = ry;
}

/// Reflect centres of cloned arcs/conics.
#[allow(clippy::too_many_arguments)]
pub(crate) fn reflect_cloned_centres(
    sketch: &mut Sketch,
    seed: &[ElementHandle],
    remap: &crate::pattern::HandleRemap,
    lx1: f64,
    ly1: f64,
    lx2: f64,
    ly2: f64,
) {
    let line_angle = (ly2 - ly1).atan2(lx2 - lx1);
    for elem in seed {
        match elem {
            ElementHandle::Arc(ah) => {
                if let Some(&new_ah) = remap.arcs.get(ah) {
                    reflect_arc(sketch, *ah, new_ah, lx1, ly1, lx2, ly2, line_angle);
                }
            }
            ElementHandle::Ellipse(eh) => {
                if let Some(&new_eh) = remap.ellipses.get(eh) {
                    reflect_5param(sketch, *eh, new_eh, lx1, ly1, lx2, ly2, line_angle, 4);
                }
            }
            ElementHandle::Parabola(ph) => {
                if let Some(&new_ph) = remap.parabolas.get(ph) {
                    reflect_4param(sketch, *ph, new_ph, lx1, ly1, lx2, ly2, line_angle);
                }
            }
            ElementHandle::Hyperbola(hh) => {
                if let Some(&new_hh) = remap.hyperbolas.get(hh) {
                    reflect_5param_hyp(sketch, *hh, new_hh, lx1, ly1, lx2, ly2, line_angle);
                }
            }
            _ => {}
        }
    }
}

/// Reflect arc centre and swap angles.
#[allow(clippy::too_many_arguments)]
fn reflect_arc(
    sketch: &mut Sketch,
    orig: ArcHandle,
    cloned: ArcHandle,
    lx1: f64,
    ly1: f64,
    lx2: f64,
    ly2: f64,
    line_angle: f64,
) {
    let Some(oa) = sketch.arcs.get(orig) else {
        return;
    };
    let os = oa.param_start;
    let (rcx, rcy) = reflect_point(sketch.params[os], sketch.params[os + 1], lx1, ly1, lx2, ly2);
    let rsa = 2.0 * line_angle - sketch.params[os + 4];
    let rea = 2.0 * line_angle - sketch.params[os + 3];
    let Some(ca) = sketch.arcs.get(cloned) else {
        return;
    };
    let cs = ca.param_start;
    sketch.params[cs] = rcx;
    sketch.params[cs + 1] = rcy;
    sketch.params[cs + 3] = rsa;
    sketch.params[cs + 4] = rea;
}

/// Reflect ellipse centre and rotation (rotation at offset 4).
#[allow(clippy::too_many_arguments)]
fn reflect_5param(
    sketch: &mut Sketch,
    orig: crate::EllipseHandle,
    cloned: crate::EllipseHandle,
    lx1: f64,
    ly1: f64,
    lx2: f64,
    ly2: f64,
    line_angle: f64,
    rot_offset: usize,
) {
    let Some(oe) = sketch.ellipses.get(orig) else {
        return;
    };
    let os = oe.param_start;
    let (rcx, rcy) = reflect_point(sketch.params[os], sketch.params[os + 1], lx1, ly1, lx2, ly2);
    let rrot = 2.0 * line_angle - sketch.params[os + rot_offset];
    let Some(ce) = sketch.ellipses.get(cloned) else {
        return;
    };
    let cs = ce.param_start;
    sketch.params[cs] = rcx;
    sketch.params[cs + 1] = rcy;
    sketch.params[cs + rot_offset] = rrot;
}

/// Reflect parabola vertex and rotation (rotation at offset 3).
#[allow(clippy::too_many_arguments)]
fn reflect_4param(
    sketch: &mut Sketch,
    orig: crate::ParabolaHandle,
    cloned: crate::ParabolaHandle,
    lx1: f64,
    ly1: f64,
    lx2: f64,
    ly2: f64,
    line_angle: f64,
) {
    let Some(op) = sketch.parabolas.get(orig) else {
        return;
    };
    let os = op.param_start;
    let (rcx, rcy) = reflect_point(sketch.params[os], sketch.params[os + 1], lx1, ly1, lx2, ly2);
    let rrot = 2.0 * line_angle - sketch.params[os + 3];
    let Some(cp) = sketch.parabolas.get(cloned) else {
        return;
    };
    let cs = cp.param_start;
    sketch.params[cs] = rcx;
    sketch.params[cs + 1] = rcy;
    sketch.params[cs + 3] = rrot;
}

/// Reflect hyperbola centre and rotation (rotation at offset 4).
#[allow(clippy::too_many_arguments)]
fn reflect_5param_hyp(
    sketch: &mut Sketch,
    orig: crate::HyperbolaHandle,
    cloned: crate::HyperbolaHandle,
    lx1: f64,
    ly1: f64,
    lx2: f64,
    ly2: f64,
    line_angle: f64,
) {
    let Some(oh) = sketch.hyperbolas.get(orig) else {
        return;
    };
    let os = oh.param_start;
    let (rcx, rcy) = reflect_point(sketch.params[os], sketch.params[os + 1], lx1, ly1, lx2, ly2);
    let rrot = 2.0 * line_angle - sketch.params[os + 4];
    let Some(ch) = sketch.hyperbolas.get(cloned) else {
        return;
    };
    let cs = ch.param_start;
    sketch.params[cs] = rcx;
    sketch.params[cs + 1] = rcy;
    sketch.params[cs + 4] = rrot;
}

/// Sync phantom point positions to match cloned arc/conic centres.
pub(crate) fn sync_phantom_positions(
    sketch: &mut Sketch,
    seed: &[ElementHandle],
    remap: &crate::pattern::HandleRemap,
    phantom_points: &[PointHandle],
) {
    let mut idx = 0;
    for elem in seed {
        let center_start = match elem {
            ElementHandle::Arc(ah) => remap
                .arcs
                .get(ah)
                .and_then(|&na| sketch.arcs.get(na).map(|a| a.param_start)),
            ElementHandle::Ellipse(eh) => remap
                .ellipses
                .get(eh)
                .and_then(|&ne| sketch.ellipses.get(ne).map(|e| e.param_start)),
            ElementHandle::Parabola(ph) => remap
                .parabolas
                .get(ph)
                .and_then(|&np| sketch.parabolas.get(np).map(|p| p.param_start)),
            ElementHandle::Hyperbola(hh) => remap
                .hyperbolas
                .get(hh)
                .and_then(|&nh| sketch.hyperbolas.get(nh).map(|h| h.param_start)),
            _ => None,
        };
        if let Some(cs) = center_start {
            if idx < phantom_points.len() {
                let cx = sketch.params[cs];
                let cy = sketch.params[cs + 1];
                set_point_params(sketch, phantom_points[idx], cx, cy);
                idx += 1;
            }
        }
    }
}
