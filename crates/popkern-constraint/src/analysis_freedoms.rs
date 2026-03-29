//! Helper functions for mapping null-space vectors and zero-coverage
//! parameters back to sketch elements as [`FreedomInfo`] records.
//!
//! Split from the main `analysis` module to stay within the 500-line limit.

use crate::elements::ElementHandle;
use crate::error::FreedomInfo;
use crate::Sketch;
use faer::Mat;

/// Map all parameters to freedom entries (used when there are no equations).
pub(super) fn params_to_freedoms(sketch: &Sketch, _n: usize) -> Vec<FreedomInfo> {
    let mut freedoms = Vec::new();
    for (handle, point) in sketch.points.iter() {
        let _ = point;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Point(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, arc) in sketch.arcs.iter() {
        let _ = arc;
        // Arcs have 5 DOF but we report as element-level freedom.
        freedoms.push(FreedomInfo {
            element: ElementHandle::Arc(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, ellipse) in sketch.ellipses.iter() {
        let _ = ellipse;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Ellipse(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, parabola) in sketch.parabolas.iter() {
        let _ = parabola;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Parabola(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, hyperbola) in sketch.hyperbolas.iter() {
        let _ = hyperbola;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Hyperbola(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    freedoms
}

/// Extract freedom info from a single null-space vector column of V.
pub(super) fn add_freedom_from_vector(
    sketch: &Sketch,
    v: faer::MatRef<'_, f64>,
    col_idx: usize,
    threshold: f64,
    freedoms: &mut Vec<FreedomInfo>,
) {
    // Check which elements have significant components in this null-space
    // vector.
    for (handle, point) in sketch.points.iter() {
        let [ix, iy] = point.param_indices();
        let vx = if ix < v.nrows() {
            v.read(ix, col_idx)
        } else {
            0.0
        };
        let vy = if iy < v.nrows() {
            v.read(iy, col_idx)
        } else {
            0.0
        };
        let mag = (vx * vx + vy * vy).sqrt();
        if mag > threshold {
            let inv = 1.0 / mag;
            freedoms.push(FreedomInfo {
                element: ElementHandle::Point(handle),
                directions: vec![(vx * inv, vy * inv)],
            });
        }
    }
    for (handle, arc) in sketch.arcs.iter() {
        let [icx, icy, _ir, _is, _ie] = arc.param_indices();
        let vcx = if icx < v.nrows() {
            v.read(icx, col_idx)
        } else {
            0.0
        };
        let vcy = if icy < v.nrows() {
            v.read(icy, col_idx)
        } else {
            0.0
        };
        let mag = (vcx * vcx + vcy * vcy).sqrt();
        if mag > threshold {
            let inv = 1.0 / mag;
            freedoms.push(FreedomInfo {
                element: ElementHandle::Arc(handle),
                directions: vec![(vcx * inv, vcy * inv)],
            });
        }
    }
    add_conic_freedom_from_vector(sketch, v, col_idx, threshold, freedoms);
}

/// Extract conic element freedom info from a null-space vector column.
fn add_conic_freedom_from_vector(
    sketch: &Sketch,
    v: faer::MatRef<'_, f64>,
    col_idx: usize,
    threshold: f64,
    freedoms: &mut Vec<FreedomInfo>,
) {
    for (handle, ellipse) in sketch.ellipses.iter() {
        let [icx, icy, _, _, _] = ellipse.param_indices();
        let vcx = if icx < v.nrows() {
            v.read(icx, col_idx)
        } else {
            0.0
        };
        let vcy = if icy < v.nrows() {
            v.read(icy, col_idx)
        } else {
            0.0
        };
        let mag = (vcx * vcx + vcy * vcy).sqrt();
        if mag > threshold {
            let inv = 1.0 / mag;
            freedoms.push(FreedomInfo {
                element: ElementHandle::Ellipse(handle),
                directions: vec![(vcx * inv, vcy * inv)],
            });
        }
    }
    for (handle, parabola) in sketch.parabolas.iter() {
        let [ivx, ivy, _, _] = parabola.param_indices();
        let vvx = if ivx < v.nrows() {
            v.read(ivx, col_idx)
        } else {
            0.0
        };
        let vvy = if ivy < v.nrows() {
            v.read(ivy, col_idx)
        } else {
            0.0
        };
        let mag = (vvx * vvx + vvy * vvy).sqrt();
        if mag > threshold {
            let inv = 1.0 / mag;
            freedoms.push(FreedomInfo {
                element: ElementHandle::Parabola(handle),
                directions: vec![(vvx * inv, vvy * inv)],
            });
        }
    }
    for (handle, hyperbola) in sketch.hyperbolas.iter() {
        let [icx, icy, _, _, _] = hyperbola.param_indices();
        let vcx = if icx < v.nrows() {
            v.read(icx, col_idx)
        } else {
            0.0
        };
        let vcy = if icy < v.nrows() {
            v.read(icy, col_idx)
        } else {
            0.0
        };
        let mag = (vcx * vcx + vcy * vcy).sqrt();
        if mag > threshold {
            let inv = 1.0 / mag;
            freedoms.push(FreedomInfo {
                element: ElementHandle::Hyperbola(handle),
                directions: vec![(vcx * inv, vcy * inv)],
            });
        }
    }
}

/// Add freedoms for parameters that have zero Jacobian coverage (implicit
/// null-space directions when n > min(m, n)).
pub(super) fn add_implicit_freedoms(
    sketch: &Sketch,
    j_mat: &Mat<f64>,
    n: usize,
    freedoms: &mut Vec<FreedomInfo>,
) {
    // Find parameters with all-zero Jacobian columns.
    let mut col_has_entry = vec![false; n];
    for (col, has_entry) in col_has_entry.iter_mut().enumerate() {
        for row in 0..j_mat.nrows() {
            if j_mat.read(row, col).abs() > 1e-15 {
                *has_entry = true;
                break;
            }
        }
    }

    for (handle, point) in sketch.points.iter() {
        let [ix, iy] = point.param_indices();
        let free_x = ix < n && !col_has_entry[ix];
        let free_y = iy < n && !col_has_entry[iy];
        if free_x || free_y {
            let mut dirs = Vec::new();
            if free_x {
                dirs.push((1.0, 0.0));
            }
            if free_y {
                dirs.push((0.0, 1.0));
            }
            freedoms.push(FreedomInfo {
                element: ElementHandle::Point(handle),
                directions: dirs,
            });
        }
    }
    for (handle, arc) in sketch.arcs.iter() {
        let [icx, icy, _ir, _is, _ie] = arc.param_indices();
        let free_cx = icx < n && !col_has_entry[icx];
        let free_cy = icy < n && !col_has_entry[icy];
        if free_cx || free_cy {
            let mut dirs = Vec::new();
            if free_cx {
                dirs.push((1.0, 0.0));
            }
            if free_cy {
                dirs.push((0.0, 1.0));
            }
            freedoms.push(FreedomInfo {
                element: ElementHandle::Arc(handle),
                directions: dirs,
            });
        }
    }
    add_conic_implicit_freedoms(sketch, &col_has_entry, n, freedoms);
}

/// Add implicit freedoms for conic elements.
fn add_conic_implicit_freedoms(
    sketch: &Sketch,
    col_has_entry: &[bool],
    n: usize,
    freedoms: &mut Vec<FreedomInfo>,
) {
    for (handle, ellipse) in sketch.ellipses.iter() {
        let [icx, icy, _, _, _] = ellipse.param_indices();
        let free_cx = icx < n && !col_has_entry[icx];
        let free_cy = icy < n && !col_has_entry[icy];
        if free_cx || free_cy {
            let mut dirs = Vec::new();
            if free_cx {
                dirs.push((1.0, 0.0));
            }
            if free_cy {
                dirs.push((0.0, 1.0));
            }
            freedoms.push(FreedomInfo {
                element: ElementHandle::Ellipse(handle),
                directions: dirs,
            });
        }
    }
    for (handle, parabola) in sketch.parabolas.iter() {
        let [ivx, ivy, _, _] = parabola.param_indices();
        let free_vx = ivx < n && !col_has_entry[ivx];
        let free_vy = ivy < n && !col_has_entry[ivy];
        if free_vx || free_vy {
            let mut dirs = Vec::new();
            if free_vx {
                dirs.push((1.0, 0.0));
            }
            if free_vy {
                dirs.push((0.0, 1.0));
            }
            freedoms.push(FreedomInfo {
                element: ElementHandle::Parabola(handle),
                directions: dirs,
            });
        }
    }
    for (handle, hyperbola) in sketch.hyperbolas.iter() {
        let [icx, icy, _, _, _] = hyperbola.param_indices();
        let free_cx = icx < n && !col_has_entry[icx];
        let free_cy = icy < n && !col_has_entry[icy];
        if free_cx || free_cy {
            let mut dirs = Vec::new();
            if free_cx {
                dirs.push((1.0, 0.0));
            }
            if free_cy {
                dirs.push((0.0, 1.0));
            }
            freedoms.push(FreedomInfo {
                element: ElementHandle::Hyperbola(handle),
                directions: dirs,
            });
        }
    }
}

/// Merge [`FreedomInfo`] entries that reference the same element.
pub(super) fn merge_freedoms(freedoms: &mut Vec<FreedomInfo>) {
    if freedoms.len() <= 1 {
        return;
    }
    let mut merged: Vec<FreedomInfo> = Vec::new();
    for f in freedoms.drain(..) {
        if let Some(existing) = merged.iter_mut().find(|e| e.element == f.element) {
            for dir in f.directions {
                // Only add if not a near-duplicate direction.
                let is_dup = existing.directions.iter().any(|&(dx, dy)| {
                    let dot = (dx * dir.0 + dy * dir.1).abs();
                    dot > 0.99
                });
                if !is_dup {
                    existing.directions.push(dir);
                }
            }
        } else {
            merged.push(f);
        }
    }
    *freedoms = merged;
}
