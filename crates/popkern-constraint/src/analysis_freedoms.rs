//! Helper functions for mapping null-space vectors and zero-coverage
//! parameters back to sketch elements as [`FreedomInfo`] records.
//!
//! Split from the main `analysis` module to stay within the 500-line limit.

use crate::elements::ElementHandle;
use crate::error::FreedomInfo;
use crate::rigid_set;
use crate::Sketch;
use faer::Mat;

/// Map all parameters to freedom entries (used when there are no equations).
pub(super) fn params_to_freedoms(sketch: &Sketch, _n: usize) -> Vec<FreedomInfo> {
    let mut freedoms = Vec::new();
    let rigid_params = rigid_set::rigid_set_point_params(sketch);

    for (handle, point) in sketch.points.iter() {
        if point.rigid_set.is_some() {
            continue; // reported via the rigid set
        }
        let _ = point;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Point(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, arc) in sketch.arcs.iter() {
        if arc.rigid_set.is_some() {
            continue;
        }
        let _ = arc;
        // Arcs have 5 DOF but we report as element-level freedom.
        freedoms.push(FreedomInfo {
            element: ElementHandle::Arc(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, ellipse) in sketch.ellipses.iter() {
        if ellipse.rigid_set.is_some() {
            continue;
        }
        let _ = ellipse;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Ellipse(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, parabola) in sketch.parabolas.iter() {
        if parabola.rigid_set.is_some() {
            continue;
        }
        let _ = parabola;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Parabola(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    for (handle, hyperbola) in sketch.hyperbolas.iter() {
        if hyperbola.rigid_set.is_some() {
            continue;
        }
        let _ = hyperbola;
        freedoms.push(FreedomInfo {
            element: ElementHandle::Hyperbola(handle),
            directions: vec![(1.0, 0.0), (0.0, 1.0)],
        });
    }
    // Report rigid-set freedoms as the first element in the set.
    for rs in sketch.rigid_sets.values() {
        if let Some(&first) = rs.elements.first() {
            let _ = &rigid_params; // keep the binding alive for clarity
                                   // Rigid set has 3 DOF: tx, ty, theta.
            freedoms.push(FreedomInfo {
                element: first,
                directions: vec![(1.0, 0.0), (0.0, 1.0)],
            });
        }
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
    let rigid_params = rigid_set::rigid_set_point_params(sketch);

    // Check which elements have significant components in this null-space
    // vector.
    for (handle, point) in sketch.points.iter() {
        if point.rigid_set.is_some() {
            continue; // skip rigid-set members
        }
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
        if arc.rigid_set.is_some() {
            continue;
        }
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

    // Report rigid-set freedoms based on their (tx, ty, theta) params.
    add_rigid_set_freedoms(sketch, v, col_idx, threshold, &rigid_params, freedoms);
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
        if ellipse.rigid_set.is_some() {
            continue;
        }
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
        if parabola.rigid_set.is_some() {
            continue;
        }
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
        if hyperbola.rigid_set.is_some() {
            continue;
        }
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

/// Report rigid-set DOF as freedoms on the first element.
fn add_rigid_set_freedoms(
    sketch: &Sketch,
    v: faer::MatRef<'_, f64>,
    col_idx: usize,
    threshold: f64,
    _rigid_params: &std::collections::HashSet<usize>,
    freedoms: &mut Vec<FreedomInfo>,
) {
    for rs in sketch.rigid_sets.values() {
        let first = match rs.elements.first() {
            Some(e) => *e,
            None => continue,
        };
        // Check the tx, ty columns of V.
        let itx = rs.param_start;
        let ity = rs.param_start + 1;
        let vtx = if itx < v.nrows() {
            v.read(itx, col_idx)
        } else {
            0.0
        };
        let vty = if ity < v.nrows() {
            v.read(ity, col_idx)
        } else {
            0.0
        };
        let mag = (vtx * vtx + vty * vty).sqrt();
        if mag > threshold {
            let inv = 1.0 / mag;
            freedoms.push(FreedomInfo {
                element: first,
                directions: vec![(vtx * inv, vty * inv)],
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
        if point.rigid_set.is_some() {
            continue;
        }
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
        if arc.rigid_set.is_some() {
            continue;
        }
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

    // Report rigid-set implicit freedoms.
    for rs in sketch.rigid_sets.values() {
        let first = match rs.elements.first() {
            Some(e) => *e,
            None => continue,
        };
        let itx = rs.param_start;
        let ity = rs.param_start + 1;
        let ith = rs.param_start + 2;
        let free_tx = itx < n && !col_has_entry[itx];
        let free_ty = ity < n && !col_has_entry[ity];
        let free_th = ith < n && !col_has_entry[ith];
        if free_tx || free_ty || free_th {
            let mut dirs = Vec::new();
            if free_tx {
                dirs.push((1.0, 0.0));
            }
            if free_ty {
                dirs.push((0.0, 1.0));
            }
            if free_th {
                // Theta is a rotational DOF; represent as (0, 0) or
                // skip — for simplicity report as translational proxy.
                dirs.push((1.0, 0.0));
            }
            freedoms.push(FreedomInfo {
                element: first,
                directions: dirs,
            });
        }
    }
}

/// Add implicit freedoms for conic elements.
fn add_conic_implicit_freedoms(
    sketch: &Sketch,
    col_has_entry: &[bool],
    n: usize,
    freedoms: &mut Vec<FreedomInfo>,
) {
    for (handle, ellipse) in sketch.ellipses.iter() {
        if ellipse.rigid_set.is_some() {
            continue;
        }
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
        if parabola.rigid_set.is_some() {
            continue;
        }
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
        if hyperbola.rigid_set.is_some() {
            continue;
        }
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
