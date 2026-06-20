pub(super) fn unit_domain_loop_candidates(edge_count: usize) -> Vec<Vec<(f64, f64)>> {
    let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let bases = if edge_count == 4 {
        vec![corners.to_vec()]
    } else if edge_count == 5 {
        (0..4)
            .map(|split| {
                let mut out = Vec::with_capacity(5);
                for index in 0..4 {
                    out.push(corners[index]);
                    if index == split {
                        let next = corners[(index + 1) % 4];
                        out.push((
                            (corners[index].0 + next.0) * 0.5,
                            (corners[index].1 + next.1) * 0.5,
                        ));
                    }
                }
                out
            })
            .collect()
    } else if edge_count > 5 {
        boundary_aligned_unit_square_loops(edge_count)
    } else {
        Vec::new()
    };
    oriented_loop_candidates(bases)
}

fn oriented_loop_candidates(bases: Vec<Vec<(f64, f64)>>) -> Vec<Vec<(f64, f64)>> {
    let mut candidates = Vec::new();
    for base in bases {
        for shift in 0..base.len() {
            let rotated = (0..base.len())
                .map(|index| base[(index + shift) % base.len()])
                .collect::<Vec<_>>();
            candidates.push(rotated.clone());
            candidates.push(rotated.into_iter().rev().collect());
        }
    }
    candidates
}

fn boundary_aligned_unit_square_loops(edge_count: usize) -> Vec<Vec<(f64, f64)>> {
    let extra_count = edge_count.saturating_sub(4);
    let mut loops = Vec::new();
    for bottom_extra in 0..=extra_count {
        for right_extra in 0..=extra_count - bottom_extra {
            for top_extra in 0..=extra_count - bottom_extra - right_extra {
                let left_extra = extra_count - bottom_extra - right_extra - top_extra;
                loops.push(boundary_aligned_unit_square_loop([
                    bottom_extra,
                    right_extra,
                    top_extra,
                    left_extra,
                ]));
            }
        }
    }
    loops
}

fn boundary_aligned_unit_square_loop(side_extra_counts: [usize; 4]) -> Vec<(f64, f64)> {
    let corners = [(0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)];
    let mut out = Vec::with_capacity(4 + side_extra_counts.iter().sum::<usize>());
    for side in 0..4 {
        let start = corners[side];
        let end = corners[(side + 1) % 4];
        out.push(start);
        let extra = side_extra_counts[side];
        for split in 0..extra {
            let t = (split + 1) as f64 / (extra + 1) as f64;
            out.push(lerp_uv(start, end, t));
        }
    }
    out
}

fn lerp_uv(a: (f64, f64), b: (f64, f64), t: f64) -> (f64, f64) {
    (a.0 + (b.0 - a.0) * t, a.1 + (b.1 - a.1) * t)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn high_valence_slots_have_boundary_aligned_unit_domain_candidates() {
        let candidates = unit_domain_loop_candidates(6);

        assert!(!candidates.is_empty());
        assert!(candidates.iter().all(|candidate| candidate.len() == 6));
        assert!(candidates.iter().flatten().all(|&(u, v)| {
            (0.0..=1.0).contains(&u)
                && (0.0..=1.0).contains(&v)
                && (u == 0.0 || u == 1.0 || v == 0.0 || v == 1.0)
        }));
    }
}
