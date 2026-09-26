//! Enclosure of a span's partial nets over a parameter box (PK-562).
//!
//! Every partial net is restricted to the box across u, then across v,
//! and hulled — except along an iso line: a seam on a patch boundary
//! holds u (or v) at one exact value on every interval, so that
//! direction is evaluated once per span into a cached line per net
//! and each interval only restricts the line across the other
//! direction (measured on the G2 pyramid/boss corpus: 92% of
//! certification time was cap-cap seams on bi-nonic patch boundaries).
//! Both orders enclose rigorously; the cached order evaluates the
//! point direction first.

use std::collections::HashMap;

use super::super::bernstein::{self, Coef};
use super::super::interval::Iv;
use super::Net;

/// `(span_u, span_v, bits of the local box ends in the point direction)`.
type IsoKey = (usize, usize, [u64; 4]);

/// Iso lines cached per source and direction before the cache resets
/// (a seam holds one or two; the bound only caps memory).
const ISO_LINES_CACHED: usize = 16;

/// Reused buffers for restricting span nets to a box (PK-562: the
/// per-interval restriction of 28 partial nets allocated ~7 vectors
/// per row and column).
#[derive(Debug, Clone, Default)]
struct BoxScratch {
    restrict: bernstein::RestrictScratch<4>,
    line: Vec<Coef<4>>,
    restricted: Vec<Coef<4>>,
    /// Column-major `[column][restricted u index]`.
    across_u: Vec<Coef<4>>,
}

/// Box enclosure state of one [`super::NurbsJetSource`].
#[derive(Debug, Clone, Default)]
pub(super) struct SpanBoxEncloser {
    scratch: BoxScratch,
    /// Per span and exact point `u` box: every partial net's columns
    /// evaluated there (one row per net). A seam along a `u = const`
    /// patch boundary meets the same point `u` at every interval, so
    /// only the row restriction across `v` is per-interval work.
    iso_u_rows: HashMap<IsoKey, Vec<Vec<Coef<4>>>>,
    /// The same for an exact point `v` box: every partial net's rows
    /// evaluated there (one column per net), restricted across `u`.
    iso_v_columns: HashMap<IsoKey, Vec<Vec<Coef<4>>>>,
}

impl SpanBoxEncloser {
    /// Hull of every partial net of span `span` restricted to the box
    /// with local ends `u = (au, bu)`, `v = (av, bv)`, indexed like the
    /// nets.
    pub(super) fn enclose(
        &mut self,
        span: (usize, usize),
        nets: &[Net],
        (au, bu): (Iv, Iv),
        (av, bv): (Iv, Iv),
    ) -> Vec<Coef<4>> {
        let (su, sv) = span;
        let s = &mut self.scratch;
        // A seam along a patch boundary holds u (or v) at one exact
        // value on every interval: that direction is evaluated once per
        // span into a cached line per partial net, and each interval
        // only restricts the line across the other direction. Caches
        // fill from genuine iso boxes; a point box (an interval
        // midpoint) only looks them up.
        let u_point = bu.hi - au.lo <= bernstein::POINT_BOX_WIDTH;
        let v_point = bv.hi - av.lo <= bernstein::POINT_BOX_WIDTH;
        let u_key = (su, sv, [au.lo, au.hi, bu.lo, bu.hi].map(f64::to_bits));
        let v_key = (su, sv, [av.lo, av.hi, bv.lo, bv.hi].map(f64::to_bits));
        if u_point && !v_point && !self.iso_u_rows.contains_key(&u_key) {
            if self.iso_u_rows.len() >= ISO_LINES_CACHED {
                self.iso_u_rows.clear();
            }
            let rows = nets
                .iter()
                .map(|net| {
                    let mut row = Vec::new();
                    restrict_columns(net, au, bu, s, &mut row);
                    row
                })
                .collect();
            self.iso_u_rows.insert(u_key, rows);
        }
        if v_point && !u_point && !self.iso_v_columns.contains_key(&v_key) {
            if self.iso_v_columns.len() >= ISO_LINES_CACHED {
                self.iso_v_columns.clear();
            }
            let columns = nets
                .iter()
                .map(|net| {
                    let mut column = Vec::new();
                    // Empty nets (orders past the degree) stay empty.
                    for row in net.iter().filter(|row| !row.is_empty()) {
                        bernstein::restrict_into(row, av, bv, &mut s.restrict, &mut s.restricted);
                        column.extend_from_slice(&s.restricted);
                    }
                    column
                })
                .collect();
            self.iso_v_columns.insert(v_key, columns);
        }
        let iso_u = u_point.then(|| self.iso_u_rows.get(&u_key)).flatten();
        let iso_v = v_point.then(|| self.iso_v_columns.get(&v_key)).flatten();
        let homogeneous: Vec<Coef<4>> = nets
            .iter()
            .enumerate()
            .map(|(n, net)| {
                if net.is_empty() || net[0].is_empty() {
                    return [Iv::point(0.0); 4];
                }
                let mut hull: Option<Coef<4>> = None;
                let mut absorb = |coefficients: &[Coef<4>]| {
                    for c in coefficients {
                        hull = Some(match hull {
                            None => *c,
                            Some(mut acc) => {
                                for k in 0..4 {
                                    acc[k] = acc[k].hull(c[k]);
                                }
                                acc
                            }
                        });
                    }
                };
                if let Some(row) = iso_u.and_then(|rows| rows.get(n)) {
                    // Point u: the columns' values at u, across v.
                    bernstein::restrict_into(row, av, bv, &mut s.restrict, &mut s.restricted);
                    absorb(&s.restricted);
                } else if let Some(column) = iso_v.and_then(|columns| columns.get(n)) {
                    // Point v: the rows' values at v, across u.
                    bernstein::restrict_into(column, au, bu, &mut s.restrict, &mut s.restricted);
                    absorb(&s.restricted);
                } else {
                    // Each column across u, then each resulting row
                    // across v, hulled in row-major order.
                    let columns = net[0].len();
                    let mut across_u = std::mem::take(&mut s.across_u);
                    let rows = restrict_columns(net, au, bu, s, &mut across_u);
                    for i in 0..rows {
                        s.line.clear();
                        s.line.extend((0..columns).map(|j| across_u[j * rows + i]));
                        bernstein::restrict_into(
                            &s.line,
                            av,
                            bv,
                            &mut s.restrict,
                            &mut s.restricted,
                        );
                        absorb(&s.restricted);
                    }
                    s.across_u = across_u;
                }
                hull.unwrap_or([Iv::point(0.0); 4])
            })
            .collect();
        homogeneous
    }
}

/// Every column `net[·][j]` restricted across u to `[alpha, beta]`,
/// written column-major into `out`; returns the restricted length.
fn restrict_columns(
    net: &[Vec<Coef<4>>],
    alpha: Iv,
    beta: Iv,
    s: &mut BoxScratch,
    out: &mut Vec<Coef<4>>,
) -> usize {
    out.clear();
    let mut rows = 0;
    for j in 0..net.first().map_or(0, Vec::len) {
        s.line.clear();
        s.line.extend(net.iter().map(|row| row[j]));
        bernstein::restrict_into(&s.line, alpha, beta, &mut s.restrict, &mut s.restricted);
        rows = s.restricted.len();
        out.extend_from_slice(&s.restricted);
    }
    rows
}
