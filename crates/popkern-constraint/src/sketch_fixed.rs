//! [`Sketch::capture_fixed_values`] — split from `lib.rs` for the line limit.

use crate::constraint::Constraint;
use crate::Sketch;

impl Sketch {
    /// Snapshot the current parameter values that a fixed constraint needs
    /// to lock in place.
    ///
    /// Returns an empty vector for all non-fixed constraint kinds.
    pub fn capture_fixed_values(&self, constraint: &Constraint) -> Vec<f64> {
        match constraint {
            Constraint::Fixed(ph) => {
                if let Some(p) = self.points.get(*ph) {
                    let [xi, yi] = p.param_indices();
                    vec![self.params[xi], self.params[yi]]
                } else {
                    vec![]
                }
            }
            Constraint::FixedLine(lh) => {
                if let Some(l) = self.lines.get(*lh) {
                    let mut vals = Vec::with_capacity(4);
                    for ph in [l.start, l.end] {
                        if let Some(p) = self.points.get(ph) {
                            let [xi, yi] = p.param_indices();
                            vals.push(self.params[xi]);
                            vals.push(self.params[yi]);
                        }
                    }
                    vals
                } else {
                    vec![]
                }
            }
            Constraint::FixedEllipse(eh) => {
                if let Some(e) = self.ellipses.get(*eh) {
                    let idx = e.param_indices();
                    idx.iter().map(|&i| self.params[i]).collect()
                } else {
                    vec![]
                }
            }
            Constraint::FixedParabola(ph) => {
                if let Some(p) = self.parabolas.get(*ph) {
                    let idx = p.param_indices();
                    idx.iter().map(|&i| self.params[i]).collect()
                } else {
                    vec![]
                }
            }
            Constraint::FixedHyperbola(hh) => {
                if let Some(h) = self.hyperbolas.get(*hh) {
                    let idx = h.param_indices();
                    idx.iter().map(|&i| self.params[i]).collect()
                } else {
                    vec![]
                }
            }
            Constraint::SplineFixedEndTangent(_s, _end_idx, tx, ty) => {
                vec![*tx, *ty]
            }
            _ => vec![],
        }
    }
}
