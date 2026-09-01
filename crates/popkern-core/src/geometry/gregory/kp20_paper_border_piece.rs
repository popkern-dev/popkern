//! The KP20 paper border piece: one sector half at native degree 6
//! along the rim (PK-398).

use super::kp20_border_degree::{elevate_bi3_rows_to_sextic, reduce_sextic_rows_to_bi3};
use crate::geometry::Point3;

/// One border piece: rows `[depth][along]`, depth 0 on the cap
/// boundary, degree 6 along (seven coefficients) and three depth
/// rows (the bi-3 depth convention). PK-398: the along degree is
/// native 6 so the setback fan's exact strip section reaches the
/// band machinery directly; cubic sources elevate exactly through
/// [`Self::from_bi3`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Kp20PaperGuideBorderPiece {
    /// Degree-6 BB rows, `[depth][along]`.
    pub rows: [[Point3; 7]; 3],
}

impl Kp20PaperGuideBorderPiece {
    /// Exact elevation of a bi-3 piece.
    pub fn from_bi3(rows: [[Point3; 4]; 3]) -> Self {
        Self {
            rows: elevate_bi3_rows_to_sextic(rows),
        }
    }

    /// End-pinned least-squares cubic reduction of the rows (exact
    /// for elevated cubics); for consumers that keep the bi-3 form.
    pub fn bi3_rows(&self) -> [[Point3; 4]; 3] {
        reduce_sextic_rows_to_bi3(&self.rows)
    }

    /// Cubic reduction that reproduces this piece's 2-jet (value,
    /// first and second derivative) at one end exactly — `at_end`
    /// selects the far end — with the other endpoint interpolated.
    /// Two pieces meeting C² at a spoke point reduce to cubics that
    /// still meet C² there (the KP15 backstop's eq (2)/(3) spoke
    /// gates consume exactly that), where independent least-squares
    /// reductions would not.
    pub fn bi3_rows_matching_end_jet(&self, at_end: bool) -> [[Point3; 4]; 3] {
        let reduce = |row: &[Point3; 7]| -> [Point3; 4] {
            let (a, b, c, far) = if at_end {
                (row[6], row[5], row[4], row[0])
            } else {
                (row[0], row[1], row[2], row[6])
            };
            // Degree-6 end jet: f' = 6(b − a), f'' = 30(c − 2b + a);
            // cubic coefficients: p1 = p0 + f'/3, p2 = 2p1 − p0 + f''/6.
            let first = (b - a) * 2.0;
            let second = ((c - b) - (b - a)) * 5.0;
            let p1 = a + first;
            let p2 = p1 + first + second;
            if at_end {
                [far, p2, p1, a]
            } else {
                [a, p1, p2, far]
            }
        };
        [
            reduce(&self.rows[0]),
            reduce(&self.rows[1]),
            reduce(&self.rows[2]),
        ]
    }
}
