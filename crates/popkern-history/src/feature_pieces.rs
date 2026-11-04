//! Which pieces of a body a feature keeps (design: operations on the
//! cell complex, 5.4): a split's, and `KeepPieces` / `SeparatePieces`
//! on any body of several solids.

use serde::{Deserialize, Serialize};

use crate::failure::FeatureFailure;

/// Which pieces of a [`crate::feature::Feature::Split`] (or of a body,
/// for `KeepPieces` / `SeparatePieces`) are kept, recorded as data so a
/// replay with moved tools selects again.
#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub enum PieceRule {
    /// Every piece.
    #[default]
    All,
    /// The pieces whose membership in the split's tools (tool k as body
    /// k) the expression accepts.
    InsideTools(popkern_bool::SetExpr),
    /// The piece containing a point (model units).
    Containing([f64; 3]),
    /// The `n` largest pieces by volume.
    Largest(usize),
}

/// The pieces (by index into `pieces`) `rule` selects, in order.
/// `membership` gives each piece's tools and the tool count where the
/// pieces come from a split; without it `InsideTools` has no tools to
/// ask and is rejected.
///
/// # Errors
///
/// [`FeatureFailure::InvalidParameter`] for `InsideTools` without
/// tools; [`FeatureFailure::PieceRuleMatchedNothing`] when no piece is
/// selected.
pub(crate) fn select_pieces(
    feature: &str,
    rule: &PieceRule,
    pieces: &[&popkern_core::topology::BRep],
    membership: Option<(&[popkern_bool::BodySet], usize)>,
    tolerance: &popkern_core::tolerance::Tolerance,
) -> Result<Vec<usize>, Box<FeatureFailure>> {
    let count = pieces.len();
    let kept: Vec<usize> = match rule {
        PieceRule::All => (0..count).collect(),
        PieceRule::InsideTools(expr) => {
            let Some((inside, tools)) = membership else {
                return Err(FeatureFailure::InvalidParameter {
                    parameter: "rule".into(),
                    value: format!("{rule:?}"),
                    reason: format!(
                        "{feature} has no tools to be inside: select by tool on the split"
                    ),
                }
                .into());
            };
            (0..count)
                .filter(|&i| expr.holds(&inside[i], tools))
                .collect()
        }
        PieceRule::Containing([x, y, z]) => {
            let p = popkern_core::geometry::Point3::new(*x, *y, *z);
            (0..count)
                .filter(|&i| {
                    let boxes = popkern_bool::adaptive_precision::body_face_aabbs(pieces[i]);
                    popkern_bool::adaptive_precision::robust_point_in_solid(
                        &p, pieces[i], &boxes, tolerance,
                    ) == popkern_bool::PointClassification::Inside
                })
                .collect()
        }
        PieceRule::Largest(n) => {
            let mut by_volume: Vec<(usize, f64)> = (0..count)
                .map(|i| {
                    let v = popkern_measure::body_volume(pieces[i], tolerance).unwrap_or(0.0);
                    (i, v)
                })
                .collect();
            // Largest first; equal volumes keep the pieces' order.
            by_volume.sort_by(|a, b| b.1.total_cmp(&a.1));
            let mut kept: Vec<usize> = by_volume.into_iter().take(*n).map(|(i, _)| i).collect();
            kept.sort_unstable();
            kept
        }
    };
    if kept.is_empty() {
        return Err(FeatureFailure::PieceRuleMatchedNothing {
            feature: feature.into(),
            rule: format!("{rule:?}"),
            pieces: count,
        }
        .into());
    }
    Ok(kept)
}
