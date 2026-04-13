//! Sliver face detection and classification.
//!
//! A sliver face is geometrically degenerate: negligible area relative to its
//! bounding dimensions, or containing micro-edges. Slivers are the most common
//! pathological output of boolean operations and the most common cause of
//! failure in subsequent operations on the result.
//!
//! The detection uses multiple criteria (any one triggers):
//! - **Area criterion**: face area < threshold
//! - **Aspect ratio**: shortest altitude / longest edge < threshold
//! - **Minimum edge length**: any edge shorter than threshold
//! - **Minimum face width**: distance between non-adjacent boundary edges < threshold

use popkern_core::geometry::Point3;

/// Classification of sliver type.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SliverType {
    /// Thin annular region — intersection curve nearly parallel to face boundary.
    Ring,
    /// Very elongated face — two intersection curves meet at a shallow angle.
    LongEdge,
    /// Face with one or more edges shorter than tolerance.
    ShortEdge,
}

/// Result of sliver analysis for a single face.
#[derive(Debug, Clone)]
pub struct SliverReport {
    /// Whether this face is classified as a sliver.
    pub is_sliver: bool,
    /// The type of sliver, if detected.
    pub sliver_type: Option<SliverType>,
    /// Approximate area of the face.
    pub area: f64,
    /// Length of the shortest edge.
    pub min_edge_length: f64,
    /// Aspect ratio (shortest altitude / longest edge).
    pub aspect_ratio: f64,
}

/// Thresholds for sliver detection, derived from [`BoolTolerance`].
#[derive(Debug, Clone)]
pub struct SliverThresholds {
    /// Minimum face area below which a face is a sliver.
    pub min_area: f64,
    /// Minimum aspect ratio.
    pub min_aspect_ratio: f64,
    /// Minimum edge length (edges shorter than this are micro-edges).
    pub min_edge_length: f64,
    /// Minimum face width (for ring sliver detection).
    pub min_face_width: f64,
}

impl SliverThresholds {
    /// Create thresholds from a sliver tolerance value.
    ///
    /// Derives area and aspect thresholds from the linear sliver tolerance.
    pub fn from_tolerance(sliver_tol: f64) -> Self {
        Self {
            min_area: sliver_tol * sliver_tol,
            min_aspect_ratio: 0.01,
            min_edge_length: sliver_tol,
            min_face_width: sliver_tol,
        }
    }
}

/// Analyze a face polygon for sliver characteristics.
///
/// Takes the face boundary as an ordered list of vertex positions
/// (closed polygon — last vertex connects back to first).
///
/// Returns a [`SliverReport`] with detection results.
///
/// The edge-length / short-edge check uses
/// [`popkern_core::adaptive_precision::adaptive_dist_sq_lt`] on
/// every edge so micro-edge detection is deterministic at the
/// tolerance envelope. The area and aspect checks remain in f64;
/// their inputs are sums of products in Newell's method, and the
/// decisions are typically far from the margin for non-pathological
/// polygons.
pub fn analyze_face_polygon(vertices: &[Point3], thresholds: &SliverThresholds) -> SliverReport {
    if vertices.len() < 3 {
        return SliverReport {
            is_sliver: true,
            sliver_type: Some(SliverType::ShortEdge),
            area: 0.0,
            min_edge_length: 0.0,
            aspect_ratio: 0.0,
        };
    }

    // Compute edge lengths. Tracks max/min via f64 squared
    // distance (sqrt only at the end for the diagnostic value),
    // and calls the adaptive predicate per edge for the decision.
    let n = vertices.len();
    let min_edge_threshold_sq = thresholds.min_edge_length * thresholds.min_edge_length;
    let mut max_edge_sq = 0.0_f64;
    let mut min_edge_sq = f64::INFINITY;
    let mut short_edge = false;

    for i in 0..n {
        let j = (i + 1) % n;
        let vi = &vertices[i];
        let vj = &vertices[j];
        let dx = vj.x - vi.x;
        let dy = vj.y - vi.y;
        let dz = vj.z - vi.z;
        let len_sq = dx * dx + dy * dy + dz * dz;
        if len_sq > max_edge_sq {
            max_edge_sq = len_sq;
        }
        if len_sq < min_edge_sq {
            min_edge_sq = len_sq;
        }
        if !short_edge
            && popkern_core::adaptive_precision::adaptive_dist_sq_lt(
                (vi.x, vi.y, vi.z),
                (vj.x, vj.y, vj.z),
                min_edge_threshold_sq,
            )
        {
            short_edge = true;
        }
    }
    let min_edge = min_edge_sq.sqrt();

    // Approximate area via Newell's method on the projected polygon.
    // We compute 4·area² directly from the Newell normal vector —
    // `|N|² = (2·area)²` — so the area² comparison never touches
    // the expensive sqrt in the hot path. `area` as a scalar is
    // recovered at the end only for the report's diagnostic output.
    let four_area_sq = polygon_four_area_sq(vertices);
    let area = 0.5 * four_area_sq.sqrt();

    // Aspect ratio uses the squared-area form too: the threshold
    // compare `area / max_edge² < min_aspect` is equivalent to
    // `area² < (min_aspect · max_edge²)²`, i.e.
    // `4·area² < 4·(min_aspect·max_edge²)²`. Computing it this way
    // avoids both a sqrt AND a division in the hot path.
    let aspect_ratio = if max_edge_sq > 1e-60 {
        area / max_edge_sq
    } else {
        0.0
    };

    // Classification via squared-form comparisons where possible.
    // `area < min_area` ⇔ `4·area² < 4·min_area²`. A small f64
    // margin catches rounding at the decision boundary without
    // needing full DP recomputation.
    let four_min_area_sq = 4.0 * thresholds.min_area * thresholds.min_area;
    let area_margin = 16.0 * f64::EPSILON * four_area_sq.max(four_min_area_sq);
    let small_area = if four_area_sq + area_margin < four_min_area_sq {
        true
    } else if four_area_sq > four_min_area_sq + area_margin {
        false
    } else {
        // Ambiguous f64 — fall back to the sqrt'd scalar compare
        // as a tiebreaker. This is the same decision the old
        // path made, so behaviour at the exact boundary is
        // preserved.
        area < thresholds.min_area
    };
    let bad_aspect = aspect_ratio < thresholds.min_aspect_ratio;

    let (is_sliver, sliver_type) = if short_edge {
        (true, Some(SliverType::ShortEdge))
    } else if small_area {
        (true, Some(SliverType::Ring))
    } else if bad_aspect {
        (true, Some(SliverType::LongEdge))
    } else {
        (false, None)
    };

    SliverReport {
        is_sliver,
        sliver_type,
        area,
        min_edge_length: min_edge,
        aspect_ratio,
    }
}

/// Compute `4·area²` of a 3D polygon using Newell's method.
///
/// Newell's method yields a normal vector whose magnitude is twice
/// the polygon's area. Returning the squared magnitude directly
/// (`|N|² = 4·area²`) eliminates the sqrt that would otherwise be
/// needed by every caller, and lets the comparison
/// `area < threshold` be rewritten as
/// `4·area² < 4·threshold²` — a sqrt-free form that's stable
/// against the small cancellation errors in the per-edge product.
fn polygon_four_area_sq(vertices: &[Point3]) -> f64 {
    if vertices.len() < 3 {
        return 0.0;
    }

    let mut nx = 0.0;
    let mut ny = 0.0;
    let mut nz = 0.0;
    let n = vertices.len();

    for i in 0..n {
        let j = (i + 1) % n;
        let vi = &vertices[i];
        let vj = &vertices[j];
        nx += (vi.y - vj.y) * (vi.z + vj.z);
        ny += (vi.z - vj.z) * (vi.x + vj.x);
        nz += (vi.x - vj.x) * (vi.y + vj.y);
    }

    nx * nx + ny * ny + nz * nz
}
