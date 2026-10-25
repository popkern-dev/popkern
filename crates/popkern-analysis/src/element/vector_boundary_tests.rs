use super::*;
use crate::clip::FacetPiece;
use crate::clip::polygon::Polygon;
use crate::element::scalar::surface_moments;
use crate::grid::CELL_CORNERS;
use crate::model::FaceId;

/// For a linear field the flux Gram form is the traction's square times
/// the area (all of it for a fixed support, its normal part for a
/// roller).
#[test]
fn the_flux_gram_of_a_linear_field_is_its_traction_squared() {
    let (h, lo) = (0.5, [1.0, -0.5, 2.0]);
    let (a, b, c) = ([1.1, -0.4, 2.1], [1.45, -0.2, 2.2], [1.2, -0.05, 2.45]);
    let polygon = Polygon::triangle(a, b, c);
    let u = [0, 1, 2].map(|k| b[k] - a[k]);
    let w = [0, 1, 2].map(|k| c[k] - a[k]);
    let x = [
        u[1] * w[2] - u[2] * w[1],
        u[2] * w[0] - u[0] * w[2],
        u[0] * w[1] - u[1] * w[0],
    ];
    let len = (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]).sqrt();
    let n = x.map(|v| v / len);
    let piece = FacetPiece {
        cell: 0,
        triangle: 0,
        face: FaceId(0),
        area: 0.5 * len,
        polygon,
        normal: n,
    };
    let s = surface_moments(std::iter::once(&piece), lo, h);
    let (lambda, mu) = (1.7, 0.9);
    let g = [
        [1e-3, 2e-3, -4e-4],
        [-1e-3, 5e-4, 1e-3],
        [3e-4, 1e-3, -2e-3],
    ];
    let dofs: [f64; NODAL_DOFS] = std::array::from_fn(|r| {
        let corner = CELL_CORNERS[r / 3];
        let p: [f64; 3] = std::array::from_fn(|k| lo[k] + corner[k] as f64 * h);
        (0..3).map(|j| g[r % 3][j] * p[j]).sum::<f64>()
    });
    let trace = g[0][0] + g[1][1] + g[2][2];
    let traction: [f64; 3] = std::array::from_fn(|i| {
        (0..3)
            .map(|j| {
                let sigma = if i == j { lambda * trace } else { 0.0 } + mu * (g[i][j] + g[j][i]);
                sigma * n[j]
            })
            .sum()
    });
    let bt = BoundaryTables::new();
    let form = |t: &Matrix24| -> f64 {
        (0..NODAL_DOFS)
            .map(|r| {
                (0..NODAL_DOFS)
                    .map(|c| dofs[r] * t[r][c] * dofs[c])
                    .sum::<f64>()
            })
            .sum()
    };
    let fixed = PieceSupport::Fixed {
        u0: [0.0; 3],
        gradient: [[0.0; 3]; 3],
    };
    let full = form(&bt.flux_gram(fixed, &s, n, h, lambda, mu, None));
    let want = traction.iter().map(|t| t * t).sum::<f64>() * piece.area;
    assert!((full - want).abs() <= 1e-12 * want, "{full} vs {want}");
    let normal = form(&bt.flux_gram(PieceSupport::Roller, &s, n, h, lambda, mu, None));
    let tn: f64 = (0..3).map(|k| traction[k] * n[k]).sum();
    let want = tn * tn * piece.area;
    assert!((normal - want).abs() <= 1e-12 * want, "{normal} vs {want}");
}
