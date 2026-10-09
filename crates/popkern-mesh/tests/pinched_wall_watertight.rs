//! A face whose hole lobes touch at pinch vertices, beside a face that
//! straddles its cylinder's seam, tessellates watertight (PK-732): the
//! equal-radius tee with its branch end inside the main cylinder. The
//! main wall's hole is a lens with a lobe touching it at each singular
//! point (0, ±40, 60), and a lobe of the branch wall runs across the
//! branch cylinder's seam at negative u.
//!
//! - The boundary-first mesher cut the wall's ring at each pinch, but
//!   the cocircular jitter moved the two copies of the pinch point
//!   apart and the lobes' first segments crossed; held together, the
//!   triangulation still kept one copy and lost the segments at the
//!   other. Either way the wall fell back to a mesh off the canonical
//!   edge samples.
//! - Meshed face by face (the trim-grid path), the seam-straddling
//!   lobe's grid was clamped at u = 0, and a periodic face's grid was
//!   cut to the principal branch.
//!
//! Measured before: 893 open mesh edges.

use popkern_core::geometry::Transform3;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;

#[test]
fn the_equal_radius_tee_with_pinched_holes_meshes_watertight() {
    let tol = Tolerance::default();
    // The branch cylinder's seam along its bottom, then through the
    // singular point S+ = (0, 40, 60): there a lobe of the branch wall
    // straddles the seam.
    for seam_turn in [0.0, std::f64::consts::FRAC_PI_2] {
        let main = BRep::cylinder(40.0, 120.0, tol).expect("main");
        // The branch along +x from x = −20, its axis at z = 60.
        let branch = BRep::cylinder(40.0, 120.0, tol)
            .expect("branch")
            .transform(&Transform3::rotation_z(seam_turn))
            .expect("seam")
            .transform(&Transform3::rotation_y(std::f64::consts::FRAC_PI_2))
            .expect("turn")
            .transform(&Transform3::translation(-20.0, 0.0, 60.0))
            .expect("place");
        let tee = popkern_bool::union(main, branch, tol, None).expect("union");
        let mesh = popkern_mesh::tessellate(&tee, &popkern_mesh::TessellationOptions::default())
            .expect("tessellation");
        let open = popkern_mesh::welded_open_edges(&mesh);
        assert_eq!(open, 0, "seam turn {seam_turn}: open mesh edges");
        // Mesh volume against the analytic union (chord-height accuracy).
        let volume: f64 = mesh
            .triangles
            .iter()
            .map(|t| {
                let [a, b, c] = t.map(|i| mesh.vertices[i].position);
                (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                    + a.z * (b.x * c.y - b.y * c.x))
                    / 6.0
            })
            .sum();
        let want = 936_581.667;
        assert!(
            (volume.abs() - want).abs() / want < 2e-3,
            "seam turn {seam_turn}: mesh volume {volume} vs {want}"
        );
        // Each face meshed on its own (the trim-grid path) covers what
        // the whole-body mesh does: the seam-straddling lobe's grid was
        // clamped at u = 0 and covered half its area.
        let area = |triangles: &[[usize; 3]],
                    at: &dyn Fn(usize) -> popkern_core::geometry::Point3| {
            triangles
                .iter()
                .map(|t| {
                    let [a, b, c] = t.map(at);
                    (b - a).cross(&(c - a)).length() / 2.0
                })
                .sum::<f64>()
        };
        let whole = area(&mesh.triangles, &|i| mesh.vertices[i].position);
        let per_face: f64 = tee
            .faces
            .keys()
            .map(|f| {
                let fm = popkern_mesh::tessellate_face(
                    &tee,
                    f,
                    &popkern_mesh::TessellationOptions::default(),
                )
                .expect("face mesh");
                area(&fm.triangles, &|i| fm.vertices[i].position)
            })
            .sum();
        assert!(
            (per_face - whole).abs() / whole < 1e-3,
            "seam turn {seam_turn}: per-face area {per_face} vs whole {whole}"
        );
    }
}
