//! Offset-campaign fixture corpus (PK-323, design doc §9).
//!
//! Ten bodies spanning the behaviors the offset/shell campaign must
//! master: G1 fillet joins, bores and countersinks, concave edges,
//! collapse-width ribs, free-form NURBS, global self-intersection
//! under modest offsets (spindled torus and thin-wall collision),
//! open rims, mixed analytic junctions, and high-valence corners. Each fixture carries known-correct topology
//! counts (derived, not transcribed from output), an analytic volume
//! where one is derivable, and Steiner/tube constants where the body
//! is convex analytic (plus the torus, whose tube formula holds below
//! the reach) — the constants feed the SDF-oracle calibration battery.

use popkern_blend::{FilletOptions, FilletRadius, FilletSpec};
use popkern_core::error::CoreError;
use popkern_core::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use popkern_core::geometry::transform::Transform3;
use popkern_core::geometry::{NurbsSurface, Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_part::{Part, PartError};
use std::f64::consts::PI;

/// Closed-form Steiner/tube constants for
/// [`popkern_mesh::offset_oracle::steiner_offset_volume`].
#[derive(Debug, Clone, Copy)]
pub struct SteinerConstants {
    /// Body volume.
    pub volume: f64,
    /// Boundary surface area.
    pub surface_area: f64,
    /// Total mean-curvature integral `∫H dA + Σ (θₑ/2)·ℓₑ`.
    pub mean_curvature_integral: f64,
    /// Euler characteristic of the boundary (2 sphere-like, 0 torus).
    pub euler_characteristic: f64,
}

/// One corpus body with its ground-truth expectations.
#[derive(Debug, Clone)]
pub struct CorpusFixture {
    /// Stable fixture name (used in reports and STEP exports).
    pub name: &'static str,
    /// The body.
    pub brep: BRep,
    /// Whether the fixture bounds a closed solid (the open sheet does
    /// not — watertightness checks are waived for it, nothing else).
    pub watertight: bool,
    /// Expected face count.
    pub expected_faces: usize,
    /// Expected edge count.
    pub expected_edges: usize,
    /// Expected vertex count.
    pub expected_vertices: usize,
    /// Closed-form volume, where derivable.
    pub analytic_volume: Option<f64>,
    /// Relative band for the volume pin (tessellation + construction
    /// approximation error, fixture-specific).
    pub volume_band: f64,
    /// Steiner/tube constants where the closed form exists.
    pub steiner: Option<SteinerConstants>,
    /// Expected face count per surface kind, sorted by kind name.
    pub expected_census: Vec<(&'static str, usize)>,
}

/// Corpus construction failure.
#[derive(Debug, thiserror::Error)]
pub enum CorpusError {
    /// A kernel constructor or Euler operation failed.
    #[error("core construction failed: {0}")]
    Core(#[from] CoreError),
    /// A Part-facade operation (boolean, serving validation) failed.
    #[error("part operation failed: {0}")]
    Part(#[from] PartError),
    /// A fillet failed (the filleted-cube fixture requires full
    /// success, not partial).
    #[error("fillet failed: {0}")]
    Fillet(String),
}

fn part(brep: BRep, tol: Tolerance, label: &str) -> Result<Part, CorpusError> {
    Ok(Part::from_brep(brep, tol, label)?)
}

/// Count the faces of `brep` per surface kind, sorted by kind name —
/// the corpus's census pins compare against this.
pub fn surface_census(brep: &BRep) -> Vec<(&'static str, usize)> {
    let mut counts = std::collections::BTreeMap::new();
    for (_, face) in brep.faces.iter() {
        let kind = match face.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => "plane",
            Some(SurfaceGeometry::Cylinder { .. }) => "cylinder",
            Some(SurfaceGeometry::Cone { .. }) => "cone",
            Some(SurfaceGeometry::Sphere { .. }) => "sphere",
            Some(SurfaceGeometry::Torus { .. }) => "torus",
            Some(SurfaceGeometry::Nurbs(_)) => "nurbs",
            Some(SurfaceGeometry::Gregory(_)) => "gregory",
            Some(SurfaceGeometry::PolyhedralSpline(_)) => "polyhedral_spline",
            None => "untyped",
        };
        *counts.entry(kind).or_insert(0) += 1;
    }
    counts.into_iter().collect()
}

/// Cube `[0,1]³` with all 12 edges filleted at `r = 0.1`: cylindrical
/// edge strips and corner caps promoted to exact sphere octants
/// (corner-sphere promotion), so the body IS the rounded cube and its
/// closed forms are exact up to tessellation error.
pub fn filleted_cube(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let rf = 0.1;
    let base = BRep::unit_cube(tol)?;
    let specs: Vec<FilletSpec> = base
        .edges
        .keys()
        .map(|edge| FilletSpec {
            edge,
            radius: FilletRadius::Constant(rf),
        })
        .collect();
    let blend = popkern_blend::fillet_edges(&base, &specs, &FilletOptions::default(), tol, None)
        .map_err(|e| CorpusError::Fillet(format!("{e:?}")))?;
    // Rounded-cube closed forms with core side c = 1 − 2·rf:
    // V = c³ + 6c²·rf + 3πc·rf² + (4π/3)rf³, A = 6c² + 6πc·rf + 4πrf²,
    // M = 3πc + 4π·rf (cube edge term shrunk to the core skeleton,
    // plus the full sphere contribution of the corners).
    let c = 1.0 - 2.0 * rf;
    let volume =
        c.powi(3) + 6.0 * c * c * rf + 3.0 * PI * c * rf * rf + 4.0 * PI / 3.0 * rf.powi(3);
    let area = 6.0 * c * c + 6.0 * PI * c * rf + 4.0 * PI * rf * rf;
    let m = 3.0 * PI * c + 4.0 * PI * rf;
    Ok(CorpusFixture {
        name: "filleted_cube",
        brep: blend.brep,
        watertight: true,
        // 6 shrunk squares + 12 quarter-cylinder strips + 8 corner
        // caps (one patch per corner). Edges: 24 square↔strip + 24
        // strip↔cap; vertices: 3 per cap, unshared between caps.
        // Euler: 24 − 48 + 26 = 2.
        expected_faces: 26,
        expected_edges: 48,
        expected_vertices: 24,
        analytic_volume: Some(volume),
        volume_band: 0.01,
        steiner: Some(SteinerConstants {
            volume,
            surface_area: area,
            mean_curvature_integral: m,
            euler_characteristic: 2.0,
        }),
        expected_census: vec![("cylinder", 12), ("plane", 6), ("sphere", 8)],
    })
}

/// Block `2×1×0.5` with a through-bore (r = 0.2) and a 45° countersink
/// meeting the bore at `z = 0.4`. Exercises holes, coaxial
/// cylinder×cone junctions, and inner loops.
pub fn bracket_bore_countersink(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(2.0, 0.0, 0.0),
        Point3::new(2.0, 1.0, 0.0),
        Point3::new(0.0, 1.0, 0.0),
    ];
    let block = BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), 0.5, tol)?;
    // Bore cutter overshoots both faces (transversal intersections).
    let bore =
        BRep::cylinder(0.2, 0.7, tol)?.transform(&Transform3::translation(0.5, 0.5, -0.1))?;
    // Countersink cutter: cone(0.35, 0.35) flipped apex-down, apex at
    // (0.5, 0.5, 0.2); its 45° flank crosses the bore wall at z = 0.4
    // and the top face plane (z = 0.5) at radius 0.3.
    let sink = BRep::cone(0.35, 0.35, tol)?
        .transform(&Transform3::rotation_x(PI))?
        .transform(&Transform3::translation(0.5, 0.5, 0.55))?;
    let result = part(block, tol, "bracket block")?
        .subtract(&part(bore, tol, "bore cutter")?)?
        .subtract(&part(sink, tol, "countersink cutter")?)?;
    // Removed material: bore π·0.2²·0.5, countersink ring beyond the
    // bore between z = 0.4 and 0.5: π∫((z−0.2)² − 0.2²)dz.
    let bore_v = PI * 0.04 * 0.5;
    let sink_v = PI * (((0.3_f64.powi(3) - 0.2_f64.powi(3)) / 3.0) - 0.04 * 0.1);
    Ok(CorpusFixture {
        name: "bracket_bore_countersink",
        brep: result.brep().clone(),
        watertight: true,
        // 6 block faces + bore wall + countersink flank.
        expected_faces: 8,
        // Block 12 + bore rim at bottom + cone/bore junction circle +
        // countersink rim at top + 2 seams (bore wall, sink flank).
        expected_edges: 17,
        // Block 8 + one seam-anchor vertex per circular edge ring.
        expected_vertices: 11,
        analytic_volume: Some(1.0 - bore_v - sink_v),
        volume_band: 0.01,
        steiner: None, // not convex (hole)
        expected_census: vec![("cone", 1), ("cylinder", 1), ("plane", 6)],
    })
}

/// L-section prism (concave edge along the inner corner) — the gap /
/// mixed-corner fixture for outward offsets.
pub fn l_bracket(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let profile = [
        Point3::new(0.0, 0.0, 0.0),
        Point3::new(1.5, 0.0, 0.0),
        Point3::new(1.5, 0.4, 0.0),
        Point3::new(0.4, 0.4, 0.0),
        Point3::new(0.4, 1.2, 0.0),
        Point3::new(0.0, 1.2, 0.0),
    ];
    let brep = BRep::prism(&profile, Vec3::new(0.0, 0.0, 1.0), 0.6, tol)?;
    Ok(CorpusFixture {
        name: "l_bracket",
        brep,
        watertight: true,
        // 6 side walls + 2 caps; V−E+F = 12−18+8 = 2.
        expected_faces: 8,
        expected_edges: 18,
        expected_vertices: 12,
        analytic_volume: Some((1.5 * 0.4 + 0.4 * 0.8) * 0.6),
        volume_band: 1e-6,
        steiner: None, // not convex (inner corner)
        expected_census: vec![("plane", 8)],
    })
}

/// Thin vertical rib (0.06 wide) standing on a base plate — the
/// face-collapse fixture: inward offsets beyond 0.03 must collapse the
/// rib, outward offsets must not weld its flanks.
pub fn thin_rib(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let plate = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(2.0, 0.0, 0.0),
            Point3::new(2.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        0.2,
        tol,
    )?;
    let rib = BRep::prism(
        &[
            Point3::new(0.2, 0.47, 0.2),
            Point3::new(1.8, 0.47, 0.2),
            Point3::new(1.8, 0.53, 0.2),
            Point3::new(0.2, 0.53, 0.2),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        0.5,
        tol,
    )?;
    let result = part(plate, tol, "rib base plate")?.union(&part(rib, tol, "rib")?)?;
    Ok(CorpusFixture {
        name: "thin_rib",
        brep: result.brep().clone(),
        watertight: true,
        // Plate 6 (top face gains the rib-footprint inner loop) + rib
        // walls 4 + rib top 1.
        expected_faces: 11,
        // Plate 12 + rib top ring 4 + footprint ring 4 + verticals 4.
        expected_edges: 24,
        expected_vertices: 16,
        analytic_volume: Some(2.0 * 1.0 * 0.2 + 1.6 * 0.06 * 0.5),
        volume_band: 1e-6,
        steiner: None, // not convex
        expected_census: vec![("plane", 11)],
    })
}

/// Plate whose top face is a free-form quadratic NURBS bulge (only the
/// center control point raised) — the certified-approximation fixture.
/// The patch boundary rows stay straight, so the plate's original
/// edges lie on the grafted surface exactly.
pub fn wavy_plate(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let (t, h) = (0.3, 0.25);
    let mut brep = BRep::prism(
        &[
            Point3::new(0.0, 0.0, 0.0),
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(1.0, 1.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        t,
        tol,
    )?;
    let handles: Vec<_> = brep.faces.keys().collect();
    let top = handles
        .into_iter()
        .find(|&fh| {
            matches!(
                brep.faces[fh].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Plane { .. })
            ) && {
                let ps = brep.face_outer_loop_positions(fh);
                !ps.is_empty() && ps.iter().all(|p| (p.z - t).abs() < 1e-9)
            }
        })
        .ok_or_else(|| CorpusError::Core(CoreError::PreconditionFailed("no top face".into())))?;
    let z = |raised: bool| if raised { t + h } else { t };
    let mut control_points = Vec::with_capacity(9);
    for yi in 0..3 {
        for xi in 0..3 {
            control_points.push(Point3::new(
                xi as f64 * 0.5,
                yi as f64 * 0.5,
                z(xi == 1 && yi == 1),
            ));
        }
    }
    brep.faces[top].surface = Some(SurfaceOnFace {
        geometry: SurfaceGeometry::Nurbs(NurbsSurface {
            degree_u: 2,
            degree_v: 2,
            knots_u: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            knots_v: vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0],
            control_points,
            num_u: 3,
            num_v: 3,
            weights: vec![1.0; 9],
        }),
        orientation: true,
    });
    brep.recompute_face_pcurves(top)?;
    // Bulge volume: h·∫B₁²(u)du·∫B₁²(v)dv = h·(1/3)² over the unit
    // square, on top of the t-thick plate.
    Ok(CorpusFixture {
        name: "wavy_plate",
        brep,
        watertight: true,
        expected_faces: 6,
        expected_edges: 12,
        expected_vertices: 8,
        analytic_volume: Some(t + h / 9.0),
        volume_band: 0.01,
        steiner: None, // free-form top
        expected_census: vec![("nurbs", 1), ("plane", 5)],
    })
}

/// Torus with a narrow center hole (major 0.6, minor 0.45, hole radius
/// 0.15): outward offsets ≥ 0.15 self-intersect globally at the hole —
/// the global-SI fixture. Tube-formula constants hold below the reach.
pub fn torus_lens(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let (major, minor) = (0.6, 0.45);
    let brep = BRep::torus(major, minor, tol)?;
    Ok(CorpusFixture {
        name: "torus_lens",
        brep,
        watertight: true,
        // Constructor chart layout: 4×4 grid of tube patches.
        expected_faces: 16,
        expected_edges: 32,
        expected_vertices: 16,
        analytic_volume: Some(2.0 * PI * PI * major * minor * minor),
        volume_band: 0.005,
        steiner: Some(SteinerConstants {
            volume: 2.0 * PI * PI * major * minor * minor,
            surface_area: 4.0 * PI * PI * major * minor,
            mean_curvature_integral: 2.0 * PI * PI * major,
            euler_characteristic: 0.0,
        }),
        expected_census: vec![("torus", 16)],
    })
}

/// Cylinder (r 0.5, h 1) with a circumferential rectangular groove
/// (width 0.1 at mid-height, floor at r 0.3) — the thin-wall global-SI
/// fixture. Outward offsets `d > 0.05` collide the groove walls: the
/// wall/floor offsets are unrealized and the two mouth-rim join tori
/// cross in an exact coaxial circle at `z = 0.5`,
/// `r = 0.5 + sqrt(d² − 0.05²)`. Rotational symmetry keeps every
/// self-intersection curve closed (no curve-endpoint termination on
/// corner faces), so the fixture isolates the trim/re-close mechanism.
pub fn grooved_cylinder(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    // Stacked coaxial unions (coincident-cap booleans): band + core +
    // band. The subtractive construction (cylinder minus tube) trips
    // PK-341 — the boolean drops the island-split annular groove walls
    // and serves a non-closed body — so the fixture builds additively.
    let bottom = BRep::cylinder(0.5, 0.45, tol)?;
    let core =
        BRep::cylinder(0.3, 0.1, tol)?.transform(&Transform3::translation(0.0, 0.0, 0.45))?;
    let top =
        BRep::cylinder(0.5, 0.45, tol)?.transform(&Transform3::translation(0.0, 0.0, 0.55))?;
    let result = part(bottom, tol, "grooved cylinder lower band")?
        .union(&part(core, tol, "groove core")?)?
        .union(&part(top, tol, "grooved cylinder upper band")?)?;
    let brep = result.brep().clone();
    Ok(CorpusFixture {
        name: "grooved_cylinder",
        brep,
        watertight: true,
        // Outer wall split into two bands + groove floor (3 cylinders),
        // two caps + two groove walls (4 planes).
        expected_faces: 7,
        // Six rim circles (two cap rims, two groove mouth rims, two
        // floor rims) + three wall seams (two bands + floor).
        expected_edges: 9,
        // One anchor vertex per rim circle; seams reuse the anchors.
        expected_vertices: 6,
        // Cylinder minus the groove ring: π(0.25·1 − (0.25 − 0.09)·0.1).
        analytic_volume: Some(PI * (0.25 - 0.016)),
        // Inscribed-chord tessellation deficit on the cylinder walls
        // measures 0.64% — a mesh-density calibration, not geometry.
        volume_band: 0.01,
        steiner: None, // not convex (groove)
        expected_census: vec![("cylinder", 3), ("plane", 4)],
    })
}

/// Open five-face box (one face removed, rim exposed) — the
/// thicken-rims fixture. Not watertight by construction.
pub fn open_sheet(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let brep = BRep::open_box(tol)?;
    Ok(CorpusFixture {
        name: "open_sheet",
        brep,
        watertight: false,
        expected_faces: 5,
        expected_edges: 12,
        expected_vertices: 8,
        analytic_volume: None,
        volume_band: 0.0,
        steiner: None,
        expected_census: vec![("plane", 5)],
    })
}

/// Cylinder rod (r 0.3, h 1.2) with a ball end (r 0.4 centered at
/// z 1.1): one transversal junction circle at z = 1.1 − √0.07, the
/// rod's top cap consumed inside the ball. Coaxial mixed analytic
/// junction under union.
pub fn sphere_capped_rod(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let rod = BRep::cylinder(0.3, 1.2, tol)?;
    let ball = BRep::sphere(0.4, tol)?.transform(&Transform3::translation(0.0, 0.0, 1.1))?;
    let result = part(rod, tol, "rod")?.union(&part(ball, tol, "ball end")?)?;
    // Coaxial union: cylinder below the junction plane z₀ = 1.1 − √0.07
    // plus the spherical cap above it (cap height H = 0.4 + √0.07).
    let a = 0.07_f64.sqrt();
    let z0 = 1.1 - a;
    let cap_h = 0.4 + a;
    let volume = PI * 0.09 * z0 + PI * cap_h * cap_h * (3.0 * 0.4 - cap_h) / 3.0;
    Ok(CorpusFixture {
        name: "sphere_capped_rod",
        brep: result.brep().clone(),
        watertight: true,
        // Bottom cap + cylinder wall + sphere cap; edges: bottom rim,
        // wall seam, junction circle, and the cap seam from the junction
        // to the pole. The pins were 6/5 until PK-921: the cap seam was
        // split where the rod's consumed top plane (r 0.387, z 1.2) and
        // its unbounded cylinder (z 1.365, above the rod) met it, hits
        // on no face's trim.
        expected_faces: 3,
        expected_edges: 4,
        expected_vertices: 3,
        analytic_volume: Some(volume),
        // PK-339: the cap tessellates through the structured
        // revolution path at exactly its chord tolerance; the
        // inscribed-mesh volume deficit (measured 1.28%) is faceting,
        // not geometry — the band is a mesh-density calibration.
        volume_band: 0.02,
        steiner: None, // convex only piecewise; junction is G0
        expected_census: vec![("cylinder", 1), ("plane", 1), ("sphere", 1)],
    })
}

/// Hexagonal pyramid — the high-valence apex fixture (valence 6).
pub fn high_valence_pyramid(tol: Tolerance) -> Result<CorpusFixture, CorpusError> {
    let (n, r, h) = (6_usize, 1.0, 1.2);
    let brep = BRep::n_gonal_pyramid(n, r, h, tol)?;
    // Regular n-gon of circumradius r: area = (n/2)r²·sin(2π/n).
    let base_area = 0.5 * n as f64 * r * r * (2.0 * PI / n as f64).sin();
    Ok(CorpusFixture {
        name: "high_valence_pyramid",
        brep,
        watertight: true,
        // n lateral triangles + base; V−E+F = 7−12+7 = 2.
        expected_faces: 7,
        expected_edges: 12,
        expected_vertices: 7,
        analytic_volume: Some(base_area * h / 3.0),
        volume_band: 1e-6,
        steiner: None, // convex but constants add little over the volume pin
        expected_census: vec![("plane", 7)],
    })
}

/// Build the whole corpus in a stable order.
pub fn corpus(tol: Tolerance) -> Result<Vec<CorpusFixture>, CorpusError> {
    Ok(vec![
        filleted_cube(tol)?,
        bracket_bore_countersink(tol)?,
        l_bracket(tol)?,
        thin_rib(tol)?,
        wavy_plate(tol)?,
        torus_lens(tol)?,
        grooved_cylinder(tol)?,
        open_sheet(tol)?,
        sphere_capped_rod(tol)?,
        high_valence_pyramid(tol)?,
    ])
}
