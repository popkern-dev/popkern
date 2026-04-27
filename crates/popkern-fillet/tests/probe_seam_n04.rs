//! Diagnostic: dump per-side seam trajectories on the worst Cyl×Gregory
//! edge of pyramid_n04. We want to know whether:
//!   (a) cyl.surf(cyl.pcurve(t)) tracks the 3-D edge curve;
//!   (b) gregory.surf(gregory.pcurve(t)) tracks it;
//!   (c) where the 0.153 divergence comes from.

use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::{SurfaceEval, SurfaceGeometry};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::BRep;
use popkern_fillet::{FilletOptions, FilletRadius, FilletSpec, fillet_edges};

#[test]
#[ignore]
fn probe_seam_n04() {
    let pyr = BRep::n_gonal_pyramid(4, 1.0, 1.2, Tolerance::default()).unwrap();
    let apex = pyr
        .vertices
        .iter()
        .find(|(_, v)| v.position.x.abs() < 1e-9 && v.position.y.abs() < 1e-9 && v.position.z > 0.0)
        .map(|(h, _)| h)
        .expect("apex");
    let edges: Vec<_> = pyr
        .half_edges
        .values()
        .filter(|he| he.origin == apex)
        .map(|he| he.edge)
        .collect::<std::collections::BTreeSet<_>>()
        .into_iter()
        .collect();
    let radii = [0.07, 0.09, 0.11, 0.13];
    let specs: Vec<FilletSpec> = edges
        .iter()
        .zip(radii.iter())
        .map(|(&e, &r)| FilletSpec {
            edge: e,
            radius: FilletRadius::Constant(r),
        })
        .collect();
    let opts = FilletOptions::default();
    let result = fillet_edges(&pyr, &specs, &opts, Tolerance::default(), None).unwrap();
    let brep = &result.brep;

    // Find the worst Cyl×Gregory tess-nm edge.
    let reports = brep.find_tessellation_non_manifold_edges(32, 1e-9);
    let mut worst: Option<&popkern_core::validation::NonManifoldEdgeReport> = None;
    for r in &reports {
        let kind_a = brep.faces[r.face_a]
            .surface
            .as_ref()
            .map(|s| variant(&s.geometry))
            .unwrap_or("?");
        let kind_b = brep.faces[r.face_b]
            .surface
            .as_ref()
            .map(|s| variant(&s.geometry))
            .unwrap_or("?");
        if (kind_a == "Cyl" && kind_b == "Gregory") || (kind_a == "Gregory" && kind_b == "Cyl") {
            if worst.is_none() || r.max_divergence > worst.unwrap().max_divergence {
                worst = Some(r);
            }
        }
    }
    let worst = worst.expect("at least one Cyl×Gregory edge");
    eprintln!(
        "WORST Cyl×Gregory edge: {:?}, faces {:?}/{:?}, div={:.4}, at {:?}",
        worst.edge, worst.face_a, worst.face_b, worst.max_divergence, worst.location
    );

    // Identify which face is cylinder and which is Gregory.
    let surf_a = &brep.faces[worst.face_a].surface.as_ref().unwrap().geometry;
    let surf_b = &brep.faces[worst.face_b].surface.as_ref().unwrap().geometry;
    let (cyl_face, cyl_surf, greg_face, greg_surf) =
        if matches!(surf_a, SurfaceGeometry::Cylinder { .. }) {
            (worst.face_a, surf_a, worst.face_b, surf_b)
        } else {
            (worst.face_b, surf_b, worst.face_a, surf_a)
        };
    eprintln!("  cyl_face={:?}, greg_face={:?}", cyl_face, greg_face);

    // Get the two half-edges on this edge and identify which side each is.
    let e = &brep.edges[worst.edge];
    let he_a = e.half_edge;
    let he_b = brep.half_edges[he_a].twin;
    let face_of = |he| brep.loops[brep.half_edges[he].loop_].face;
    let (cyl_he, greg_he) = if face_of(he_a) == cyl_face {
        (he_a, he_b)
    } else {
        (he_b, he_a)
    };
    let cyl_pc = brep.half_edges[cyl_he].pcurve.as_ref().unwrap();
    let greg_pc = brep.half_edges[greg_he].pcurve.as_ref().unwrap();
    let coe = e.curve.as_ref().unwrap();
    eprintln!("  edge curve t in [{:.4}, {:.4}]", coe.t_start, coe.t_end);
    let (cpa, cpa2) = cyl_pc.parameter_range();
    let (gpa, gpa2) = greg_pc.parameter_range();
    eprintln!(
        "  cyl pcurve param [{:.4}, {:.4}], greg pcurve param [{:.4}, {:.4}]",
        cpa, cpa2, gpa, gpa2
    );

    // For each of 11 samples along the edge, dump:
    //   p_3d   = edge.point_at(t)
    //   p_cyl  = cyl.point_at(cyl_pc(t_cyl))
    //   p_greg = greg.point_at(greg_pc(t_greg))
    //   gap_cyl  = ||p_3d - p_cyl||
    //   gap_greg = ||p_3d - p_greg||
    //   gap_cg   = ||p_cyl - p_greg||  (this is what tess-nm measures)
    eprintln!(
        "\n  s        p_3d                         p_cyl   p_greg  gap_cyl   gap_greg  gap_cg"
    );

    // Direction handling — same logic as `find_tessellation_non_manifold_edges`.
    let p3_a = coe.geometry.point_at(coe.t_start);
    let p3_b = coe.geometry.point_at(coe.t_end);
    let cyl_p_lo = cyl_surf.point_at(cyl_pc.point_at(cpa).x, cyl_pc.point_at(cpa).y);
    let cyl_p_hi = cyl_surf.point_at(cyl_pc.point_at(cpa2).x, cyl_pc.point_at(cpa2).y);
    let cyl_forward = p3_a.distance_to(&cyl_p_lo) + p3_b.distance_to(&cyl_p_hi)
        <= p3_a.distance_to(&cyl_p_hi) + p3_b.distance_to(&cyl_p_lo);
    let greg_p_lo = greg_surf.point_at(greg_pc.point_at(gpa).x, greg_pc.point_at(gpa).y);
    let greg_p_hi = greg_surf.point_at(greg_pc.point_at(gpa2).x, greg_pc.point_at(gpa2).y);
    let greg_forward = p3_a.distance_to(&greg_p_lo) + p3_b.distance_to(&greg_p_hi)
        <= p3_a.distance_to(&greg_p_hi) + p3_b.distance_to(&greg_p_lo);

    for i in 0..=10 {
        let s = i as f64 / 10.0;
        let t_curve = coe.t_start + s * (coe.t_end - coe.t_start);
        let p3 = coe.geometry.point_at(t_curve);
        let s_cyl = if cyl_forward { s } else { 1.0 - s };
        let s_greg = if greg_forward { s } else { 1.0 - s };
        let t_c = cpa + s_cyl * (cpa2 - cpa);
        let t_g = gpa + s_greg * (gpa2 - gpa);
        let uv_c = cyl_pc.point_at(t_c);
        let uv_g = greg_pc.point_at(t_g);
        let p_c = cyl_surf.point_at(uv_c.x, uv_c.y);
        let p_g = greg_surf.point_at(uv_g.x, uv_g.y);
        eprintln!(
            "  {:.2}   ({:+.4},{:+.4},{:+.4})  cylUV=({:+.3},{:+.3})  gregUV=({:+.3},{:+.3})  {:.3e}  {:.3e}  {:.3e}",
            s,
            p3.x,
            p3.y,
            p3.z,
            uv_c.x,
            uv_c.y,
            uv_g.x,
            uv_g.y,
            p3.distance_to(&p_c),
            p3.distance_to(&p_g),
            p_c.distance_to(&p_g),
        );
    }

    // Also check Gregory's boundary geometry directly: at the side
    // it should be matching, evaluate Gregory(u, 0) or similar at
    // many samples and compare to the 3D edge curve. The sub-patch
    // boundary should geometrically equal the edge curve.
    let ((u_lo, u_hi), (v_lo, v_hi)) = greg_surf.parameter_range();
    eprintln!(
        "\n  greg UV box: u[{:.3},{:.3}] v[{:.3},{:.3}]",
        u_lo, u_hi, v_lo, v_hi
    );
    eprintln!("  Sample greg boundaries vs edge endpoints:");
    let endpoints = [
        (
            "u=0",
            greg_surf.point_at(u_lo, v_lo),
            greg_surf.point_at(u_lo, v_hi),
        ),
        (
            "u=1",
            greg_surf.point_at(u_hi, v_lo),
            greg_surf.point_at(u_hi, v_hi),
        ),
        (
            "v=0",
            greg_surf.point_at(u_lo, v_lo),
            greg_surf.point_at(u_hi, v_lo),
        ),
        (
            "v=1",
            greg_surf.point_at(u_lo, v_hi),
            greg_surf.point_at(u_hi, v_hi),
        ),
    ];
    for (label, p_a, p_b) in endpoints.iter() {
        let d_a_to_start = p_a.distance_to(&p3_a);
        let d_a_to_end = p_a.distance_to(&p3_b);
        let d_b_to_start = p_b.distance_to(&p3_a);
        let d_b_to_end = p_b.distance_to(&p3_b);
        eprintln!(
            "    {} corners: d_to_edge_start={:.4}/{:.4}, d_to_edge_end={:.4}/{:.4}",
            label, d_a_to_start, d_b_to_start, d_a_to_end, d_b_to_end
        );
    }

    // Disambiguation #4 — what IS the BRep edge curve? Dump the
    // CurveGeometry variant + its parameter range vs the edge's
    // (t_start, t_end). If t range is a small slice (e.g. 0.92..0.94)
    // of a larger [0,1] curve, the edge is a TRIM of a pre-existing
    // longer curve — likely the cylinder fillet's apex-side trim curve.
    eprintln!("\n  === BRep edge curve structure ===");
    use popkern_core::geometry::curve::CurveGeometry;
    eprintln!(
        "  edge t range: [{:.4}, {:.4}] (span = {:.4})",
        coe.t_start,
        coe.t_end,
        coe.t_end - coe.t_start
    );
    match &coe.geometry {
        CurveGeometry::Nurbs(nc) => {
            let (full_lo, full_hi) = (nc.knots[0], *nc.knots.last().unwrap());
            eprintln!(
                "  full curve param range: [{:.4}, {:.4}]  → edge slice covers {:.1}% of full curve",
                full_lo,
                full_hi,
                100.0 * (coe.t_end - coe.t_start) / (full_hi - full_lo).max(1e-15)
            );
            eprintln!(
                "  NURBS: degree={}, n_cps={}, n_knots={}",
                nc.degree,
                nc.control_points.len(),
                nc.knots.len()
            );
            eprintln!("  endpoints of full curve:");
            let p_lo = nc.point_at(full_lo);
            let p_hi = nc.point_at(full_hi);
            eprintln!(
                "    full(t={:.4}) = ({:+.4},{:+.4},{:+.4})",
                full_lo, p_lo.x, p_lo.y, p_lo.z
            );
            eprintln!(
                "    full(t={:.4}) = ({:+.4},{:+.4},{:+.4})",
                full_hi, p_hi.x, p_hi.y, p_hi.z
            );
        }
        other => eprintln!(
            "  curve geometry variant: {:?}",
            std::mem::discriminant(other)
        ),
    }

    // Compare to the trim curves on the fillet (any fillet whose
    // FilletEnd would correspond to this edge).
    // We need access to the FilletSurfaceResult. Re-build it using
    // the same path the production code would.
    // (Skip — easier to just see the t range slice in the BRep.)

    // Disambiguation #3 — verbatim corner CPs vs BRep edge endpoint vertices.
    //
    // What we want to know: did the Gregory cap construction set
    // its (u=0, v=0) and (u=0, v=1) corner CPs at the SAME 3-D
    // positions as the BRep edge endpoint vertices?
    //
    // If YES — corner CPs match BRep vertices → boundary geometry
    //          IS the edge curve, and the bug is elsewhere
    //          (orientation? pcurve assumes wrong UV side?)
    // If NO  — corner CPs differ from BRep vertices → either the
    //          construction's hex_vertices drift from BRep authoring,
    //          or a later pass overwrote the corner CPs.
    eprintln!("\n  === Verbatim corner CPs vs BRep edge endpoint vertices ===");
    let v_origin = brep.vertices[brep.half_edges[he_a].origin].position;
    let v_dest_he = brep.half_edges[he_a].twin;
    let v_dest = brep.vertices[brep.half_edges[v_dest_he].origin].position;
    eprintln!(
        "  BRep edge endpoint vertices: ({:+.4},{:+.4},{:+.4}) → ({:+.4},{:+.4},{:+.4})",
        v_origin.x, v_origin.y, v_origin.z, v_dest.x, v_dest.y, v_dest.z
    );
    eprintln!(
        "  Edge curve t={:.4} ({:+.4},{:+.4},{:+.4})",
        coe.t_start, p3_a.x, p3_a.y, p3_a.z,
    );
    eprintln!(
        "  Edge curve t={:.4} ({:+.4},{:+.4},{:+.4})",
        coe.t_end, p3_b.x, p3_b.y, p3_b.z,
    );

    // Read all 4 corner positions of the Gregory directly from
    // its boundary_cps grid. Variant-aware.
    let greg_patch = match greg_surf {
        SurfaceGeometry::Gregory(gp) => gp,
        _ => panic!("greg_surf is not Gregory"),
    };
    use popkern_core::geometry::gregory::GregoryPatch;
    let (c00, c05, c55, c50, variant) = match greg_patch {
        GregoryPatch::Rectangular(r) => (
            r.boundary_cps[0][0],
            r.boundary_cps[0][3],
            r.boundary_cps[3][3],
            r.boundary_cps[3][0],
            "Rectangular(cubic)",
        ),
        GregoryPatch::RectangularDeg5(r) => (
            r.boundary_cps[0][0],
            r.boundary_cps[0][5],
            r.boundary_cps[5][5],
            r.boundary_cps[5][0],
            "RectangularDeg5",
        ),
        GregoryPatch::RectangularBiquinticFull(r) => (
            r.grid[0][0],
            r.grid[0][5],
            r.grid[5][5],
            r.grid[5][0],
            "RectangularBiquinticFull",
        ),
        GregoryPatch::RectangularBiNonic(r) => (
            r.grid[0][0],
            r.grid[0][9],
            r.grid[9][9],
            r.grid[9][0],
            "RectangularBiNonic",
        ),
        other => {
            eprintln!(
                "  unhandled Gregory variant: {:?}",
                std::mem::discriminant(other)
            );
            return;
        }
    };
    eprintln!("  Gregory variant: {variant}");
    eprintln!(
        "  Corner CP (0,0) = ({:+.4},{:+.4},{:+.4})",
        c00.x, c00.y, c00.z
    );
    eprintln!(
        "  Corner CP (0,1) = ({:+.4},{:+.4},{:+.4})",
        c05.x, c05.y, c05.z
    );
    eprintln!(
        "  Corner CP (1,1) = ({:+.4},{:+.4},{:+.4})",
        c55.x, c55.y, c55.z
    );
    eprintln!(
        "  Corner CP (1,0) = ({:+.4},{:+.4},{:+.4})",
        c50.x, c50.y, c50.z
    );
    eprintln!();
    eprintln!("  Distances from each Gregory corner CP to each BRep edge vertex:");
    let corners = [
        ("(0,0)", c00),
        ("(0,1)", c05),
        ("(1,1)", c55),
        ("(1,0)", c50),
    ];
    for (label, cp) in &corners {
        let d_to_origin = cp.distance_to(&v_origin);
        let d_to_dest = cp.distance_to(&v_dest);
        eprintln!(
            "    cp{}: d_to_v_origin={:.4}, d_to_v_dest={:.4}",
            label, d_to_origin, d_to_dest
        );
    }

    // Disambiguation #2: collect ALL Gregory cap faces in the BRep.
    // For each of the 4 cap faces, find the minimum distance from
    // the worst-edge midpoint to that Gregory's geometric boundary.
    // If a *different* Gregory than the one assigned to this face
    // is closer, the assignment-to-face mapping is the bug.
    eprintln!("\n  All Gregory cap faces vs the worst edge midpoint:");
    let mid_t = (coe.t_start + coe.t_end) * 0.5;
    let p_mid = coe.geometry.point_at(mid_t);
    eprintln!(
        "    edge midpoint: ({:+.4},{:+.4},{:+.4})",
        p_mid.x, p_mid.y, p_mid.z
    );
    let n_bdry_samples = 200;
    let cap_faces: Vec<_> = brep
        .faces
        .iter()
        .filter_map(|(fh, f)| {
            f.surface.as_ref().and_then(|sof| {
                if matches!(sof.geometry, SurfaceGeometry::Gregory(_)) {
                    Some((fh, &sof.geometry))
                } else {
                    None
                }
            })
        })
        .collect();
    for (fh, gs) in cap_faces.iter() {
        let ((u_lo, u_hi), (v_lo, v_hi)) = gs.parameter_range();
        let mut min_d = f64::INFINITY;
        let mut min_uv = (0.0, 0.0);
        for j in 0..=n_bdry_samples {
            let p = j as f64 / n_bdry_samples as f64;
            for (uu, vv) in &[
                (u_lo, v_lo + p * (v_hi - v_lo)),
                (u_hi, v_lo + p * (v_hi - v_lo)),
                (u_lo + p * (u_hi - u_lo), v_lo),
                (u_lo + p * (u_hi - u_lo), v_hi),
            ] {
                let q = gs.point_at(*uu, *vv);
                let d = p_mid.distance_to(&q);
                if d < min_d {
                    min_d = d;
                    min_uv = (*uu, *vv);
                }
            }
        }
        let marker = if *fh == greg_face {
            " ← assigned to this edge"
        } else {
            ""
        };
        eprintln!(
            "    Gregory face {:?}: min boundary distance to edge midpoint = {:.4} at uv=({:.3},{:.3}){}",
            fh, min_d, min_uv.0, min_uv.1, marker
        );
    }
}

fn variant(g: &SurfaceGeometry) -> &'static str {
    match g {
        SurfaceGeometry::Plane { .. } => "Plane",
        SurfaceGeometry::Cylinder { .. } => "Cyl",
        SurfaceGeometry::Sphere { .. } => "Sph",
        SurfaceGeometry::Torus { .. } => "Tor",
        SurfaceGeometry::Cone { .. } => "Cone",
        SurfaceGeometry::Nurbs(_) => "Nurbs",
        SurfaceGeometry::Gregory(_) => "Gregory",
    }
}
