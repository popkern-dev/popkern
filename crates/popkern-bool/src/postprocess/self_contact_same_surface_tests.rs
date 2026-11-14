use super::*;
use popkern_core::geometry::curve::{CurveEval, CurveGeometry};
use popkern_core::geometry::{Point3, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::FaceHandle;

/// A small deterministic generator in [0, 1).
struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn vec(&mut self, s: f64) -> Vec3 {
        Vec3::new(
            s * (self.next() - 0.5),
            s * (self.next() - 0.5),
            s * (self.next() - 0.5),
        )
    }
}

/// A segment or an arc of up to a full turn, its parameter range.
fn curve(g: &mut Lcg, k: usize) -> (CurveGeometry, (f64, f64)) {
    let at = Point3::new(0.0, 0.0, 0.0) + g.vec(2.0);
    if k.is_multiple_of(3) {
        let d = g.vec(3.0);
        (
            CurveGeometry::Line {
                origin: at,
                direction: d,
            },
            (0.0, 1.0),
        )
    } else {
        let normal = g.vec(1.0).normalized().unwrap_or(Vec3::new(0.0, 0.0, 1.0));
        let span = 0.3 + 6.0 * g.next();
        (
            CurveGeometry::Circle {
                center: at,
                normal,
                radius: 0.2 + 1.5 * g.next(),
            },
            (0.4, 0.4 + span.min(std::f64::consts::TAU)),
        )
    }
}

#[test]
fn curves_cleared_lie_farther_than_the_pad() {
    let mut g = Lcg(31);
    let pad = 0.3;
    let (mut cleared, mut kept) = (0, 0);
    let n = 300;
    for k in 0..400 {
        let (a, ra) = curve(&mut g, k);
        let (b, rb) = curve(&mut g, k + 1);
        let (na, nb) = (a.to_nurbs_in_range(ra), b.to_nurbs_in_range(rb));
        if curves_apart(&na, &nb, pad, DEPTH) {
            cleared += 1;
            let sample = |c: &CurveGeometry, (t0, t1): (f64, f64)| -> Vec<Point3> {
                (0..=n)
                    .map(|i| c.point_at(t0 + (t1 - t0) * i as f64 / n as f64))
                    .collect()
            };
            let (pa, pb) = (sample(&a, ra), sample(&b, rb));
            for p in &pa {
                for q in &pb {
                    assert!(p.distance_to(q) > pad, "pair {k} cleared within the pad");
                }
            }
        } else {
            kept += 1;
        }
    }
    assert!(cleared > 40 && kept > 40, "cleared {cleared}, kept {kept}");
    // Two parallel unit segments 0.36 apart are cleared, 0.24 apart not.
    let seg = |y: f64| {
        CurveGeometry::Line {
            origin: Point3::new(0.0, y, 0.0),
            direction: Vec3::new(1.0, 0.0, 0.0),
        }
        .to_nurbs_in_range((0.0, 1.0))
    };
    assert!(curves_apart(&seg(0.0), &seg(0.36), pad, DEPTH));
    assert!(!curves_apart(&seg(0.0), &seg(0.24), pad, DEPTH));
}

#[test]
fn the_reach_is_known_only_where_it_is_positive() {
    let tol = Tolerance::default();
    let surface = |b: BRep| -> Vec<SurfaceGeometry> {
        b.faces
            .values()
            .filter_map(|f| f.surface.as_ref().map(|s| s.geometry.clone()))
            .collect()
    };
    let ring = surface(BRep::torus(5.0, 1.5, tol).expect("torus"));
    assert!(ring.iter().all(|s| reach(s) == Some(1.5)));
    let ball = surface(BRep::sphere(3.0, tol).expect("ball"));
    assert!(ball.iter().all(|s| reach(s) == Some(3.0)));
    let spindle = SurfaceGeometry::Torus {
        center: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        major_radius: 1.0,
        minor_radius: 1.5,
    };
    assert_eq!(reach(&spindle), None);
    let cone = SurfaceGeometry::Cone {
        apex: Point3::new(0.0, 0.0, 0.0),
        axis: Vec3::new(0.0, 0.0, 1.0),
        half_angle: 0.5,
    };
    assert_eq!(reach(&cone), None);
}

/// A torus (major 5, tube 1.5, about z) less a slab of thickness `gap`
/// across it from x = 0.2 (clear of its meridian edges at x = 0), with the pair of its torus faces sharing no
/// vertex whose edges come nearest (sampled), and that distance.
fn halved_torus(gap: f64) -> (BRep, FaceHandle, FaceHandle, f64) {
    let tol = Tolerance::default();
    let ring = BRep::torus(5.0, 1.5, tol).expect("torus");
    let (x0, x1) = (0.2, 0.2 + gap);
    let slab = BRep::prism(
        &[
            Point3::new(x0, -10.0, -5.0),
            Point3::new(x1, -10.0, -5.0),
            Point3::new(x1, 10.0, -5.0),
            Point3::new(x0, 10.0, -5.0),
        ],
        Vec3::new(0.0, 0.0, 1.0),
        10.0,
        tol,
    )
    .expect("slab");
    let body = crate::boolean(crate::BoolOperation::Subtraction, ring, slab, tol, None)
        .expect("torus less slab");
    let tori: Vec<FaceHandle> = body
        .faces
        .keys()
        .filter(|&f| {
            matches!(
                body.faces[f].surface.as_ref().map(|s| &s.geometry),
                Some(SurfaceGeometry::Torus { .. })
            )
        })
        .collect();
    let samples = |f: FaceHandle| -> Vec<Point3> {
        super::super::face_boundary(&body, f, &[])
            .edges
            .iter()
            .flat_map(|&e| {
                let c = body.edges[e].curve.clone().expect("curve");
                (0..=40).map(move |i| {
                    c.geometry
                        .point_at(c.t_start + (c.t_end - c.t_start) * i as f64 / 40.0)
                })
            })
            .collect()
    };
    let mut best = (tori[0], tori[0], f64::INFINITY);
    for (k, &f) in tori.iter().enumerate() {
        for &g in &tori[k + 1..] {
            let (vf, vg) = (
                super::super::face_boundary(&body, f, &[]).vertices,
                super::super::face_boundary(&body, g, &[]).vertices,
            );
            if vf.iter().any(|v| vg.contains(v)) {
                continue;
            }
            let (sf, sg) = (samples(f), samples(g));
            let d = sf
                .iter()
                .flat_map(|p| sg.iter().map(move |q| p.distance_to(q)))
                .fold(f64::INFINITY, f64::min);
            if d < best.2 {
                best = (f, g, d);
            }
        }
    }
    (body, best.0, best.1, best.2)
}

#[test]
fn pieces_of_one_torus_are_cleared_by_their_edges() {
    let (body, f, g, d) = halved_torus(0.01);
    assert!((d - 0.01).abs() < 1e-4, "the nearest rims lie {d} apart");
    let edges = |b: &BRep, h: FaceHandle| super::super::face_boundary(b, h, &[]).edges;
    let (ef, eg) = (edges(&body, f), edges(&body, g));
    let surface = body.faces[f]
        .surface
        .as_ref()
        .expect("surface")
        .geometry
        .clone();
    assert_eq!(
        Some(&surface),
        body.faces[g].surface.as_ref().map(|s| &s.geometry),
        "one surface"
    );
    // Cleared under the rims' distance, kept over it (the reach, 1.5,
    // over a hundred times the band in both).
    assert!(same_surface_clear(&body, &surface, (&ef, &eg), 1e-6));
    assert!(same_surface_clear(&body, &surface, (&ef, &eg), 0.009));
    assert!(!same_surface_clear(&body, &surface, (&ef, &eg), 0.011));
    assert!(super::super::self_contacts(&body, &Tolerance::default()).is_empty());
    // Rims 0.5 apart under a band of 0.02 the reach does not exceed a
    // hundred times: kept, though the rims are clear of it.
    let (wide, f, g, d) = halved_torus(0.5);
    assert!(d > 0.4, "{d}");
    let (ef, eg) = (edges(&wide, f), edges(&wide, g));
    assert!(same_surface_clear(&wide, &surface, (&ef, &eg), 0.014));
    assert!(!same_surface_clear(&wide, &surface, (&ef, &eg), 0.02));
}
