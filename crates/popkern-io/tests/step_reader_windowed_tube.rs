//! PK-652: a seamless cylinder wall with a hole meshes through the
//! annular-band slit, and the hole stays a hole. The slit's band polygon
//! is the face's outer boundary; appended after the hole rings, it left
//! the holes a period off its chart. Their region then meshed as
//! material and the interior lattice was sized from a hole's box
//! (corpus bracket1, a cross-drilled shaft: facets 2.11 off at a 0.005
//! chord, still 0.0907 after deviation refinement).
//!
//! The fixture is a thick tube (outer radius R, bore r, height h, axis
//! z) with a rectangular window cut through one side (|x| < a, y > 0,
//! |z - zc| < b): the window walls are planes, so every edge is exact
//! (rulings and arcs on the cylinders, radial lines between them), and
//! both cylinder faces are bounded only by their rim circles and the
//! window: wrapped rings, no seam edge.

use popkern_io::{StepImportOptions, import_step};

/// STEP entity list with sequential ids.
struct Step {
    lines: Vec<String>,
}

impl Step {
    fn add(&mut self, entity: String) -> usize {
        self.lines.push(entity);
        self.lines.len() + 9
    }
    fn point(&mut self, p: [f64; 3]) -> usize {
        self.add(format!(
            "CARTESIAN_POINT('',({:.17e},{:.17e},{:.17e}))",
            p[0], p[1], p[2]
        ))
    }
    fn dir(&mut self, d: [f64; 3]) -> usize {
        self.add(format!(
            "DIRECTION('',({:.17e},{:.17e},{:.17e}))",
            d[0], d[1], d[2]
        ))
    }
    fn placement(&mut self, origin: [f64; 3], axis: [f64; 3], reference: [f64; 3]) -> usize {
        let (o, a, r) = (self.point(origin), self.dir(axis), self.dir(reference));
        self.add(format!("AXIS2_PLACEMENT_3D('',#{o},#{a},#{r})"))
    }
    fn vertex(&mut self, p: [f64; 3]) -> usize {
        let pt = self.point(p);
        self.add(format!("VERTEX_POINT('',#{pt})"))
    }
    fn circle_edge(&mut self, z: f64, radius: f64, from: usize, to: usize) -> usize {
        let pl = self.placement([0.0, 0.0, z], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
        let c = self.add(format!("CIRCLE('',#{pl},{radius:.17e})"));
        self.add(format!("EDGE_CURVE('',#{from},#{to},#{c},.T.)"))
    }
    fn line_edge(&mut self, from: (usize, [f64; 3]), to: (usize, [f64; 3])) -> usize {
        let d = [
            to.1[0] - from.1[0],
            to.1[1] - from.1[1],
            to.1[2] - from.1[2],
        ];
        let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
        let p = self.point(from.1);
        let dir = self.dir([d[0] / len, d[1] / len, d[2] / len]);
        let v = self.add(format!("VECTOR('',#{dir},{len:.17e})"));
        let l = self.add(format!("LINE('',#{p},#{v})"));
        self.add(format!("EDGE_CURVE('',#{},#{},#{l},.T.)", from.0, to.0))
    }
    /// A face on `surface` whose bounds are loops of (edge, forward).
    fn face(&mut self, surface: usize, sense: bool, loops: &[&[(usize, bool)]]) -> usize {
        let mut bounds = Vec::new();
        for lp in loops {
            let oes: Vec<String> = lp
                .iter()
                .map(|&(e, fwd)| {
                    let id = self.add(format!(
                        "ORIENTED_EDGE('',*,*,#{e},{})",
                        if fwd { ".T." } else { ".F." }
                    ));
                    format!("#{id}")
                })
                .collect();
            let l = self.add(format!("EDGE_LOOP('',({}))", oes.join(",")));
            bounds.push(format!("#{}", self.add(format!("FACE_BOUND('',#{l},.T.)"))));
        }
        let s = if sense { ".T." } else { ".F." };
        self.add(format!(
            "ADVANCED_FACE('',({}),#{surface},{s})",
            bounds.join(",")
        ))
    }
}

/// The windowed tube as STEP text.
fn windowed_tube((big, small, h): (f64, f64, f64), (a, zc, b): (f64, f64, f64)) -> String {
    let mut s = Step { lines: Vec::new() };
    let (yo, yi) = ((big * big - a * a).sqrt(), (small * small - a * a).sqrt());
    let (z1, z2) = (zc - b, zc + b);
    // Rim vertices and circles.
    let o0 = s.vertex([big, 0.0, 0.0]);
    let o1 = s.vertex([big, 0.0, h]);
    let i0 = s.vertex([small, 0.0, 0.0]);
    let i1 = s.vertex([small, 0.0, h]);
    let co0 = s.circle_edge(0.0, big, o0, o0);
    let co1 = s.circle_edge(h, big, o1, o1);
    let ci0 = s.circle_edge(0.0, small, i0, i0);
    let ci1 = s.circle_edge(h, small, i1, i1);
    // Window corners: outer/inner, x = +a / -a, z = z1 / z2.
    let corner = |s: &mut Step, x: f64, y: f64, z: f64| ([x, y, z], s.vertex([x, y, z]));
    let (wo_p1, wo_p1v) = corner(&mut s, a, yo, z1);
    let (wo_p2, wo_p2v) = corner(&mut s, a, yo, z2);
    let (wo_m1, wo_m1v) = corner(&mut s, -a, yo, z1);
    let (wo_m2, wo_m2v) = corner(&mut s, -a, yo, z2);
    let (wi_p1, wi_p1v) = corner(&mut s, a, yi, z1);
    let (wi_p2, wi_p2v) = corner(&mut s, a, yi, z2);
    let (wi_m1, wi_m1v) = corner(&mut s, -a, yi, z1);
    let (wi_m2, wi_m2v) = corner(&mut s, -a, yi, z2);
    // Rulings (upward), arcs (counter-clockwise about z: x = +a to -a),
    // radial lines (inner to outer).
    let lo_p = s.line_edge((wo_p1v, wo_p1), (wo_p2v, wo_p2));
    let lo_m = s.line_edge((wo_m1v, wo_m1), (wo_m2v, wo_m2));
    let li_p = s.line_edge((wi_p1v, wi_p1), (wi_p2v, wi_p2));
    let li_m = s.line_edge((wi_m1v, wi_m1), (wi_m2v, wi_m2));
    let ao_1 = s.circle_edge(z1, big, wo_p1v, wo_m1v);
    let ao_2 = s.circle_edge(z2, big, wo_p2v, wo_m2v);
    let ai_1 = s.circle_edge(z1, small, wi_p1v, wi_m1v);
    let ai_2 = s.circle_edge(z2, small, wi_p2v, wi_m2v);
    let rd_p1 = s.line_edge((wi_p1v, wi_p1), (wo_p1v, wo_p1));
    let rd_p2 = s.line_edge((wi_p2v, wi_p2), (wo_p2v, wo_p2));
    let rd_m1 = s.line_edge((wi_m1v, wi_m1), (wo_m1v, wo_m1));
    let rd_m2 = s.line_edge((wi_m2v, wi_m2), (wo_m2v, wo_m2));
    // Surfaces.
    let axis = s.placement([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let outer = s.add(format!("CYLINDRICAL_SURFACE('',#{axis},{big:.17e})"));
    let bore = s.add(format!("CYLINDRICAL_SURFACE('',#{axis},{small:.17e})"));
    let bottom_pl = s.placement([0.0; 3], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let bottom = s.add(format!("PLANE('',#{bottom_pl})"));
    let top_pl = s.placement([0.0, 0.0, h], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let top = s.add(format!("PLANE('',#{top_pl})"));
    let wall = |s: &mut Step, origin: [f64; 3], normal: [f64; 3], reference: [f64; 3]| {
        let pl = s.placement(origin, normal, reference);
        s.add(format!("PLANE('',#{pl})"))
    };
    let wall_xp = wall(&mut s, [a, 0.0, 0.0], [-1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let wall_xm = wall(&mut s, [-a, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]);
    let wall_z2 = wall(&mut s, [0.0, 0.0, z2], [0.0, 0.0, -1.0], [1.0, 0.0, 0.0]);
    let wall_z1 = wall(&mut s, [0.0, 0.0, z1], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    // Faces, every loop counter-clockwise about the outward normal.
    let faces = [
        s.face(
            outer,
            true,
            &[
                &[(co0, true)],
                &[(co1, false)],
                &[(lo_p, true), (ao_2, true), (lo_m, false), (ao_1, false)],
            ],
        ),
        s.face(
            bore,
            false,
            &[
                &[(ci0, false)],
                &[(ci1, true)],
                &[(ai_1, true), (li_m, true), (ai_2, false), (li_p, false)],
            ],
        ),
        s.face(bottom, false, &[&[(co0, false)], &[(ci0, true)]]),
        s.face(top, true, &[&[(co1, true)], &[(ci1, false)]]),
        s.face(
            wall_xp,
            true,
            &[&[(li_p, true), (rd_p2, true), (lo_p, false), (rd_p1, false)]],
        ),
        s.face(
            wall_xm,
            true,
            &[&[(rd_m1, true), (lo_m, true), (rd_m2, false), (li_m, false)]],
        ),
        s.face(
            wall_z2,
            true,
            &[&[(ai_2, true), (rd_m2, true), (ao_2, false), (rd_p2, false)]],
        ),
        s.face(
            wall_z1,
            true,
            &[&[(rd_p1, true), (ao_1, true), (rd_m1, false), (ai_1, false)]],
        ),
    ];
    let refs: Vec<String> = faces.iter().map(|f| format!("#{f}")).collect();
    let shell = s.add(format!("CLOSED_SHELL('',({}))", refs.join(",")));
    let solid = s.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
    let mm = s.add("(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.))".into());
    let rad = s.add("(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.))".into());
    let ctx = s.add(format!(
        "(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{mm},#{rad})) REPRESENTATION_CONTEXT('',''))"
    ));
    s.add(format!(
        "ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{solid}),#{ctx})"
    ));
    let mut text = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('tube'),'2;1');\nFILE_NAME('tube.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n",
    );
    for (k, line) in s.lines.iter().enumerate() {
        text += &format!("#{}={line};\n", k + 10);
    }
    text + "ENDSEC;\nEND-ISO-10303-21;\n"
}

#[test]
fn a_windowed_tube_meshes_its_window_as_a_hole_within_twice_the_chord() {
    use popkern_core::geometry::surface::SurfaceEval;
    let (big, small, h) = (10.0, 6.0, 20.0);
    let (a, zc, b) = (2.0, 10.0, 3.0);
    let imported = import_step(
        windowed_tube((big, small, h), (a, zc, b)).as_bytes(),
        &StepImportOptions::default(),
    )
    .expect("import");
    let tube = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .next()
        .expect("solid");
    tube.validate().expect("valid");
    let options = popkern_mesh::TessellationOptions::default();
    let mesh = popkern_mesh::tessellate(&tube, &options).expect("mesh");
    assert_eq!(popkern_mesh::welded_open_edges(&mesh), 0, "watertight");
    let mut worst: f64 = 0.0;
    let mut volume = 0.0;
    for (ti, t) in mesh.triangles.iter().enumerate() {
        let surface = &tube.faces[mesh.face_map[ti]]
            .surface
            .as_ref()
            .expect("surface")
            .geometry;
        let [p, q, r] = t.map(|i| mesh.vertices[i].position);
        volume += (p.x * (q.y * r.z - q.z * r.y) - p.y * (q.x * r.z - q.z * r.x)
            + p.z * (q.x * r.y - q.y * r.x))
            / 6.0;
        let g = popkern_core::geometry::Point3::new(
            (p.x + q.x + r.x) / 3.0,
            (p.y + q.y + r.y) / 3.0,
            (p.z + q.z + r.z) / 3.0,
        );
        let (u, v) = surface.closest_point(&g).expect("projection");
        worst = worst.max(surface.point_at(u, v).distance_to(&g));
    }
    assert!(
        worst <= 2.0 * options.chord_height,
        "a facet lies {worst} off its surface"
    );
    // The window: |x| < a between the cylinders, over 2b of height.
    let segment = |rho: f64| {
        // Area under sqrt(rho^2 - x^2) for |x| < a.
        a * (rho * rho - a * a).sqrt() + rho * rho * (a / rho).asin()
    };
    let window = 2.0 * b * (segment(big) - segment(small));
    let exact = std::f64::consts::PI * (big * big - small * small) * h - window;
    assert!(
        (volume - exact).abs() < 2e-3 * exact,
        "signed volume {volume} vs {exact} (window {window})"
    );
}
