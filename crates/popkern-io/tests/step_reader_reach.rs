//! PK-631/PK-635 STEP reader conventions, each a variant of the unit
//! cube fixture as another writer spells it (measured on the CAx-IF /
//! STEP Tools corpus, docs/reports/step-assembly-corpus.md):
//!
//! - the product's representation reaches the solid through a plain
//!   `SHAPE_REPRESENTATION_RELATIONSHIP` (UG/NX/SolidWorks style);
//! - through a `MAPPED_ITEM` with a placement;
//! - a file with solids but no product structure;
//! - a shell listing a reversed `ORIENTED_FACE`;
//! - a face with an extra `VERTEX_LOOP` bound (a cone apex);
//! - outer loops marked only `FACE_BOUND(.F.)` (OpenCASCADE style).
//!
//! Every variant imports one valid solid with the cube's six faces.

use popkern_io::{StepImportOptions, import_step};

const CUBE: &str = include_str!("fixtures/cube.step");

fn solids(text: &str) -> Vec<popkern_core::topology::BRep> {
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter().cloned())
        .collect()
}

fn one_valid_cube(text: &str) -> popkern_core::topology::BRep {
    let found = solids(text);
    assert_eq!(found.len(), 1, "one solid");
    let cube = found.into_iter().next().expect("solid");
    assert_eq!(cube.faces.len(), 6);
    cube.validate().expect("valid");
    let mesh = popkern_mesh::tessellate(&cube, &popkern_mesh::TessellationOptions::default())
        .expect("mesh");
    assert_eq!(popkern_mesh::welded_open_edges(&mesh), 0, "watertight");
    cube
}

fn replace(text: &str, from: &str, to: &str) -> String {
    assert!(text.contains(from), "fixture lacks {from}");
    text.replacen(from, to, 1)
}

#[test]
fn the_fixture_itself_imports() {
    one_valid_cube(CUBE);
}

#[test]
fn a_solid_reached_through_a_shape_representation_relationship_is_attached() {
    let text = replace(
        CUBE,
        "#501 = SHAPE_DEFINITION_REPRESENTATION(#305, #500);",
        "#502 = SHAPE_REPRESENTATION('product_rep', (), #403);\n\
         #501 = SHAPE_DEFINITION_REPRESENTATION(#305, #502);\n\
         #503 = SHAPE_REPRESENTATION_RELATIONSHIP('', '', #502, #500);",
    );
    one_valid_cube(&text);
}

#[test]
fn a_mapped_item_places_its_representation() {
    let text = replace(
        CUBE,
        "#501 = SHAPE_DEFINITION_REPRESENTATION(#305, #500);",
        "#510 = CARTESIAN_POINT('', (0.0, 0.0, 0.0));\n\
         #511 = AXIS2_PLACEMENT_3D('', #510, #22, #20);\n\
         #512 = CARTESIAN_POINT('', (10.0, 0.0, 0.0));\n\
         #513 = AXIS2_PLACEMENT_3D('', #512, #22, #20);\n\
         #514 = REPRESENTATION_MAP(#511, #500);\n\
         #515 = MAPPED_ITEM('', #514, #513);\n\
         #516 = SHAPE_REPRESENTATION('product_rep', (#515), #403);\n\
         #501 = SHAPE_DEFINITION_REPRESENTATION(#305, #516);",
    );
    let cube = one_valid_cube(&text);
    let min_x = cube
        .vertices
        .values()
        .map(|v| v.position.x)
        .fold(f64::INFINITY, f64::min);
    assert!(
        (min_x - 10.0).abs() < 1e-12,
        "moved by the mapping: {min_x}"
    );
}

#[test]
fn solids_without_a_product_import_as_one_product() {
    let text = replace(
        CUBE,
        "#501 = SHAPE_DEFINITION_REPRESENTATION(#305, #500);",
        "",
    );
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let holder = imported
        .document
        .products
        .iter()
        .find(|p| p.shapes.iter().any(|s| !s.solids.is_empty()))
        .expect("a product holding the solid");
    assert_eq!(holder.name, "geometry");
    one_valid_cube(&text);
}

#[test]
fn a_reversed_oriented_face_in_the_shell_is_the_face_reversed() {
    // The top face written reversed (sense and bound flipped), then
    // wrapped in ORIENTED_FACE(.F.): the shell's face is the original.
    let text = replace(
        CUBE,
        "#115 = FACE_OUTER_BOUND('', #114, .T.);",
        "#115 = FACE_OUTER_BOUND('', #114, .F.);",
    );
    let text = replace(
        &text,
        "#116 = ADVANCED_FACE('', (#115), #83, .T.);",
        "#116 = ADVANCED_FACE('', (#115), #83, .F.);\n#117 = ORIENTED_FACE('', *, #116, .F.);",
    );
    let text = replace(
        &text,
        "#200 = CLOSED_SHELL('', (#106, #116,",
        "#200 = CLOSED_SHELL('', (#106, #117,",
    );
    one_valid_cube(&text);
}

#[test]
fn a_vertex_loop_bound_is_skipped() {
    let text = replace(
        CUBE,
        "#106 = ADVANCED_FACE('', (#105), #81, .T.);",
        "#107 = VERTEX_LOOP('', #11);\n#108 = FACE_BOUND('', #107, .T.);\n#106 = ADVANCED_FACE('', (#105, #108), #81, .T.);",
    );
    one_valid_cube(&text);
}

#[test]
fn outer_loops_marked_only_as_reversed_face_bounds_import() {
    // Every outer bound written FACE_BOUND(.F.) over its loop reversed
    // (order and edge senses): the same faces, OpenCASCADE's spelling.
    let mut text = CUBE.to_string();
    for face in ["#10", "#11", "#12", "#13", "#14", "#15"] {
        let bound = format!("{face}5 = FACE_OUTER_BOUND('', {face}4, .T.);");
        text = replace(
            &text,
            &bound,
            &format!("{face}5 = FACE_BOUND('', {face}4, .F.);"),
        );
        let e = |k: u32| format!("{face}{k}");
        // Reverse the loop: (#x0,#x1,#x2,#x3) -> (#x3,#x2,#x1,#x0), each edge flipped.
        let old_loop = format!(
            "{face}4 = EDGE_LOOP('', ({}, {}, {}, {}));",
            e(0),
            e(1),
            e(2),
            e(3)
        );
        let new_loop = format!(
            "{face}4 = EDGE_LOOP('', ({}, {}, {}, {}));",
            e(3),
            e(2),
            e(1),
            e(0)
        );
        text = replace(&text, &old_loop, &new_loop);
        for k in 0..4u32 {
            let id = e(k);
            let line_start = text
                .find(&format!("{id} = ORIENTED_EDGE("))
                .expect("oriented edge");
            let line_end = line_start + text[line_start..].find(';').expect("end");
            let line = text[line_start..line_end].to_string();
            let flipped = if line.ends_with(".T.)") {
                line.replace(".T.)", ".F.)")
            } else {
                line.replace(".F.)", ".T.)")
            };
            text.replace_range(line_start..line_end, &flipped);
        }
    }
    one_valid_cube(&text);
}

#[test]
fn placements_without_a_reference_direction_complete_their_frame() {
    // Every plane written with its ref_direction omitted: the X-normal
    // planes need ISO 10303-42's first_proj_axis (+Y), not +X.
    let mut text = CUBE.to_string();
    for (id, point, axis) in [
        ("#80", "#1", "#25"),
        ("#82", "#5", "#22"),
        ("#84", "#1", "#24"),
        ("#86", "#4", "#21"),
        ("#88", "#1", "#23"),
        ("#90", "#2", "#20"),
    ] {
        let line_start = text
            .find(&format!("{id} = AXIS2_PLACEMENT_3D("))
            .expect("placement");
        let line_end = line_start + text[line_start..].find(';').expect("end");
        text.replace_range(
            line_start..line_end,
            &format!("{id} = AXIS2_PLACEMENT_3D('', {point}, {axis}, $)"),
        );
    }
    one_valid_cube(&text);
}

#[test]
fn intersection_curve_edges_carry_their_3d_curve() {
    // Every edge's line wrapped in an INTERSECTION_CURVE (a
    // SURFACE_CURVE subtype naming the two surfaces it lies on).
    let mut text = CUBE.to_string();
    for k in 0..12u32 {
        let edge = format!("#{}", 60 + k);
        let line_start = text.find(&format!("{edge} = EDGE_CURVE(")).expect("edge");
        let line_end = line_start + text[line_start..].find(';').expect("end");
        let line = text[line_start..line_end].to_string();
        let curve = line
            .split(',')
            .nth(3)
            .expect("curve ref")
            .trim()
            .to_string();
        let wrapper = format!("#{}", 700 + k);
        text.replace_range(
            line_start..line_end,
            &line.replacen(&format!("{curve},"), &format!("{wrapper},"), 1),
        );
        text = text.replacen(
            &format!("{edge} = EDGE_CURVE("),
            &format!("{wrapper} = INTERSECTION_CURVE('', {curve}, (#81, #83), .CURVE_3D.);\n{edge} = EDGE_CURVE("),
            1,
        );
    }
    assert_eq!(text.matches("INTERSECTION_CURVE(").count(), 12);
    one_valid_cube(&text);
}

#[test]
fn an_inch_file_imports_in_millimetres() {
    let text = replace(
        CUBE,
        "#401 = (LENGTH_UNIT() NAMED_UNIT(#400) SI_UNIT(.MILLI., .METRE.));",
        "#401 = (CONVERSION_BASED_UNIT('INCH', #404) LENGTH_UNIT() NAMED_UNIT(#400));\n\
         #404 = LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4), #405);\n\
         #405 = (LENGTH_UNIT() NAMED_UNIT(#400) SI_UNIT(.MILLI., .METRE.));",
    );
    let span = |b: &popkern_core::topology::BRep| {
        let xs: Vec<f64> = b.vertices.values().map(|v| v.position.x).collect();
        xs.iter().cloned().fold(f64::MIN, f64::max) - xs.iter().cloned().fold(f64::MAX, f64::min)
    };
    let inch = one_valid_cube(&text);
    let mm = one_valid_cube(CUBE);
    assert!(
        (span(&inch) - 25.4 * span(&mm)).abs() < 1e-9,
        "{} vs {}",
        span(&inch),
        span(&mm)
    );
}

#[test]
fn quasi_uniform_edge_curves_carry_their_implicit_knots() {
    // Every edge as a linear QUASI_UNIFORM_CURVE between its vertex
    // points (knots implicit in the type, flags written .U.).
    let mut text = CUBE.to_string();
    for k in 0..12u32 {
        let edge = format!("#{}", 60 + k);
        let start = text.find(&format!("{edge} = EDGE_CURVE(")).expect("edge");
        let end = start + text[start..].find(';').expect("end");
        let line = text[start..end].to_string();
        let refs: Vec<u32> = line
            .split('#')
            .skip(2)
            .filter_map(|r| r.split(|c: char| !c.is_ascii_digit()).next()?.parse().ok())
            .collect();
        let (va, vb, curve) = (refs[0], refs[1], refs[2]);
        let spline = 720 + k;
        let new_line = line.replacen(&format!("#{curve},"), &format!("#{spline},"), 1);
        text.replace_range(start..end, &format!(
            "#{spline} = QUASI_UNIFORM_CURVE('', 1, (#{}, #{}), .UNSPECIFIED., .F., .U.);\n{new_line}",
            va - 10,
            vb - 10
        ));
    }
    assert_eq!(text.matches("QUASI_UNIFORM_CURVE(").count(), 12);
    one_valid_cube(&text);
}
