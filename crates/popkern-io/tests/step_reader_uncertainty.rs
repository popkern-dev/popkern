//! PK-639, PK-644: the file's geometric uncertainty is read (writers type it,
//! `LENGTH_MEASURE(1.E-06)`, which once never parsed, so every import
//! ran at the default) and is the trim-recovery allowance when looser
//! than the kernel default (reported); the imported body keeps the default
//! tolerance with per-edge tolerances where imprecision is measured.
//! Angular uncertainties are skipped and a value in another length
//! unit is converted.

use popkern_io::{StepImportOptions, import_step};

const CUBE: &str = include_str!("fixtures/cube.step");

/// The mapping allowance the import reported (mm), or `None` when it
/// kept the default; the body tolerance is always the default.
fn linear_tolerance(text: &str) -> Option<f64> {
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let solid = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .next()
        .expect("a solid");
    solid.validate().expect("valid");
    let default = popkern_core::tolerance::Tolerance::default().linear();
    assert!(
        (solid.tolerance.linear() - default).abs() < 1e-18,
        "the body keeps the default"
    );
    imported.diagnostics.entries().iter().find_map(|e| {
        e.message
            .strip_prefix("trim allowance widened to the file's uncertainty ")?
            .trim_end_matches(" mm")
            .parse::<f64>()
            .ok()
    })
}

fn with_uncertainty(entity: &str) -> String {
    let old = "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-06), #401, 'tolerance');";
    assert!(CUBE.contains(old));
    CUBE.replacen(old, entity, 1)
}

#[test]
fn a_looser_file_uncertainty_is_the_import_tolerance() {
    let text = with_uncertainty(
        "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-03), #401, 'tolerance');",
    );
    assert!((linear_tolerance(&text).expect("loosened") - 1e-3).abs() < 1e-15);
}

#[test]
fn a_tighter_file_uncertainty_keeps_the_kernel_default() {
    let text = with_uncertainty(
        "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-12), #401, 'tolerance');",
    );
    assert_eq!(linear_tolerance(&text), None);
}

#[test]
fn an_angular_uncertainty_is_not_a_length() {
    let text = CUBE
        .replacen(
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#402))",
            "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#404, #402))",
            1,
        )
        .replacen(
            "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-06), #401, 'tolerance');",
            "#404 = UNCERTAINTY_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.5), #405, 'angle');\n\
             #405 = (NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($, .RADIAN.));\n\
             #402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(2.0E-03), #401, 'tolerance');",
            1,
        );
    assert!((linear_tolerance(&text).expect("loosened") - 2e-3).abs() < 1e-15);
}

#[test]
fn an_uncertainty_in_another_length_unit_is_converted() {
    // 1e-6 metre = 1e-3 mm.
    let text = with_uncertainty(
        "#406 = (LENGTH_UNIT() NAMED_UNIT(#400) SI_UNIT($, .METRE.));\n\
         #402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.0E-06), #406, 'tolerance');",
    );
    assert!((linear_tolerance(&text).expect("loosened") - 1e-3).abs() < 1e-12);
}

/// The cube with edge E0's line moved 1e-5 off both its vertices (and
/// off one of its faces by as much).
fn cube_with_an_edge_off_its_vertices(uncertainty: &str) -> String {
    cube_with_an_edge_moved(uncertainty, "1.0E-05")
}

/// The cube with edge E0's line moved `offset` off its vertices and
/// one of its faces.
fn cube_with_an_edge_moved(uncertainty: &str, offset: &str) -> String {
    let text = with_uncertainty(&format!(
        "#402 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({uncertainty}), #401, 'tolerance');"
    ));
    let text = text.replacen(
        "#40 = LINE('', #1, #30);",
        &format!("#409 = CARTESIAN_POINT('', (0.0, {offset}, 0.0));\n#40 = LINE('', #409, #30);"),
        1,
    );
    assert!(text.contains("#409"));
    text
}

/// A radius-1 cylinder of height 2 whose two rim circles are written at
/// radius `1 + miss` (the caps hold them exactly; the cylinder misses
/// them by `miss`), under the file uncertainty `uncertainty`.
fn cylinder_off_its_rims(miss: f64, uncertainty: f64) -> String {
    let r = 1.0 + miss;
    format!(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('rod'),'2;1');\nFILE_NAME('rod.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n\
#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n\
#5=CYLINDRICAL_SURFACE('',#4,1.);\n#6=CIRCLE('',#4,{r:.17e});\n#7=CARTESIAN_POINT('',(0.,0.,2.));\n#8=AXIS2_PLACEMENT_3D('',#7,#2,#3);\n#9=CIRCLE('',#8,{r:.17e});\n\
#10=CARTESIAN_POINT('',({r:.17e},0.,0.));\n#11=VERTEX_POINT('',#10);\n#12=EDGE_CURVE('',#11,#11,#6,.T.);\n\
#13=CARTESIAN_POINT('',({r:.17e},0.,2.));\n#14=VERTEX_POINT('',#13);\n#15=EDGE_CURVE('',#14,#14,#9,.T.);\n\
#20=ORIENTED_EDGE('',*,*,#12,.T.);\n#21=EDGE_LOOP('',(#20));\n#22=FACE_BOUND('',#21,.T.);\n\
#23=ORIENTED_EDGE('',*,*,#15,.F.);\n#24=EDGE_LOOP('',(#23));\n#25=FACE_BOUND('',#24,.T.);\n#26=ADVANCED_FACE('',(#22,#25),#5,.T.);\n\
#30=DIRECTION('',(0.,0.,-1.));\n#31=AXIS2_PLACEMENT_3D('',#1,#30,#3);\n#32=PLANE('',#31);\n#33=ORIENTED_EDGE('',*,*,#12,.F.);\n#34=EDGE_LOOP('',(#33));\n#35=FACE_OUTER_BOUND('',#34,.T.);\n#36=ADVANCED_FACE('',(#35),#32,.T.);\n\
#40=PLANE('',#8);\n#41=ORIENTED_EDGE('',*,*,#15,.T.);\n#42=EDGE_LOOP('',(#41));\n#43=FACE_OUTER_BOUND('',#42,.T.);\n#44=ADVANCED_FACE('',(#43),#40,.T.);\n\
#50=CLOSED_SHELL('',(#26,#36,#44));\n#51=MANIFOLD_SOLID_BREP('',#50);\n\
#52=(LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.));\n#53=(NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.));\n\
#54=UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({uncertainty:.17e}),#52,'tolerance');\n\
#55=(GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#54)) GLOBAL_UNIT_ASSIGNED_CONTEXT((#52,#53)) REPRESENTATION_CONTEXT('',''));\n\
#56=ADVANCED_BREP_SHAPE_REPRESENTATION('',(#51),#55);\nENDSEC;\nEND-ISO-10303-21;\n"
    )
}

/// PK-644: a face whose surface misses its boundary at the file's
/// precision records that as its face tolerance, and the mesher's
/// boundary residual gate honours it — the face no longer falls to
/// the fallback grid, whose boundary samples the surface and cracks
/// against its neighbours (measured: weldment_asm_solid body 1's
/// cylinders, 199 open mesh edges).
#[test]
fn a_face_off_its_edges_within_the_file_precision_meshes_watertight() {
    for miss in [3e-3, 6e-3] {
        let text = cylinder_off_its_rims(miss, 1e-2);
        let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
        let solid = imported
            .document
            .products
            .iter()
            .flat_map(|p| p.shapes.iter())
            .flat_map(|s| s.solids.iter())
            .next()
            .expect("a solid");
        solid
            .validate()
            .expect("valid within the recorded tolerances");
        assert!(
            solid.faces.values().any(|f| f.tolerance.is_some()),
            "miss {miss}: the miss is recorded on the cylinder"
        );
        let mesh = popkern_mesh::tessellate(solid, &popkern_mesh::TessellationOptions::default())
            .expect("mesh");
        assert_eq!(
            popkern_mesh::welded_open_edges(&mesh),
            0,
            "miss {miss}: watertight"
        );
        let volume: f64 = mesh
            .triangles
            .iter()
            .map(|t| {
                let (a, b, c) = (
                    mesh.vertices[t[0]].position,
                    mesh.vertices[t[1]].position,
                    mesh.vertices[t[2]].position,
                );
                (a.x * (b.y * c.z - b.z * c.y) - a.y * (b.x * c.z - b.z * c.x)
                    + a.z * (b.x * c.y - b.y * c.x))
                    / 6.0
            })
            .sum();
        let exact = 2.0 * std::f64::consts::PI;
        assert!(
            (volume.abs() - exact).abs() < 1e-2 * exact,
            "miss {miss}: volume {volume}"
        );
    }
}

/// Beyond the file's uncertainty the miss is a defect: nothing is
/// recorded on the face, and the strict gate keeps it visible.
#[test]
fn a_face_off_its_edges_beyond_the_file_precision_records_nothing() {
    let text = cylinder_off_its_rims(3e-3, 1e-3);
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let solid = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .next()
        .expect("a solid");
    assert!(
        solid.faces.values().all(|f| f.tolerance.is_none()),
        "no face tolerance beyond the declared precision"
    );
}

#[test]
fn the_same_gap_beyond_the_file_precision_is_not_absorbed() {
    // Written to 1e-6, a 1e-5 gap is a defect: nothing records it.
    let text = cube_with_an_edge_off_its_vertices("1.0E-06");
    let imported = import_step(text.as_bytes(), &StepImportOptions::default()).expect("import");
    let solid = imported
        .document
        .products
        .iter()
        .flat_map(|p| p.shapes.iter())
        .flat_map(|s| s.solids.iter())
        .next()
        .expect("a solid");
    assert!(
        solid.validate().is_err(),
        "a gap beyond the declared precision must fail validation"
    );
}
