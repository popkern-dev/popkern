//! Swept STEP surfaces and the ellipse (PK-636): exact kernel forms.

use super::*;
use crate::step::parser::entities::parse_step_file;
use popkern_core::geometry::SurfaceEval;
use popkern_core::geometry::curve::CurveEval;

fn cache_of(data: &str) -> EntityGraph {
    let text =
        format!("ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n{data}\nENDSEC;\nEND-ISO-10303-21;\n");
    let parsed = parse_step_file(text.as_bytes()).expect("parses");
    EntityGraph::build(parsed.entities, &mut DiagnosticLog::new())
}

fn on_surface(surface: &SurfaceGeometry, p: Point3) -> f64 {
    let (u, v) =
        popkern_core::geometry::uv_projection::project_point_to_uv(surface, &p).expect("projects");
    surface.point_at(u, v).distance_to(&p)
}

const FRAME: &str = "#1=CARTESIAN_POINT('',(0.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=DIRECTION('',(1.,0.,0.));\n\
#4=AXIS2_PLACEMENT_3D('',#1,#2,#3);\n";

#[test]
fn an_ellipse_follows_its_placement_frame() {
    let g = cache_of(&format!("{FRAME}#10=ELLIPSE('',#4,3.,1.);"));
    let mut c = GeometryCache::new(&g);
    let e = c
        .resolve_curve(10, &mut DiagnosticLog::new())
        .expect("ellipse");
    assert!(
        e.geometry
            .point_at(0.0)
            .distance_to(&Point3::new(3.0, 0.0, 0.0))
            < 1e-12
    );
    assert!(
        e.geometry
            .point_at(std::f64::consts::FRAC_PI_2)
            .distance_to(&Point3::new(0.0, 1.0, 0.0))
            < 1e-12
    );
}

#[test]
fn a_revolved_line_profile_is_its_cylinder() {
    let g = cache_of(
        "#1=CARTESIAN_POINT('',(2.,0.,0.));\n#2=CARTESIAN_POINT('',(2.,0.,5.));\n\
         #3=B_SPLINE_CURVE_WITH_KNOTS('',1,(#1,#2),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.);\n\
         #4=CARTESIAN_POINT('',(0.,0.,0.));\n#5=DIRECTION('',(0.,0.,1.));\n#6=AXIS1_PLACEMENT('',#4,#5);\n\
         #7=SURFACE_OF_REVOLUTION('',#3,#6);",
    );
    let mut c = GeometryCache::new(&g);
    let s = c
        .resolve_surface(7, true, &mut DiagnosticLog::new())
        .expect("revolution")
        .geometry;
    for k in 0..12 {
        let a = k as f64 * std::f64::consts::TAU / 12.0;
        let p = Point3::new(2.0 * a.cos(), 2.0 * a.sin(), 0.4 * k as f64);
        assert!(on_surface(&s, p) < 1e-9, "angle {a}");
    }
}

#[test]
fn an_extruded_ellipse_spans_every_point_of_the_file() {
    let g = cache_of(&format!(
        "{FRAME}#10=ELLIPSE('',#4,3.,1.);\n#11=DIRECTION('',(0.,0.,1.));\n#12=VECTOR('',#11,2.);\n\
         #13=SURFACE_OF_LINEAR_EXTRUSION('',#10,#12);\n#14=CARTESIAN_POINT('',(0.,0.,40.));"
    ));
    let mut c = GeometryCache::new(&g);
    let s = c
        .resolve_surface(13, true, &mut DiagnosticLog::new())
        .expect("extrusion")
        .geometry;
    for z in [-35.0, 0.0, 12.5, 39.0] {
        for k in 0..8 {
            let a = k as f64 * std::f64::consts::TAU / 8.0;
            let p = Point3::new(3.0 * a.cos(), a.sin(), z);
            assert!(on_surface(&s, p) < 1e-9, "z {z} angle {a}");
        }
    }
}

#[test]
fn a_swept_unbounded_line_declines_with_a_reason() {
    let g = cache_of(
        "#1=CARTESIAN_POINT('',(2.,0.,0.));\n#2=DIRECTION('',(0.,0.,1.));\n#3=VECTOR('',#2,1.);\n#4=LINE('',#1,#3);\n\
         #5=CARTESIAN_POINT('',(0.,0.,0.));\n#6=AXIS1_PLACEMENT('',#5,#2);\n#7=SURFACE_OF_REVOLUTION('',#4,#6);",
    );
    let mut c = GeometryCache::new(&g);
    let mut log = DiagnosticLog::new();
    assert!(c.resolve_surface(7, true, &mut log).is_none());
    assert!(
        log.entries()
            .iter()
            .any(|e| e.message.contains("unbounded LINE"))
    );
}
