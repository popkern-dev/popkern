//! Section drawings round-trip through their own SVG (PK-856).

use super::*;

fn drawing() -> SectionDrawing {
    let circle: Vec<[f64; 2]> = (0..64)
        .map(|k| {
            let a = std::f64::consts::TAU * k as f64 / 64.0;
            [5.0 + 2.0 * a.cos(), 5.0 + 2.0 * a.sin()]
        })
        .collect();
    SectionDrawing {
        bodies: vec![
            SectionBody {
                name: "plate & \"boss\"".into(),
                regions: vec![SectionArea {
                    outer: vec![[0.0, 0.0], [10.0, 0.0], [10.0, 10.0], [0.0, 10.0]],
                    holes: vec![circle.iter().rev().copied().collect()],
                }],
            },
            SectionBody {
                name: "link_6/gripper_screw_0".into(),
                regions: vec![SectionArea {
                    outer: circle.clone(),
                    holes: Vec::new(),
                }],
            },
        ],
    }
}

#[test]
fn a_drawing_read_from_its_own_svg_is_the_drawing() {
    let d = drawing();
    let svg = d.to_svg();
    assert!(svg.starts_with("<svg") && svg.contains("fill-rule=\"evenodd\""));
    assert_eq!(SectionDrawing::from_svg(&svg).expect("parses"), d);
}

#[test]
fn malformed_path_data_declines_naming_the_path() {
    let svg = "<g data-body=\"b\">\n<path d=\"M 0 0 L 1\"/>\n";
    assert!(matches!(
        SectionDrawing::from_svg(svg),
        Err(SectionDrawingError::Path { path: 0, .. })
    ));
}

#[test]
fn every_loop_is_a_closed_polyline_on_its_bodys_layer() {
    let dxf = drawing().to_dxf();
    assert_eq!(dxf.matches("0\nPOLYLINE\n").count(), 3);
    assert_eq!(dxf.matches("0\nSEQEND\n").count(), 3);
    assert!(dxf.contains("\n8\nlink_6_gripper_screw_0\n"));
    assert!(dxf.ends_with("0\nEOF\n"));
}
