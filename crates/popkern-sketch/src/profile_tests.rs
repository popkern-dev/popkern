use super::*;

fn path_profile() -> AnalyticProfile<f64> {
    AnalyticProfile {
        shape: ProfileShape::Path {
            start: [0.0, 0.0],
            segments: vec![
                PathSegment::Line { to: [2.0, 0.0] },
                PathSegment::Arc {
                    via: [3.0, 1.0],
                    to: [2.0, 2.0],
                },
                PathSegment::Line { to: [0.0, 2.0] },
            ],
        },
        center: [0.5, -0.5],
    }
}

#[test]
fn scalars_enumerate_placement_then_shape_in_order() {
    let p = path_profile();
    let s: Vec<f64> = p.scalars().into_iter().copied().collect();
    assert_eq!(
        s,
        vec![0.5, -0.5, 0.0, 0.0, 2.0, 0.0, 3.0, 1.0, 2.0, 2.0, 0.0, 2.0]
    );
    let rect = AnalyticProfile::new(ProfileShape::Rect {
        width: 2.0,
        height: 1.0,
    });
    assert_eq!(rect.scalars(), vec![&0.0, &0.0, &2.0, &1.0]);
    let poly = AnalyticProfile::new(ProfileShape::Polygon {
        points: vec![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]],
    });
    assert_eq!(poly.scalars().len(), 8);
}

#[test]
fn map_scalars_rebuilds_every_arm_and_propagates_errors() {
    let p = path_profile();
    let doubled = p.map_scalars(|x| Ok::<_, ()>(x * 2.0)).unwrap();
    assert_eq!(doubled.center, [1.0, -1.0]);
    assert_eq!(
        doubled.scalars().into_iter().copied().collect::<Vec<_>>(),
        p.scalars().into_iter().map(|x| x * 2.0).collect::<Vec<_>>()
    );
    let circle = AnalyticProfile::new(ProfileShape::Circle { radius: 3.0 });
    let named = circle.map_scalars(|x| Ok::<_, ()>(format!("{x}"))).unwrap();
    assert_eq!(named.shape, ProfileShape::Circle { radius: "3".into() });
    let failed = p.map_scalars(|x| if *x == 3.0 { Err("via") } else { Ok(*x) });
    assert_eq!(failed.unwrap_err(), "via");
}

#[test]
fn descriptions_name_shape_and_placement() {
    assert_eq!(path_profile().describe(), "path[3 segs] at (0.5, -0.5)");
    assert_eq!(
        AnalyticProfile::new(ProfileShape::Circle { radius: 0.35 }).describe(),
        "circle r=0.35"
    );
    assert_eq!(
        AnalyticProfile::new(ProfileShape::Rect {
            width: 2.0,
            height: 1.0
        })
        .describe(),
        "rect 2×1"
    );
}

#[test]
fn serde_round_trips_generic_scalars() {
    let p = path_profile();
    let json = serde_json::to_string(&p).unwrap();
    let back: AnalyticProfile<f64> = serde_json::from_str(&json).unwrap();
    assert_eq!(back, p);
    let strings = p.map_scalars(|x| Ok::<_, ()>(x.to_string())).unwrap();
    let json = serde_json::to_string(&strings).unwrap();
    let back: AnalyticProfile<String> = serde_json::from_str(&json).unwrap();
    assert_eq!(back, strings);
}

#[test]
fn errors_spell_the_defect() {
    assert_eq!(
        ProfileError::PolygonTooShort { got: 2 }.to_string(),
        "polygon profile needs at least 3 points, got 2"
    );
    assert_eq!(
        ProfileError::NonPositiveDimension {
            what: "circle radius",
            value: -1.0
        }
        .to_string(),
        "circle radius must be positive and finite, got -1"
    );
}
