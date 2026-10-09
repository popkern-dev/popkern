//! A boss whose ring straddles a finned housing's seam (PK-746).
//!
//! A fin straddling the housing cylinder's seam puts a root corner
//! just before u = 2π, and the face's vertex box stopped there: a boss
//! rod's seam edge piercing the housing in the strip beyond it was
//! rejected, the ring carried no node at the rod's seam, the rod's
//! chart walk fell back to 3D-Newell and the union served VALID
//! without the boss wall. The check is symmetry: the same boss turned
//! off the seam must enclose the same volume.

use popkern_rhai::ScriptHost;

fn finned_with_boss(angle: &str) -> f64 {
    let source = format!(
        r#"let pi = 3.141592653589793;
let m = model();
let body = m.profile(xy()).circle(37.0).sketch("c").extrude(blind(88.0));
let pts = [];
for k in 0..24 {{
    let a = 2.0 * pi * k / 24.0;
    let ca = cos(a); let sa = sin(a); let w = 1.5;
    pts.push([35.0 * ca + w * sa, 35.0 * sa - w * ca]);
    pts.push([45.0 * ca + w * sa, 45.0 * sa - w * ca]);
    pts.push([45.0 * ca - w * sa, 45.0 * sa + w * ca]);
    pts.push([35.0 * ca - w * sa, 35.0 * sa + w * ca]);
}}
body = body.union(m.tool(|t| t.profile(xy_at(17.6)).polygon(pts).sketch("fins").extrude(blind(57.2))));
let ang = {angle};
body = body.union(m.tool(|t| t.profile(plane_at([0.0, 0.0, 8.8], [cos(ang), sin(ang), 0.0]))
    .circle(5.0).sketch("boss").extrude(blind(43.0))), #{{ label: "boss" }});
m"#
    );
    let mut host = ScriptHost::new();
    let outcome = host
        .run(&source, "boss.rhai")
        .unwrap_or_else(|e| panic!("{}", e.render()));
    let mut model = outcome.model.expect("model");
    let part = model.to_part().expect("part");
    let brep = part.brep();
    assert!(brep.validate().is_ok());
    popkern_measure::body_volume(brep, &brep.tolerance).expect("volume")
}

#[test]
fn a_boss_across_the_seam_keeps_its_wall() {
    let across = finned_with_boss("pi / 24.0");
    let turned = finned_with_boss("pi / 6.0 + pi / 24.0");
    // The boss outside the housing is ~600; losing its wall lost 724.
    assert!(
        (across - turned).abs() < 1.0,
        "seam-straddling {across} vs turned {turned}"
    );
}
