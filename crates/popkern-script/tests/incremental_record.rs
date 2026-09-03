//! The record loop's append fast path (`FeatureTree::replay_appended`)
//! must be observationally identical to full replay: the recorder
//! builds each feature against the stored end-of-timeline state, and
//! a final full replay of the same document must reproduce every
//! body's state exactly. A parameter edit between records invalidates
//! the fast path (the freshness bit) — the staleness hazard this
//! suite pins is recording after `set_param` against pre-edit
//! geometry.

use popkern_script::prelude::*;
use popkern_script::scalar::p;

/// Mesh signature of every alive body: face/edge counts plus a
/// quantized signed volume.
fn signature(m: &mut Model) -> Vec<(usize, usize, i64)> {
    let bodies = m.alive_bodies();
    bodies
        .into_iter()
        .map(|b| {
            let brep = m.tree().body_brep(b.id).cloned().expect("brep");
            let mesh =
                popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::default())
                    .expect("mesh");
            let mut v6 = 0.0;
            for t in &mesh.triangles {
                let (a, bb, c) = (
                    mesh.vertices[t[0]].position,
                    mesh.vertices[t[1]].position,
                    mesh.vertices[t[2]].position,
                );
                v6 += a.x * (bb.y * c.z - bb.z * c.y) - a.y * (bb.x * c.z - bb.z * c.x)
                    + a.z * (bb.x * c.y - bb.y * c.x);
            }
            (
                brep.faces.len(),
                brep.edges.len(),
                (v6 / 6.0 * 1e9).round() as i64,
            )
        })
        .collect()
}

/// Recording (fast path) then a forced full replay reproduces every
/// body identically, two bodies and a boolean in the mix.
#[test]
fn incremental_record_equals_full_replay() {
    let mut m = Model::new();
    m.param("wall", 0.06, Unit::Meters).expect("param");
    let base = m
        .profile(Plane::XY)
        .rect(1.6, 1.0)
        .sketch("base_outline")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.6.into()))
        .expect("extrude");
    let base = base
        .fillet(&mut m, edges().straight().direction(Dir::Z), 0.15)
        .expect("fillet");
    let boss = m
        .tool(|m| {
            m.profile(Plane::xy_at(0.6))
                .circle(0.2)
                .sketch("boss_sketch")?
                .extrude(m, End::Blind(0.3.into()))
        })
        .expect("tool");
    let base = base.union(&mut m, boss).expect("union");
    let _ = base;
    let lid = m
        .profile(Plane::xy_at(1.2))
        .rect(1.6, 1.0)
        .sketch("lid_outline")
        .expect("sketch")
        .extrude(&mut m, End::Blind(0.08.into()))
        .expect("extrude");
    let _ = lid;

    let before = signature(&mut m);
    assert_eq!(before.len(), 2, "two alive bodies");
    m.replay().expect("full replay");
    let after = signature(&mut m);
    assert_eq!(before, after, "fast-path state == full replay state");
}

/// A parameter edit between records dirties the fast path: the next
/// record must see post-edit geometry — its selector only matches
/// the moved rim — and the end state must equal a full replay.
#[test]
fn a_parameter_edit_between_records_is_seen_by_the_next_record() {
    let mut m = Model::new();
    m.param("depth", 0.4, Unit::Meters).expect("param");
    let block = m
        .profile(Plane::XY)
        .rect(1.0, 1.0)
        .sketch("block")
        .expect("sketch")
        .extrude(&mut m, End::Blind(p("depth")))
        .expect("extrude");
    // Edit the parameter WITHOUT an explicit replay, then record a
    // feature whose selector only matches the post-edit geometry:
    // the top rim now sits at 0.7, not 0.4.
    m.set_param("depth", 0.7).expect("set");
    let block = block
        .fillet(&mut m, edges().straight().at_height(0.7), 0.05)
        .expect("fillet sees the edited depth");
    let _ = block;
    let before = signature(&mut m);
    m.replay().expect("full replay");
    let after = signature(&mut m);
    assert_eq!(before, after);
}
