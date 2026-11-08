use crate::geometry::transform::Transform3;
use crate::tolerance::Tolerance;
use crate::topology::{BRep, FaceHandle, HalfEdgeHandle, RegionKind};

fn cube() -> BRep {
    BRep::unit_cube(Tolerance::default()).expect("cube")
}

/// A cube is one solid region bounded by one face set, its faces
/// fronting onto one partner shell in the outside.
#[test]
fn a_cube_is_a_solid_region_inside_the_outside() {
    let cube = cube();
    let regions = cube.solid_regions();
    assert_eq!(regions.len(), 1);
    let back = cube.back_shells();
    assert_eq!(back.len(), 1);
    assert_eq!(cube.back_shell_count(), 1);
    assert_eq!(cube.region_shells(regions[0]), back);
    assert_eq!(cube.region_outer_shell(regions[0]), Some(back[0]));
    assert_eq!(cube.shell_region(back[0]), Some(regions[0]));
    assert_eq!(cube.shell_faces(back[0]).len(), 6);

    let any_face = cube.faces.keys().next().expect("a face");
    let front = cube.face_front_shell(any_face).expect("a front");
    assert_ne!(front, back[0]);
    assert_eq!(cube.shell_region(front), Some(cube.outside));
    assert_eq!(cube.region_kind(cube.outside), Some(RegionKind::Void));
    assert_eq!(cube.region_outer_shell(cube.outside), None);
    assert_eq!(cube.region_shells(cube.outside), vec![front]);
    assert_eq!(cube.shell_faces(front), cube.shell_faces(back[0]));
    for f in cube.faces.keys() {
        assert_eq!(cube.face_back_shell(f), Some(back[0]));
        assert_eq!(cube.face_front_shell(f), Some(front));
    }
    assert!(cube.cavity_regions().is_empty());
}

/// A solid region's first shell is its outer one, partnered in the
/// outside; a further shell bounds a cavity, a void region of its own.
#[test]
fn a_further_shell_of_a_solid_bounds_a_cavity() {
    let mut body = BRep::new(Tolerance::default());
    let (region, outer) = body.add_solid_region(FaceHandle::default());
    let cavity = body.add_region_shell(region, FaceHandle::default());
    assert_eq!(body.region_outer_shell(region), Some(outer.back));
    assert_eq!(body.region_shells(region), vec![outer.back, cavity.back]);
    assert_eq!(body.shell_region(cavity.back), Some(region));
    assert_eq!(body.shell_region(outer.front), Some(body.outside));
    let void = body.shell_region(cavity.front).expect("cavity region");
    assert_eq!(body.cavity_regions(), vec![void]);
    assert_eq!(body.region_outer_shell(void), Some(cavity.front));
}

/// A shell of the outside holds a sheet: both sides of its faces.
#[test]
fn a_void_regions_shell_holds_both_sides() {
    let mut body = BRep::new(Tolerance::default());
    let outside = body.outside;
    let sides = body.add_region_shell(outside, FaceHandle::default());
    assert_eq!(sides.back, sides.front);
    assert_eq!(body.shell_region(sides.back), Some(outside));
}

/// Two cubes appended share one outside, each with its own partner.
#[test]
fn appended_bodies_share_the_outside() {
    let mut a = cube();
    let b = cube();
    let map = a.append_body(&b);
    assert_eq!(map.regions.get(&b.outside), Some(&a.outside));
    assert_eq!(a.solid_regions().len(), 2);
    assert_eq!(a.regions.len(), 3);
    assert_eq!(a.region_shells(a.outside).len(), 2);
    assert_eq!(a.back_shell_count(), 2);
}

/// Keeping one of two solid regions drops the other's faces, its shell
/// and its partner, leaving one partner in the outside.
#[test]
fn retaining_a_region_drops_the_rest_with_their_partners() {
    let mut a = cube();
    let b = cube();
    a.append_body(&b);
    let keep = a.solid_regions()[0];
    a.retain_solid_regions(&[keep]);
    assert_eq!(a.solid_regions(), vec![keep]);
    assert_eq!(a.faces.len(), 6);
    assert_eq!(a.shells.len(), 2);
    assert_eq!(a.region_shells(a.outside).len(), 1);
}

/// A reversed small cube inside a big one, turned inside out into it,
/// is the big cube's cavity: its faces' backs in the big cube's region,
/// their fronts in a cavity region.
#[test]
fn an_inverted_body_becomes_the_hosts_cavity() {
    let mut host = cube();
    let mut tool = cube()
        .transform(&Transform3::scale_uniform(0.5))
        .and_then(|t| t.transform(&Transform3::translation(0.25, 0.25, 0.25)))
        .expect("small cube");
    tool.reverse_body_orientation();
    let solid = host.solid_regions()[0];
    let map = host.append_body(&tool);
    let tool_region = map.regions[&tool.solid_regions()[0]];
    host.invert_regions_into(&[tool_region], solid);

    assert_eq!(host.solid_regions(), vec![solid]);
    assert_eq!(host.cavity_regions(), vec![tool_region]);
    assert_eq!(host.region_shells(solid).len(), 2);
    assert_eq!(host.region_shells(host.outside).len(), 1);
    for &tool_face in map.faces.values() {
        let back = host.face_back_shell(tool_face).expect("back");
        let front = host.face_front_shell(tool_face).expect("front");
        assert_eq!(host.shell_region(back), Some(solid));
        assert_eq!(host.shell_region(front), Some(tool_region));
    }
    assert_eq!(host.back_shell_count(), 2);
}

/// A body written before regions were stored loads with them derived:
/// the cube's solid region and a partner in the outside.
#[test]
fn a_body_written_before_regions_loads_with_them() {
    let fresh = cube();
    let mut json = serde_json::to_value(&fresh).expect("serializes");
    let obj = json.as_object_mut().expect("an object");
    // The old form: a `solids` arena of solid regions with no kind,
    // shells owned by a `solid`, faces with one `shell`, no outside and
    // no void partners.
    let front = fresh.faces.values().next().expect("a face").front;
    let mut regions = obj.remove("regions").expect("regions");
    vacate(&mut regions, &serde_json::to_value(fresh.outside).unwrap());
    drop_in_slots(&mut regions, "kind");
    obj.insert("solids".into(), regions);
    obj.remove("outside");
    let shells = obj.get_mut("shells").expect("shells");
    vacate(shells, &serde_json::to_value(front).unwrap());
    rename_in_slots(shells, "region", "solid");
    let faces = obj.get_mut("faces").expect("faces");
    rename_in_slots(faces, "back", "shell");
    drop_in_slots(faces, "front");

    let loaded: BRep = serde_json::from_value(json).expect("old body loads");
    assert_eq!(loaded.solid_regions().len(), 1);
    assert_eq!(loaded.back_shell_count(), 1);
    assert_eq!(loaded.shells.len(), 2);
    let face = loaded.faces.keys().next().expect("a face");
    let front = loaded.face_front_shell(face).expect("front derived");
    assert_eq!(loaded.shell_region(front), Some(loaded.outside));
    loaded.validate().expect("derived regions validate");
}

/// A body round-trips with its regions as stored.
#[test]
fn a_body_round_trips_with_its_regions() {
    let json = serde_json::to_string(&cube()).expect("serializes");
    let back: BRep = serde_json::from_str(&json).expect("loads");
    assert_eq!(serde_json::to_string(&back).expect("serializes"), json);
}

/// A field the stored form does not know is refused, not dropped.
#[test]
fn an_unknown_body_field_is_refused() {
    let mut json = serde_json::to_value(cube()).expect("serializes");
    json.as_object_mut()
        .expect("an object")
        .insert("unheard_of".into(), serde_json::json!(1));
    assert!(serde_json::from_value::<BRep>(json).is_err());
}

/// A solid's shell opened by a removed face settles into the outside,
/// its solid region and partner gone; healing the hole closed would
/// fill it again.
#[test]
fn an_open_shell_settles_into_the_outside_and_fills_when_closed() {
    let mut body = cube();
    let face = body.faces.keys().next().expect("a face");
    let solid = body.solid_regions()[0];
    assert!(body.move_shell_anchor_off(face));
    let outer = body.faces[face].outer_loop;
    for fin in loop_fins(&body, body.loops[outer].half_edge) {
        let edge = body.half_edges[fin].edge;
        let survivor = body.mate(fin).expect("a cube edge has a mate");
        body.edges[edge].half_edge = survivor;
        body.remove_half_edge_reanchored(fin);
    }
    body.loops.remove(outer);
    body.faces.remove(face);
    let settled = body.settle_open_shells();
    assert_eq!(settled.len(), 1);
    assert!(!body.regions.contains_key(solid));
    assert!(body.solid_regions().is_empty());
    assert_eq!(body.shells.len(), 1);
    assert_eq!(body.shell_region(settled[0]), Some(body.outside));
    assert!(
        body.faces
            .values()
            .all(|f| f.back == settled[0] && f.front == settled[0])
    );
    // Still open: nothing to fill.
    assert!(body.fill_closed_shells().is_empty());
}

/// A closed shell lying in the outside with both sides of its faces is
/// filled into a solid region with a partner.
#[test]
fn a_closed_sheet_fills_into_a_solid() {
    let mut body = cube();
    let shell = body.back_shells()[0];
    let solid = body.solid_regions()[0];
    let front = body.face_front_shell(body.faces.keys().next().unwrap());
    // Make the cube a closed sheet in the outside.
    body.shells.remove(front.unwrap());
    body.regions.remove(solid);
    body.shells[shell].region = body.outside;
    for f in body.faces.values_mut() {
        f.front = shell;
    }
    let made = body.fill_closed_shells();
    assert_eq!(made.len(), 1);
    assert_eq!(body.solid_regions(), made);
    assert_eq!(body.region_outer_shell(made[0]), Some(shell));
    assert_eq!(body.shells.len(), 2);
    body.validate().expect("filled cube validates");
}

fn loop_fins(body: &BRep, start: HalfEdgeHandle) -> Vec<HalfEdgeHandle> {
    let mut fins = vec![start];
    let mut at = body.half_edges[start].next;
    while at != start {
        fins.push(at);
        at = body.half_edges[at].next;
    }
    fins
}

fn slots(arena: &mut serde_json::Value) -> Vec<&mut serde_json::Value> {
    arena
        .as_array_mut()
        .map(|a| a.iter_mut().collect())
        .unwrap_or_default()
}

/// Mark the slot `key` names vacant, the way a slot map writes one.
fn vacate(arena: &mut serde_json::Value, key: &serde_json::Value) {
    let index = key["idx"].as_u64().expect("a slot index") as usize;
    let slot = &mut arena[index];
    let version = slot["version"].as_u64().expect("a version");
    *slot = serde_json::json!({ "value": null, "version": version + 1 });
}

fn rename_in_slots(arena: &mut serde_json::Value, from: &str, to: &str) {
    for slot in slots(arena) {
        if let Some(value) = slot.get_mut("value").and_then(|v| v.as_object_mut())
            && let Some(x) = value.remove(from)
        {
            value.insert(to.into(), x);
        }
    }
}

fn drop_in_slots(arena: &mut serde_json::Value, field: &str) {
    for slot in slots(arena) {
        if let Some(value) = slot.get_mut("value").and_then(|v| v.as_object_mut()) {
            value.remove(field);
        }
    }
}
