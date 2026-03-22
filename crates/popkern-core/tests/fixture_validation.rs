use popkern_core::topology_v0::BRepModel;
use popkern_core::validation_v0::{euler_characteristic, validate_model};
use serde::Deserialize;
use std::fs;

#[derive(Deserialize)]
struct FixtureExpected {
    vertices: usize,
    edges: usize,
    faces: usize,
    shells: usize,
    solids: usize,
    euler_characteristic: i64,
    #[allow(dead_code)]
    genus: u32,
}

#[derive(Deserialize)]
struct Fixture {
    #[allow(dead_code)]
    name: String,
    expected: FixtureExpected,
    model: BRepModel,
}

fn load_fixture(name: &str) -> Fixture {
    let path = format!("tests/fixtures/{}.json", name);
    let data = fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("Failed to load fixture {}: {}", path, e));
    serde_json::from_str(&data)
        .unwrap_or_else(|e| panic!("Failed to parse fixture {}: {}", path, e))
}

#[test]
fn unit_cube_topology_counts() {
    let fixture = load_fixture("unit_cube");
    assert_eq!(fixture.model.vertices.len(), fixture.expected.vertices);
    assert_eq!(fixture.model.edges.len(), fixture.expected.edges);
    assert_eq!(fixture.model.faces.len(), fixture.expected.faces);
    assert_eq!(fixture.model.shells.len(), fixture.expected.shells);
    assert_eq!(fixture.model.solids.len(), fixture.expected.solids);
}

#[test]
fn unit_cube_euler_characteristic() {
    let fixture = load_fixture("unit_cube");
    assert_eq!(
        euler_characteristic(&fixture.model),
        fixture.expected.euler_characteristic
    );
}

#[test]
fn unit_cube_passes_full_validation() {
    let fixture = load_fixture("unit_cube");
    validate_model(&fixture.model, 1e-6).expect("Unit cube should pass all validation checks");
}
