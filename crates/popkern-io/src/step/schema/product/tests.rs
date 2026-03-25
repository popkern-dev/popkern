// #[cfg(test)] — included only in test mode; .unwrap() allowed in tests.
use super::*;
use crate::step::parser::types::Parameter;

fn p_str(s: &str) -> Parameter {
    Parameter::String(s.into())
}
fn p_ref(id: u64) -> Parameter {
    Parameter::EntityRef(id)
}
fn p_enum(s: &str) -> Parameter {
    Parameter::Enum(s.into())
}

#[test]
fn step_product_parse() {
    let params = vec![
        p_str("P001"),
        p_str("Widget"),
        p_str("A test widget"),
        Parameter::List(vec![p_ref(10)]),
    ];
    let p = StepProduct::from_params(&params).unwrap();
    assert_eq!(p.id, "P001");
    assert_eq!(p.name, "Widget");
    assert_eq!(p.frame_of_reference, vec![10]);
}

#[test]
fn product_definition_parse() {
    let params = vec![p_str("def1"), p_str("desc"), p_ref(5), p_ref(6)];
    let pd = ProductDefinition::from_params(&params).unwrap();
    assert_eq!(pd.id, "def1");
    assert_eq!(pd.formation, 5);
    assert_eq!(pd.frame_of_reference, 6);
}

#[test]
fn shape_representation_parse() {
    let params = vec![
        p_str("shape"),
        Parameter::List(vec![p_ref(1), p_ref(2)]),
        p_ref(3),
    ];
    let sr = StepShapeRepresentation::from_params(&params).unwrap();
    assert_eq!(sr.name, "shape");
    assert_eq!(sr.items, vec![1, 2]);
    assert_eq!(sr.context_of_items, 3);
}

#[test]
fn next_assembly_usage_occurrence_parse() {
    let params = vec![
        p_str("rel1"),
        p_str("asm"),
        p_str("desc"),
        p_ref(10),
        p_ref(20),
        p_str("R1"),
    ];
    let nauo = NextAssemblyUsageOccurrence::from_params(&params).unwrap();
    assert_eq!(nauo.relating, 10);
    assert_eq!(nauo.related, 20);
    assert_eq!(nauo.reference_designator, Some("R1".to_string()));
}

#[test]
fn product_definition_formation_parse() {
    let params = vec![p_str("v1"), p_str("first version"), p_ref(100)];
    let pdf = ProductDefinitionFormation::from_params(&params).unwrap();
    assert_eq!(pdf.id, "v1");
    assert_eq!(pdf.description, "first version");
    assert_eq!(pdf.of_product, 100);
}

#[test]
fn product_definition_shape_parse() {
    let params = vec![p_str("shape1"), p_str("main shape"), p_ref(50)];
    let pds = ProductDefinitionShape::from_params(&params).unwrap();
    assert_eq!(pds.name, "shape1");
    assert_eq!(pds.description, "main shape");
    assert_eq!(pds.definition, 50);
}

#[test]
fn shape_definition_representation_parse() {
    let params = vec![p_ref(10), p_ref(20)];
    let sdr = ShapeDefinitionRepresentation::from_params(&params).unwrap();
    assert_eq!(sdr.definition, 10);
    assert_eq!(sdr.used_representation, 20);
}

#[test]
fn advanced_brep_shape_representation_parse() {
    let params = vec![
        p_str("abrep"),
        Parameter::List(vec![p_ref(1), p_ref(2)]),
        p_ref(3),
    ];
    let ab = AdvancedBrepShapeRepresentation::from_params(&params).unwrap();
    assert_eq!(ab.name, "abrep");
    assert_eq!(ab.items, vec![1, 2]);
    assert_eq!(ab.context_of_items, 3);
}

#[test]
fn faceted_brep_shape_representation_parse() {
    let params = vec![p_str("fbrep"), Parameter::List(vec![p_ref(4)]), p_ref(5)];
    let fb = FacetedBrepShapeRepresentation::from_params(&params).unwrap();
    assert_eq!(fb.name, "fbrep");
    assert_eq!(fb.items, vec![4]);
    assert_eq!(fb.context_of_items, 5);
}

#[test]
fn nauo_omitted_ref_designator() {
    let params = vec![
        p_str("rel2"),
        p_str("asm"),
        p_str("desc"),
        p_ref(10),
        p_ref(20),
        Parameter::Omitted,
    ];
    let nauo = NextAssemblyUsageOccurrence::from_params(&params).unwrap();
    assert_eq!(nauo.reference_designator, None);
}

#[test]
fn representation_relationship_parse() {
    let params = vec![p_str("rr"), p_str("link"), p_ref(1), p_ref(2)];
    let rr = RepresentationRelationship::from_params(&params).unwrap();
    assert_eq!(rr.name, "rr");
    assert_eq!(rr.description, "link");
    assert_eq!(rr.rep_1, 1);
    assert_eq!(rr.rep_2, 2);
}

#[test]
fn representation_relationship_with_transformation_parse() {
    let params = vec![p_str("rrt"), p_str("xform"), p_ref(1), p_ref(2), p_ref(3)];
    let rrt = RepresentationRelationshipWithTransformation::from_params(&params).unwrap();
    assert_eq!(rrt.name, "rrt");
    assert_eq!(rrt.rep_1, 1);
    assert_eq!(rrt.rep_2, 2);
    assert_eq!(rrt.transformation, 3);
}

#[test]
fn item_defined_transformation_parse() {
    let params = vec![p_str("idt"), p_str("transform"), p_ref(10), p_ref(20)];
    let idt = ItemDefinedTransformation::from_params(&params).unwrap();
    assert_eq!(idt.name, "idt");
    assert_eq!(idt.description, "transform");
    assert_eq!(idt.transform_item_1, 10);
    assert_eq!(idt.transform_item_2, 20);
}

#[test]
fn si_unit_parse() {
    let params = vec![Parameter::Derived, Parameter::Omitted, p_enum("METRE")];
    let su = SiUnit::from_params(&params).unwrap();
    assert!(su.dimensions.is_none());
    assert!(su.prefix.is_none());
    assert_eq!(su.name, "METRE");
}
