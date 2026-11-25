//! Pinned scenes of record ends at face-boundary junctions: an end on
//! one operand's boundary edge kept off the other operand's junctions
//! nearby.

use super::passes;

/// Failed at 0e00e593 (W3 V cells, tangency design §4.8): a cylinder and a frustum a thousandth across, A's seam passing B's rim 1.6e-5 off; the clip's end refinement asked A's boundary first and dragged both records ending at B's rim pierce onto the seam's pierces 1.9e-5 and 5.4e-5 away (one off B's face), and the walk dropped them as grazes (`clip_endpoint_refine`: the nearer of the two faces' junctions).
#[test]
fn a_seam_passing_a_rim_keeps_the_records_on_the_rim_pierce() {
    passes(
        r#"{"cell":{"a":"Y","b":"F","config":"V"},"seed":1,"index":10116234,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[0.8786776923474384,0.02175520311131477,0.476919515336252,0.14734929290841348],[0.16693448085234833,-0.949902929195906,-0.26422964293885465,-0.30951255893519114],[0.44727887505843267,0.3117870046083113,-0.8382902073171556,0.4160878458174818],[0.0,0.0,0.0,1.0]],"scene":[[-0.00019521720172968937,0.0005583357327099852,-0.0008063196969738702,-0.0005765021684618761],[-0.00012299287541996768,0.0008017103134072905,0.0005849216408821891,-0.0020784574290143],[0.0009730174699072767,0.00021335834400273263,-0.00008783632676548792,0.0021594024201553523],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":10.0,"rung":{"offset":10.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-0.0013412531094501837,-0.0012069184194864646,0.0030300111859358683],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed with the vertex table's foot alone deciding a record's split
/// at a pave (W3 V cells, tangency design §4.8): a sphere crossing a
/// block domed by a bicubic sheet at their touch; the foot, measured to
/// the table's anchor,
/// lay 2.1e-6 from the graph node the pave found within the tolerance,
/// and the union's edge ended that far off its vertex
/// (`graph_split_node`: the nearer of the foot and a closest-point
/// search splits).
#[test]
fn a_pave_splits_its_record_at_the_point_nearest_its_node() {
    passes(
        r#"{"cell":{"a":"S","b":"W","config":"V"},"seed":1,"index":10090003,"variant_a":2,"variant_b":0,"b_scale":1.0,"place_b":[[0.6376943165683445,0.16268196495410442,0.7529146943014884,0.7656263127041572],[-0.5793281661307904,0.7454833539735719,0.3295958811572889,0.21160835227467162],[-0.5076660659764128,-0.6463661092809445,0.5696367423446688,0.1919811874615851],[0.0,0.0,0.0,1.0]],"scene":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-2.0,"rung":{"offset":-2.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[0.4347128742552636,0.36911947774329024,1.2785834067855126],"curve":[],"curve_closed":false,"third":null}"#,
    );
}
