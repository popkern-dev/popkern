//! Pinned scenes at a thousandfold scale: lengths judged against an
//! absolute tolerance a millionth of the scene.

use super::passes;

/// Failed at 0e00e593 (W3 V cells, tangency design §4.8): a cylinder and a frustum crossing at their touch a thousand units across. The march's step ceiling (16 x 1000 tolerances) ran a 7331-long branch out of its step budget mid-face; two record ends at the crossing 8.5e-7 apart read 4.3e-6 apart along the cone through their uvs' misses and stayed two vertices; and a record split at a seam pierce at a closest point found again, 6.4e-7 off the foot the vertex table measured, left an edge end 1.25e-6 off its vertex (`marching` step ceiling, `face_path`, `graph_split_node::split_at_existing_node_at`: the nearer of the foot and the search).
#[test]
fn a_cylinder_and_a_frustum_a_thousand_across_crossing_at_their_touch_close() {
    passes(
        r#"{"cell":{"a":"Y","b":"F","config":"V"},"seed":1,"index":10109172,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[0.94012702818585,-0.1446154384746749,0.3086220112519836,0.007070249638308823],[0.15753937959016975,0.987362287960649,-0.01723531814203104,-0.05431195835952918],[-0.30222924205441526,0.06482350860521319,0.951028599979687,0.20708846468356246],[0.0,0.0,0.0,1.0]],"scene":[[976.6812539682917,-176.07719332150074,-122.84360031742953,1653.6333590130673],[105.49595083789008,891.9325243304309,-439.6896364463858,1608.8147945467624],[186.98751964699412,416.47712306148134,889.7092072484676,-2846.250154004853],[0.0,0.0,0.0,1.0]],"scale":1000.0,"offset":0.0,"rung":{"offset":0.0,"scale":1000.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[2677.2408175652217,2127.642899189366,-3189.8127852779207],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 97ab1a60 (W3 V cells, tangency design §4.8): a cylinder and
/// a domed block a thousand across, half a tolerance apart at their
/// touch. The single-crossing clip read a 340-long record's end side 2.5
/// inside the rim against the boundary's coarse octagon and dropped it
/// (`clip_uv`: each side read midway between its end and the crossing);
/// the sliver pass absorbed a plane strip 0.73 wide into the dome when
/// the NURBS projection of its far end gave up
/// (`sliver_removal::sliver_lies_on_face_surface`); and the common was
/// stamped with the bare default ceiling though both operands carry a
/// thousandfold one (`assemble_edge_geometry::carry_imprecision_ceiling`).
#[test]
fn a_cylinder_and_a_domed_block_a_thousand_across_touching_close() {
    passes(
        r#"{"cell":{"a":"Y","b":"W","config":"V"},"seed":1,"index":10115225,"variant_a":1,"variant_b":1,"b_scale":1.0,"place_b":[[0.16182311549012787,-0.10448102730920583,-0.9812731496507391,-0.13793186142592184],[0.9011028299652202,-0.38971382392788917,0.19009688389916837,-0.29895753392820523],[-0.4022772291862171,-0.9149900821166945,0.031083444248416636,-0.8298613050737804],[0.0,0.0,0.0,1.0]],"scene":[[993.0699617389604,-53.578721779185116,104.6010117746887,1385.7680395385898],[93.16064028445265,901.4648850425001,-422.70812168339387,1812.9954878012952],[-71.64597820932245,429.52343545820736,900.2091269247431,-1220.13975633533],[0.0,0.0,0.0,1.0]],"scale":1000.0,"offset":0.5,"rung":{"offset":0.5,"scale":1000.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[555.0825847224912,2837.2095027268583,-1785.9208896556704],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ffcb5603 (W3 I cells, tangency design §4.8): a sphere a
/// quarter the size of another, touching it from inside a thousand units
/// across. The general policy's contact-extent test stepped its tangent
/// probes along the chart normal unnormalized (∂u × ∂v grows with the
/// radius squared), read the point touch as an area and declined A − B
/// (`shell_contact::contact_extent`: the normal normalized first).
#[test]
fn a_sphere_touching_a_sphere_from_inside_a_thousand_across_serves_general() {
    passes(
        r#"{"cell":{"a":"S","b":"S","config":"I"},"seed":1,"index":10119261,"variant_a":2,"variant_b":0,"b_scale":0.2553629584682804,"place_b":[[0.9601278716391919,-0.10451372949737354,0.2592900893750962,1.0471294784880194],[0.24643879425263204,0.754326561945957,-0.6084894071634289,0.15118705738256255],[-0.1319939043626803,0.6481267764020725,0.7500061939222593,-0.4999679674469745],[0.0,0.0,0.0,1.0]],"scene":[[999.1373855480729,-3.41933332400682,41.38590291126784,1569.8076417749594],[2.8874287192096446,999.9125545306158,12.905272079647368,239.87132568393267],[-41.42641132843919,-12.774640960808995,999.0598885916546,2166.7212737017458],[0.0,0.0,0.0,1.0]],"scale":1000.0,"offset":0.0,"rung":{"offset":0.0,"scale":1000.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[2796.1427742128517,416.6341774470351,1514.9106898776467],"curve":[],"curve_closed":false,"third":null}"#,
    );
}
