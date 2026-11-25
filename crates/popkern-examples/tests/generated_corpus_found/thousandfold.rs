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
