//! Pinned touches whose tolerance band the intersector returned as a
//! curve: a small closed (or window-cut) loop round a definite touch,
//! the gap's critical point inside it within the tolerance
//! (`popkern_core::geometry::intersection::touch_ring`).

use super::passes;

/// Failed at ed681fde (W3 E cells, tangency design §4.8): two cones a
/// thousandth across, half a tolerance apart at their touch. The loop
/// round the touch split each cone off a disc a tolerance deep; the
/// union kept both and the general policy's glue could not order the
/// ring of their rims, and A - B left a disc behind.
#[test]
fn two_cones_a_thousandth_across_touching_glue_at_a_point() {
    passes(
        r#"{"cell":{"a":"C","b":"C","config":"E"},"seed":1,"index":10091015,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[0.36436392976574283,-0.37321997772080534,0.8531973833268282,1.4374142090644468],[0.7756828961346531,-0.38534526248310225,-0.4998250427159352,1.2203382355281072],[0.5153202609347144,0.84392883403277,0.14909443905883013,0.028938533142183975],[0.0,0.0,0.0,1.0]],"scene":[[-0.0002756649462613431,0.0004904197409608709,0.0008267389642919967,0.0028243718922031755],[-0.00010272517721816338,0.0008401082651657572,-0.0005326027044295581,-0.000810327782165114],[-0.0009557491173776871,-0.0002317468025152764,-0.00018120994496913727,-0.00019742334367187818],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.5,"rung":{"offset":0.5,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[0.0025422735338809484,0.0004302612625012621,-0.0001646964653877919],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ed681fde (W3 E cells, tangency design §4.8): a cylinder
/// and a ball a thousandth across, half a tolerance into each other. The
/// glued body's tangent faces read as crossing: the self-intersection
/// check met the touch's band as a curve inside both trims.
#[test]
fn a_cylinder_and_a_ball_a_thousandth_across_touching_glue_valid() {
    passes(
        r#"{"cell":{"a":"Y","b":"S","config":"E"},"seed":1,"index":10109171,"variant_a":0,"variant_b":0,"b_scale":1.0,"place_b":[[0.875104274937379,0.11460957970652891,-0.4701671534952403,-1.215104458259633],[0.481851477079271,-0.29640223402992927,0.8245998239743878,-1.1882009112593557],[-0.044851555411699544,-0.9481615684584351,-0.31460765736960605,-0.758924167734284],[0.0,0.0,0.0,1.0]],"scene":[[-0.00017530984854535425,-0.0005597161447397168,0.0008099285735919642,0.0006561855619328539],[0.0008661686278829217,0.0003033737299142903,0.0003971351005253858,-0.0015426177149124155],[-0.0004679939797416435,0.000771156415596473,0.0004316241624492059,-0.0029066467417428886],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":-0.5,"rung":{"offset":-0.5,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[0.0004548445321694189,-0.0025091284350228413,-0.0033978541419483438],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ed681fde (W3 I cells, tangency design §4.8): a torus inside
/// a box, half a tolerance into its face at the touch. A - B's cavity
/// touching the box's boundary at a point glued with its faces read as
/// crossing at the touch's band.
#[test]
fn a_torus_touching_a_box_from_inside_leaves_a_cavity_at_a_point() {
    passes(
        r#"{"cell":{"a":"P","b":"T","config":"I"},"seed":1,"index":10100090,"variant_a":1,"variant_b":1,"b_scale":0.14817502936226024,"place_b":[[-0.5000405266224841,-0.8320098748433651,0.2402478717871122,0.024507241670304977],[0.4764493974308144,-0.49597106374413263,-0.7259536318638591,-0.7332439761362199],[0.7231565829216537,-0.24854028263370379,0.6444162354297385,0.1930270326775701],[0.0,0.0,0.0,1.0]],"scene":[[-0.9237718117178866,-0.3382209175598285,-0.1795891165979756,-0.5383337395924466],[0.12807607174577906,-0.7148389202557826,0.6874604264492009,-1.9363072329853026],[-0.3608907864183064,0.6120554950431659,0.7036660509544582,-0.6048090582096743],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-0.5,"rung":{"offset":-0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-0.20039003260104987,-1.2251169485882245,-1.139587157571246],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ed681fde (W3 E cells, tangency design §4.8): a cylinder
/// resting on a rounded box's edge strip a thousandth across, half a
/// tolerance apart. The glue's own touch detection knew closed forms for
/// a plane or ball against a ball only and found nothing, and the
/// validator's cylinder-pair sampler read the touch's band as crossing
/// points (`BRep::glue_touching` searching contacts; `touch_ring::
/// points_in_touch_band`).
#[test]
fn a_cylinder_on_a_rounded_box_edge_a_thousandth_across_glues_at_a_point() {
    passes(
        r#"{"cell":{"a":"Y","b":"L","config":"E"},"seed":1,"index":10097064,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[0.42970372465692364,0.788926111978082,-0.43926108279144593,-1.7351725741440287],[-0.8671502984131614,0.22488409725066605,-0.44438440877883895,1.7750359333239243],[-0.2518036317806637,0.5718590146556227,0.7807510476324501,-0.17748016323743526],[0.0,0.0,0.0,1.0]],"scene":[[0.0007072895214821853,-0.0006341543988778954,0.0003123935517665442,0.0012895626840962874],[0.0002509447371869873,-0.00018789635062130435,-0.0009495902802268724,0.0016254100832180924],[0.0006608844616721256,0.0007500287726527637,0.000026240589043194124,0.001541770507257832],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.5,"rung":{"offset":0.5,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[0.00040792431309453114,0.0011106114162846327,0.0016046635059707908],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ed681fde (W3 E cells, tangency design §4.8): a cone
/// touching a frustum's side 2.4e-5 from its cap's rim, a thousandth
/// across. The cap's plane meets the cone 1.7e-6 outside the rim, in the
/// trim test's parametric boundary band, and the glued body failed as
/// self-intersecting: a boundary reading is confirmed in space
/// (`BRep::trim_reading_in_space`).
#[test]
fn a_cone_touching_a_frustum_beside_its_cap_rim_glues_at_a_point() {
    passes(
        r#"{"cell":{"a":"C","b":"F","config":"E"},"seed":1,"index":10100105,"variant_a":1,"variant_b":1,"b_scale":1.0,"place_b":[[-0.8877543167149187,-0.3126737369449091,0.33782748168101095,-0.958428214580437],[0.321171005359311,0.10501589431070754,0.9411805603913694,-0.6941014907184521],[-0.32975969806893873,0.9440374972250986,0.007193424946040909,1.0985911597803204],[0.0,0.0,0.0,1.0]],"scene":[[-0.0009541998239879095,0.00029900720369030177,9.868538024908116e-6,-0.0019799149118942825],[0.00029510277646663783,0.000935294453386809,0.000195291159005091,-0.001557474744425849],[0.00004916347448181696,0.0001892590225198696,-0.0009806956587910948,-0.0021783319510702216],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.0,"rung":{"offset":0.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-0.0018379354265238082,-0.0014392915492012542,-0.0029886176305889195],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 97fc89e6 (W3 E cells, tangency design §4.8): two cones a
/// thousandth across, B half a tolerance into A, the touch 6e-6 from
/// both rims. B's rim rests on A's cone inside the touch's patch, and
/// the curve of A's cone and B's base plane runs along it, reading
/// inside A's cone and in the band of B's rim: the glued body failed as
/// self-intersecting. A point where a face meets another at its rim,
/// in the touch band of the face across the rim, is the touch
/// (`BRep::meets_across_boundary_in_touch_band`).
#[test]
fn two_cones_touching_beside_both_rims_glue_at_a_point() {
    passes(
        r#"{"cell":{"a":"C","b":"C","config":"E"},"seed":1,"index":10106164,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[-0.7414219975670058,-0.39084995399244177,0.545462863069391,-1.2898387682741608],[-0.3654714199437632,-0.4465285526495209,-0.8167269389906355,-1.4405276704766483],[0.5627824292993167,-0.804890405665784,0.18822160380553923,-0.00085467219853963],[0.0,0.0,0.0,1.0]],"scene":[[-0.00042143371851566605,-0.0004963401588341464,0.0007589730348481115,-0.00033359407184290293],[0.00012152865999324863,-0.0008602863624425492,-0.0004951142892258436,-0.00005780012337828478],[0.0008986792562968748,-0.00011642088010258806,0.0004228732351164295,0.0022459068265808426],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":-0.5,"rung":{"offset":-0.5,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-0.0006900849645552551,0.00026046978671402903,0.0008620654471717558],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 97fc89e6 (W3 E cells, tangency design §4.8): the mirror
/// of the cones beside both rims, A's rim resting on B's cone.
#[test]
fn two_cones_touching_beside_both_rims_mirrored_glue_at_a_point() {
    passes(
        r#"{"cell":{"a":"C","b":"C","config":"E"},"seed":1,"index":10092022,"variant_a":2,"variant_b":1,"b_scale":1.0,"place_b":[[-0.3418122057022926,0.020701852205474495,0.9395402329587575,-0.5061101256089366],[0.7348536599781256,-0.6172995794895445,0.28094719713637484,1.7707062722684213],[0.5857939180716092,0.7864557600255047,0.1957876989828997,-0.49193858740476126],[0.0,0.0,0.0,1.0]],"scene":[[-0.00017916158254577797,-0.0008533905750960434,-0.0004895157338400241,-0.0020187190391982323],[0.0002603643476820839,-0.0005209545867298513,0.0008129063445571162,-0.0006851711682559048],[-0.0009487420797012092,0.00001818914243098202,0.0003155275285959217,0.0017076206780945419],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.0,"rung":{"offset":0.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-0.001844416900160076,-0.002404675724810211,0.0023112453087781766],"curve":[],"curve_closed":false,"third":null}"#,
    );
}
