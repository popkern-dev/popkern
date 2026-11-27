//! Junction identity, mesh winding and chart reach at a millimetre
//! (generated corpus E, V and X cells the cone window exposed).

use super::passes;

/// Failed with the apex window alone (the cone window batch before its
/// junction and orientation fixes; generated corpus CC-E 10): two cones a
/// millimetre across touching beside one's rim. A disc piece's end sat on
/// its rim crossing and the shallow-crossing pin skipped it and moved it
/// onto the next crossing 1.2e-4 along, both ends on one point
/// (`clip_junction_pins`: the nearest junction decides); then the union's
/// mesh voted the cones' shared edge `same` past its lens and turned three
/// faces inside out, the volume read 7.18e-9 for 6.13e-9
/// (`tessellate_orientation_parity`: strongest votes first).
#[test]
fn two_cones_touching_beside_a_rim_keep_their_junctions_and_winding() {
    passes(
        r#"{"cell":{"a":"C","b":"C","config":"E"},"seed":1,"index":10093049,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[-0.6730817507392207,0.42673480211111076,0.604035069750935,-0.9673519052853453],[0.3454676777660996,-0.5407398768747742,0.766976185534044,1.4709385207548433],[0.6539212800029712,0.7249122665706833,0.216539061914082,0.2771311823621154],[0.0,0.0,0.0,1.0]],"scene":[[0.000967015654731822,0.00022820653771181746,0.00011314813144356496,0.0013518007030481449],[-0.00024499123838479944,0.000711713784057529,0.0006583637161153365,0.0023400940836439876],[0.00006971381942102854,-0.0006643683208342327,0.0007441469731535796,0.0015333046998847238],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":-2.0,"rung":{"offset":-2.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[0.00044587583105720437,0.0020251547571536656,0.0008917677278164682],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed with the apex window and the junction fix alone (generated
/// corpus YC-V 977): the union and common's volumes missed the operands'
/// sum by 7.5e-4, a weak orientation vote beside the tangency winning
/// over the cycle the stronger votes closed
/// (`tessellate_orientation_parity`).
#[test]
fn a_cylinder_crossing_a_cone_at_their_tangency_meshes_one_winding() {
    passes(
        r#"{"cell":{"a":"Y","b":"C","config":"V"},"seed":1,"index":5046,"variant_a":2,"variant_b":2,"b_scale":1.0,"place_b":[[0.7973302720353294,-0.3092448129286077,0.5182972920755211,1.1044708946580606],[-0.5961447678985267,-0.2694797028848243,0.7563015968780702,1.23838946314861],[-0.09421174556964564,-0.9120023768663419,-0.3992190020363048,-1.119517007265777],[0.0,0.0,0.0,1.0]],"scene":[[0.13733255819496185,-0.9435604660717679,-0.301385824693292,0.12196302187974206],[0.8664622480944743,-0.03300477723659353,0.4981504364212007,-1.222928313246795],[-0.47998222996971185,-0.3295517130072171,0.8130268921550576,-2.0031997480746293],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":10.0,"rung":{"offset":10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-0.8596629820394848,-1.8955562281706575,-2.387326607634134],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed with the apex window, junction and orientation fixes alone
/// (generated corpus YF-X3 1111): the common's shell open 1.4e-6 at an
/// edge; closed by the overwrap slack and the exact-rim guard
/// (`chart_branches`, `exact_rims`).
#[test]
fn three_bodies_round_a_frustum_close_their_shells() {
    passes(
        r#"{"cell":{"a":"Y","b":"F","config":"X3"},"seed":1,"index":10117244,"variant_a":2,"variant_b":1,"b_scale":1.0,"place_b":[[-0.4725205717981832,-0.6027669788663262,0.6429590021267917,0.3155024756559207],[0.8434122871047055,-0.0976137871893612,0.5283249591978025,-0.049872763739010895],[-0.2556951763102459,0.7919239343137141,0.5545055987749449,-0.7935316153932199],[0.0,0.0,0.0,1.0]],"scene":[[0.0009465966901826667,0.0002644322153627008,-0.0001844730593165159,-0.00013509005543094866],[-0.00020906556353906557,0.0009389711955306709,0.0002731751894037884,-0.0020910151277935152],[0.0002454512095657984,-0.0002200197660258663,0.0009441106959888812,-0.0028308454157583637],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.0,"rung":{"offset":0.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":{"variant":0,"scale":0.6,"place":[[0.5044196035183867,-0.7237478308590624,-0.4709033243810836,0.7845573561977096],[-0.8024532131581461,-0.19157719004670826,-0.5651256682773971,0.7050596360952348],[0.3187941409099625,0.6629383512670011,-0.6774090626356201,-0.34345985186912653],[0.0,0.0,0.0,1.0]]}}"#,
    );
}

/// Failed with the apex window, junction and orientation fixes alone
/// (generated corpus YF-X 1121): a millimetre cylinder's wall spanned its
/// period by 6.2e-4 rad, two runs reading one junction 5.6e-7 apart in
/// space, and the chart walk declined it on a relative margin; the 3D
/// walk's loops did not close (`chart_branches`: the bounds' reach is
/// the slack). Served, its frustum rim promoted to a circle 1.5e-6 off a
/// seam hit's vertex left the loop open (`exact_rims`: a vertex the
/// circle passes off its tolerance keeps the fit).
#[test]
fn a_cylinder_wall_spanning_its_period_by_a_junctions_slop_walks_its_chart() {
    passes(
        r#"{"cell":{"a":"Y","b":"F","config":"X"},"seed":1,"index":10106144,"variant_a":2,"variant_b":1,"b_scale":1.0,"place_b":[[-0.5465478800088641,0.5800071131499173,0.6040506299585453,1.765317359307774],[0.364175158588279,-0.4849216377568426,0.7951273225731985,-0.7629984703015329],[0.7540967237198541,0.6545553864500089,0.053808710672727744,-1.2702362383859271],[0.0,0.0,0.0,1.0]],"scene":[[0.00030738662621504053,0.0006270106274807048,0.0007158010443205489,0.001079033053354661],[0.0009282470622775168,-0.000363144883438056,-0.00008051822778676348,0.0006302415942564732],[0.0002094537022764117,0.0006891904429539153,-0.0006936466535230243,-0.0017367456848948284],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.0,"rung":{"offset":0.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#,
    );
}
