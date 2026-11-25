//! Pinned scenes of the tangent crossings (W3 V cells, tangency design
//! §4.8): two faces crossing at their touch.

use super::passes;

/// Failed at 4785b1bc (W3 V cells, tangency design §4.8): cylinders a thousandth across crossing at their touch; the medial sampler stepped in the raw cylinder chart (u in radians, v in length, a thousand times stretched at that scale), the common's B-wall sector at the touch sampled 1.5e-8 off A, classified outside and its loop broke (`fragment_sample_interior::chart_speeds`).
#[test]
fn cylinders_a_thousandth_across_crossing_at_their_touch_keep_all_four_sectors() {
    passes(
        r#"{"cell":{"a":"Y","b":"Y","config":"V"},"seed":1,"index":7064,"variant_a":0,"variant_b":1,"b_scale":1.0,"place_b":[[0.946760748178558,-0.3187237517247331,0.045379023732201734,-0.18479623501810183],[-0.31154046477424435,-0.942570017030637,-0.12043380672847746,-0.07677959516226207],[0.08115802188708682,0.0998845988197136,-0.9916836402815163,0.4251352833734296],[0.0,0.0,0.0,1.0]],"scene":[[-0.0004414517913640861,0.00009014376021195957,-0.0008927454387429197,-0.00001917114572989931],[-0.0002540191826330988,0.0009416836689088568,0.00022069463646578986,0.0014455735417985593],[0.0008605780445467173,0.0003242005092611723,-0.0003928096982496038,0.0006922486989175569],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.0,"rung":{"offset":0.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-0.0007497806308313036,0.0016262388759384098,0.001240698523791016],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 0e00e593 (W3 V cells, tangency design §4.8): a cylinder and a frustum crossing at their touch; the march accepted the crossing four steps ahead and pushed it straight after its last sample, so the branch records ran a chord 0.16 off both surfaces into the crossing and their ends dangled in the walk (`marching_crossing::Crossing::approach`).
#[test]
fn a_cylinder_and_a_frustum_crossing_at_their_touch_reach_the_crossing_on_both_surfaces() {
    passes(
        r#"{"cell":{"a":"Y","b":"F","config":"V"},"seed":1,"index":10090000,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[0.9200135923948783,0.2759468277015659,0.2782594797848108,0.011839699829592543],[0.2461883782422694,-0.9594230224474964,0.13747271153269716,-0.050007202012238905],[0.3049037097627914,-0.057972513134701406,-0.9506171235012203,0.6237019696881888],[0.0,0.0,0.0,1.0]],"scene":[[-0.6383288481008718,-0.19103664595479689,0.7456817562369056,-2.932244993817484],[-0.39380306794476927,0.9133938496726063,-0.10310586335191621,-2.0862547620306184],[-0.6614041316468333,-0.35946721030241147,-0.6582764611908132,2.040161450901463],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.5,"rung":{"offset":0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-3.682196813787713,-2.236433091653639,1.3860767707309969],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 0e00e593 (W3 V cells, tangency design §4.8): a cylinder and a cone crossing at their touch; a fragment beside the touch sampled within the tolerance of the other operand read its side by noise, and the loops broke where kept and dropped fragments met; its negative control (the deep sample off) still fails (`cellular/fragment_sample_deep.rs`).
#[test]
fn a_cylinder_and_a_cone_crossing_at_their_touch_sample_their_fragments_clear_of_each_other() {
    passes(
        r#"{"cell":{"a":"Y","b":"C","config":"V"},"seed":1,"index":2020,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[-0.2345802851432116,-0.9518806455138263,-0.19721898113090164,-0.44949666997381466],[0.8689928773644473,-0.11440943020723415,-0.4814165154723244,-0.4899865502292794],[0.4356873522515754,-0.28431271335597624,0.8540156978136727,1.969735405850133],[0.0,0.0,0.0,1.0]],"scene":[[0.3957520397722035,-0.9178141616396256,-0.031583028827074255,2.2861756534838626],[-0.3197704406372713,-0.1699572913769341,0.9321273434474826,2.9198031555505652],[-0.8608874423005517,-0.3587919784520934,-0.36075078362721114,2.4062792339887267],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-0.5,"rung":{"offset":-0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[3.242388143820475,3.329840176755608,2.522871516234948],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 97ab1a60 (W3 V cells, tangency design §4.8): a sphere
/// crossing a filleted box's round at their touch. Two records re-anchored
/// onto the crossing's vertex from their own ends 1.1e-3 rad along the
/// fillet's chart departed it the wrong way round inside its zone; the
/// wedge matching paired them with the wrong neighbours and three regions
/// of the round walked as one loop (`chart_embedding_zone`: a run's
/// direction is read where it leaves twice its node's zone).
#[test]
fn a_sphere_crossing_a_round_at_its_touch_reads_departures_past_the_zone() {
    passes(
        r#"{"cell":{"a":"S","b":"L","config":"V"},"seed":1,"index":1009,"variant_a":2,"variant_b":0,"b_scale":1.0,"place_b":[[0.1925428437821122,0.20040734052957243,0.9606061373789785,0.3365448632480087],[0.07897669400561796,0.9725827743755888,-0.21873597964629707,-0.5623630206818578],[-0.97810527813331,0.11798154453026247,0.17143634398228855,-0.3878391730618679],[0.0,0.0,0.0,1.0]],"scene":[[0.824785271569115,-0.0830563959228049,-0.5593128738898966,-0.27999142648171205],[0.5480394931305265,-0.12612023939189929,0.8268895931047662,1.2402065708123775],[-0.13921914297932136,-0.9885319015145866,-0.05850393077442717,2.3281736578741423],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.5,"rung":{"offset":0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[0.6993793531835348,0.24447554097452273,2.2314648570547186],"curve":[],"curve_closed":false,"third":null}"#,
    );
}
