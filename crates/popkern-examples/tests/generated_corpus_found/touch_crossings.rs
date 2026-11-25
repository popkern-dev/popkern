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

/// Failed at ffcb5603 (W3 V cells, tangency design §4.8): a domed block
/// crossing a sphere at their touch 1e-5 off the sphere's seam. The
/// crossing joined the seam's vertex, and its branch running along the
/// seam took an endpoint-refined direction a hair across it; sorted
/// through the outside half-plane, its arrival paired with its own twin
/// and the sphere's loops broke (`chart_walk_seam_side`: a ray at a seam
/// copy is read along the seam, inside); the declined differences name
/// the seam vertex, inside the touch's tolerance patch
/// (`judge_contact::on_touch`).
#[test]
fn a_domed_block_crossing_a_sphere_at_its_touch_beside_the_seam_closes() {
    passes(
        r#"{"cell":{"a":"S","b":"W","config":"V"},"seed":1,"index":15136,"variant_a":2,"variant_b":0,"b_scale":1.0,"place_b":[[-0.7660180803144425,0.6421403793274221,-0.02953021958284019,-0.00966712174222828],[0.6351580554567615,0.7631609691814385,0.11899403223033014,-0.02399352692643358],[0.09894718399106472,0.07239522329050453,-0.9924557352471477,-0.5402930007236637],[0.0,0.0,0.0,1.0]],"scene":[[0.9707613213571477,0.11433784292783923,0.2110670856186098,-2.622382300761564],[-0.11734490578307201,0.9930897056790087,0.0017348086758825591,0.9263075631544786],[-0.20941019565361843,-0.026451732438324934,0.9774700383169422,1.0767479235418982],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":{"Seam":9.999999999999999e-6},"steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-2.6093009661276017,0.8875172283809569,-0.3226534444052507],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ffcb5603 (W3 V cells, tangency design §4.8): the same pair
/// half a tolerance apart with the touch on the sphere's seam; the
/// branch along the seam turned back on itself as in the scene above
/// (`chart_walk_seam_side`).
#[test]
fn a_domed_block_crossing_a_sphere_on_its_seam_closes() {
    passes(
        r#"{"cell":{"a":"S","b":"W","config":"V"},"seed":1,"index":19175,"variant_a":2,"variant_b":0,"b_scale":1.0,"place_b":[[-0.03625402174081305,-0.8031159891985973,0.594718718219946,0.44385294582728996],[-0.9972316808840123,-0.00958622204179621,-0.07373655123623157,-0.017837043952706126],[0.06492010897144798,-0.5957455935552726,-0.8005451687510279,-0.332512035544828],[0.0,0.0,0.0,1.0]],"scene":[[-0.04612104176668141,0.8843608290478449,-0.46451993881012754,1.4018083552165033],[-0.9846623354900848,0.038081153840133206,0.1702642381400629,-2.7403356973101625],[0.16826447805038336,0.46524805186912144,0.8690404568594142,2.413656601690593],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.5,"rung":{"offset":0.5,"scale":1.0,"steer_a":{"Seam":0.0},"steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[1.784765471380901,-3.9301898913705404,1.7831085511377607],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ffcb5603 (W3 V cells, tangency design §4.8): two cones
/// crossing at a tangent touch, half a tolerance into each other, with
/// B's rim 4.4e-4 from the crossing. Two of the crossing's arms leave B's
/// face across the rim within the touch's tolerance patch; the clip kept
/// each arm's piece from the crossing to the rim, and the endpoint
/// refinement's far re-trim then moved the piece's other end, the arm's
/// own end at the crossing, onto the rim's junction: one arm became a
/// stub along the rim, dead-ending on A's cone while it split B's, and
/// every boolean declined (`clip_boundary_edges::is_open_source_end`:
/// a record's own end is no crossing to re-trim).
#[test]
fn a_cone_crossing_a_cone_beside_its_rim_closes() {
    passes(
        r#"{"cell":{"a":"C","b":"C","config":"V"},"seed":1,"index":10095045,"variant_a":2,"variant_b":1,"b_scale":1.0,"place_b":[[0.48909022849169576,0.8697944183651308,-0.06517989087751325,-0.26700692395650927],[-0.7930041636351831,0.47453979834844834,0.3820423749280669,0.619322909531233],[0.3632287775647801,-0.1351652675959234,0.921843373672837,-0.24864169190636987],[0.0,0.0,0.0,1.0]],"scene":[[0.6473771445111265,0.522207018653795,-0.5551600331736202,-0.9161053123890781],[-0.26730379638866175,0.8376897416714346,0.4762610388585802,2.3604793648239886],[0.7137587219788133,-0.15992412690760624,0.6818891115364807,-2.428106078625263],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-0.5,"rung":{"offset":-0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-1.0974600930854885,2.0440189707067415,-4.142888228331109],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ffcb5603 (W3 V cells, tangency design §4.8): a cylinder
/// crossing a ball at their touch 1e-5 from the ball's pole. The arms of
/// the crossing end at it, and the endpoint refinement's near tier
/// dragged their own ends there 1e-5 onto the junction at the pole, as
/// it corrects a clip crossing's fit slop
/// (`clip_boundary_edges::is_exact_touch_end`: an end at an exact touch
/// of the two surfaces is no crossing to refine).
#[test]
fn a_cylinder_crossing_a_ball_beside_its_pole_closes() {
    passes(
        r#"{"cell":{"a":"Y","b":"S","config":"V"},"seed":1,"index":20183,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[0.8257437675475344,-0.13022101943485498,0.548807540449064,-0.3292811310233435],[-0.5421083200707377,0.08548851249779363,0.8359487326032603,-0.501571467242693],[-0.1557748364354346,-0.9877925897313979,-2.2916917654570133e-6,0.9612752185995364],[0.0,0.0,0.0,1.0]],"scene":[[0.41006791400245285,0.05853394084970309,0.9101747544698657,-2.1470702197160643],[0.5316841989687351,0.7954898825983172,-0.29070218308552476,-1.2323324991984252],[-0.7410507529667933,0.603133173065716,0.29508330531399035,-1.035982201217023],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":{"Pole":9.999999999999999e-6}},"aligned":true,"shift":0.0,"contact":[-1.0529596754964736,-0.7463525011576977,-0.6743254517633015],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ffcb5603 (W3 V cells, tangency design §4.8): two cones
/// crossing at a tangent touch with B's rim 3e-4 from it, touching A's
/// cone. The arms' own ends at the crossing were dragged onto the rim's
/// junction (`clip_boundary_edges::is_exact_touch_end`); kept, two arms
/// of the crossing run to the rim's touch and bound a lens
/// 1e-7 deep on both cones, whose samples read their side by noise: A's
/// read Outside beside an Outside neighbour across an arm, and the union
/// kept both sides of the lens (`lens_fragments`: a lens takes the side
/// its neighbours across its records agree on).
#[test]
fn a_cone_crossing_a_cone_at_its_rims_touch_closes() {
    passes(
        r#"{"cell":{"a":"C","b":"C","config":"V"},"seed":1,"index":10090011,"variant_a":2,"variant_b":1,"b_scale":1.0,"place_b":[[-0.10503805547028794,0.9794379949888413,-0.1722446541267949,-0.4798520087985809],[-0.49165643982119356,0.09940964411194615,0.8650963344274935,1.4483270373780548],[0.8644309990307901,0.1755532301889507,0.4711052019293095,-1.0662028081044626],[0.0,0.0,0.0,1.0]],"scene":[[0.1449435143685564,-0.9198997457095992,0.364384186657307,1.0612863103163503],[0.7433686438073941,0.3442877669155191,0.57347100445993,1.6648574699672594],[-0.65298884909809,0.18775087588541425,0.7337269052977491,-2.7515287808416606],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[0.27460465418162716,0.012881397808823802,-2.832513339045187],"curve":[],"curve_closed":false,"third":null}"#,
    );
}
