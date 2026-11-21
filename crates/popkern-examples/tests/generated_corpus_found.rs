//! Generated-corpus scenes the W1.3 baseline found failing (PK-1276), one
//! per cause, replayed through the corpus's oracles. Each failed in the
//! baseline (seed 1, `CORPUS_RANDOM=30`) with an invalid cellular result
//! or a shell gap; the cause is named on each.

#[path = "generated_corpus/mod.rs"]
mod generated_corpus;

use generated_corpus::judge::Op;
use generated_corpus::judge::judge;
use generated_corpus::scene::Scene;
use popkern_bool::NonManifoldPolicy;

fn passes(json: &str) {
    let scene: Scene = serde_json::from_str(json).expect("a scene descriptor");
    let verdict = judge(&scene);
    let failures: Vec<String> = verdict
        .failures
        .iter()
        .map(|f| format!("{} {:?}: {}", f.oracle, f.op, f.detail))
        .collect();
    assert!(failures.is_empty(), "{failures:#?}");
}

/// Failed in the baseline: a record end on a box face's edge was dragged 3.2e-2 onto a torus patch iso-line by junction refinement (`ff_junction_refine`): the box face's walk never split along it.
#[test]
fn a_box_crossing_a_torus_keeps_its_box_edge_ends() {
    passes(
        r#"{"cell":{"a":"P","b":"T","config":"X"},"seed":1,"index":1,"variant_a":1,"variant_b":0,"b_scale":1.0,"place_b":[[0.9922634206316426,0.050326498285916677,0.11349250039835035,-0.2940273551765813],[-0.09340394919952683,0.9048328124455881,0.41539533434517467,1.2024350699051776],[-0.08178634574502117,-0.4227822431135027,0.9025331952662966,-0.4316324823344819],[0.0,0.0,0.0,1.0]],"scene":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline: the plane x cone section's quarter round the far end of a near-parabolic ellipse came out as the complementary arc (`plane_cone_oblique`): the box face read one crossing on a closed record and dropped it.
#[test]
fn a_box_crossing_a_cone_keeps_its_section_simple() {
    passes(
        r#"{"cell":{"a":"P","b":"C","config":"X"},"seed":1,"index":10090000,"variant_a":1,"variant_b":1,"b_scale":1.0,"place_b":[[0.4058223198769544,0.8532802525590422,-0.32744626319819914,1.16890754246454],[-0.6450043696909191,0.013558309923947102,-0.7640585941612241,-0.9677844626901481],[-0.6475164922759165,0.5212763018062336,0.5558717562585003,-0.7785213633782198],[0.0,0.0,0.0,1.0]],"scene":[[0.7477989942502219,0.45336492525703553,0.4850328944979293,2.2187798357957735],[-0.38784093947714626,0.8912399815474512,-0.23509721596986988,-2.3456075099808027],[-0.5388655396885585,-0.012310151826100527,0.842301840374445,0.5746830304802213],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline: an open record straddling the cylinder's moved seam (u = 2.31) found one of its two crossings with the wall's top edge (`clip_lift_crossings`): the wall lost its second branch.
#[test]
fn a_cylinder_crossing_a_sphere_finds_both_branches() {
    passes(
        r#"{"cell":{"a":"Y","b":"S","config":"X"},"seed":1,"index":10105135,"variant_a":0,"variant_b":1,"b_scale":1.0,"place_b":[[0.7757482609897833,0.605821412351067,-0.1766212102446519,-0.2639443582303274],[-0.6125815698632422,0.7901607399785715,0.019743993020715667,-0.19973699293182828],[0.15152047992008458,0.09287852999200004,0.9840808517755604,0.21659452919021493],[0.0,0.0,0.0,1.0]],"scene":[[-0.6979770160893686,0.4624321731963393,-0.5467948154508164,2.9078460294882102],[-0.6544380405654198,-0.7219081739683708,0.22485426262036168,-1.2826162521097486],[-0.29075580144147695,0.5147864348936978,0.8065085184779892,1.7086000915198145],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline: a record split at a seam hit 1.14e-6 beside
/// it left the vertex at the bare tolerance (`table::vertex_reaches`,
/// `vertex_reach`), and the common and A - B shells read open. Its union
/// and B - A still decline on a closed rim edge through that vertex left
/// unsplit (a pinch fan; PK-1285), so only the served operations are
/// judged here.
#[test]
fn a_cylinder_crossing_a_cone_records_its_vertex_reach() {
    let scene: Scene = serde_json::from_str(r#"{"cell":{"a":"Y","b":"C","config":"X"},"seed":1,"index":3027,"variant_a":1,"variant_b":1,"b_scale":1.0,"place_b":[[0.3680121018134219,-0.42001878328992936,0.8295488620946418,0.6642647743792004],[-0.921565765770539,-0.2833821559684715,0.26535088663583023,0.0],[0.12362698847188652,-0.8621361698492347,-0.49136299449496984,-0.09372882121497605],[0.0,0.0,0.0,1.0]],"scene":[[0.9877194112183519,0.09139907254629043,-0.12671453839299657,2.7558329668963015],[-0.03003421363251476,0.9069833040764582,0.4200943133845391,-2.566994170426126],[0.1533242013315208,-0.41112953635591404,0.8985901143579438,0.4097634172040836],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":{"Seam":0.0},"steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#).expect("a scene descriptor");
    let placed = generated_corpus::scene::place(&scene);
    let tol = popkern_core::tolerance::Tolerance::default();
    for op in [Op::Common, Op::AminusB] {
        let body = generated_corpus::judge::run(&placed, op, NonManifoldPolicy::Decline, tol)
            .unwrap_or_else(|e| panic!("{op:?}: {e}"));
        let issues = body.validate_all();
        assert!(issues.is_empty(), "{op:?}: {:?}", issues.first());
    }
}

/// Failed in the baseline: a cylinder cap's line record ended at radius
/// 1.19 inside a cap of radius 1.2, and the one-crossing clip read the
/// end against the cap's eight-sample polygon (its sides 1.11 from the
/// centre) as outside, dropping the record (`clip_uv`, the dense banded
/// verdict): both caps' walks lost their cap-to-cap edge.
#[test]
fn a_cylinder_crossing_a_cylinder_keeps_its_cap_record() {
    passes(
        r#"{"cell":{"a":"Y","b":"Y","config":"X"},"seed":1,"index":4036,"variant_a":2,"variant_b":1,"b_scale":1.0,"place_b":[[0.9949099022367337,-0.04067062031902266,-0.0921964591194204,0.812682215768239],[0.0705153908221635,0.9345891026494796,0.3486700286317108,-3.3861758996938176e-7],[0.07198517964477547,-0.35339653344608013,0.932699857435286,-0.22070143094040295],[0.0,0.0,0.0,1.0]],"scene":[[0.9835671974016756,-0.16877387377238295,0.06411667277141676,-2.2716365499999913],[0.17177274676075227,0.9841300301202586,-0.044521986543342704,2.8454943859692303],[-0.0555849949687833,0.054803862522297284,0.996948767483544,2.817215282319915],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":{"Seam":5e-7},"steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline (W3c, tangency design §4.8): a cylinder lying on another's ruling 5e-7 deep: the glued edge sat on the overlap's midline and the radial order's first probe read the overlap's order (`radial_order::tie_runs`, the lens tie): the general union declined RegionConflict.
#[test]
fn two_cylinders_on_a_ruling_overlapping_within_the_tolerance_glue() {
    passes(
        r#"{"cell":{"a":"Y","b":"Y","config":"K"},"seed":1,"index":10090008,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[1.0,0.0,0.0,1.9999995],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,-0.25079123101451406],[0.0,0.0,0.0,1.0]],"scene":[[-0.5125062128307176,-0.09468918311805455,0.8534467414023856,-2.3397804480211195],[0.01376243083057143,0.9928679358172475,0.11842236918602703,0.23242569757742837],[-0.8585732218672726,0.0724377016920987,-0.507547832293684,-2.283270647499994],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-0.5,"rung":{"offset":-0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-2.749785418285694,0.24343564224188552,-2.970129224993812],"curve":[[-3.603232159688079,0.1250132730558585,-2.462581392700128],[-3.389870474337483,0.15461886535236524,-2.589468350773549],[-3.1765087889868866,0.184224457648872,-2.71635530884697],[-2.96314710363629,0.21383004994537874,-2.843242266920391],[-2.749785418285694,0.24343564224188552,-2.970129224993812],[-2.536423732935097,0.27304123453839224,-3.097016183067233],[-2.323062047584501,0.302646826834899,-3.223903141140654],[-2.1097003622339043,0.3322524191314058,-3.350790099214075],[-1.896338676883308,0.3618580114279125,-3.477677057287496]],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline (W3c, tangency design §4.8): a cylinder lying along a box's top: its one tangent record entered the arrangement as a crossing, the split pieces kept their copies and the general glue declined a wedge in two regions (`contact_curves`).
#[test]
fn a_cylinder_lying_on_a_box_is_a_contact_curve() {
    passes(
        r#"{"cell":{"a":"P","b":"Y","config":"K"},"seed":1,"index":2018,"variant_a":1,"variant_b":0,"b_scale":0.3145764348029479,"place_b":[[6.123233995736766e-17,0.0,1.0,0.0],[0.0,1.0,0.0,0.08273375886683368],[-1.0,0.0,6.123233995736766e-17,1.0516606478423582],[0.0,0.0,0.0,1.0]],"scene":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-0.5,"rung":{"offset":-0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[4.440892098500626e-16,0.08273375886683354,0.7999999999999996],"curve":[[-0.31457643480294717,0.08273375886683354,0.7999999999999998],[-0.23593232610221104,0.08273375886683354,0.7999999999999998],[-0.15728821740147314,0.08273375886683354,0.7999999999999992],[-0.07864410870073701,0.0827337588668331,0.7999999999999996],[4.440892098500626e-16,0.08273375886683354,0.7999999999999996],[0.07864410870073746,0.08273375886683332,0.8],[0.15728821740147447,0.0827337588668336,0.7999999999999992],[0.23593232610221104,0.08273375886683315,0.7999999999999996],[0.3145764348029485,0.08273375886683337,0.7999999999999992]],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline (W3c, tangency design §4.8): a torus lying on a box's top: its closed contact circle's second half joined two vertices of one hole and the imprint declined (`imprint_edge_hole`).
#[test]
fn a_torus_lying_on_a_box_imprints_its_contact_circle() {
    passes(
        r#"{"cell":{"a":"P","b":"T","config":"K"},"seed":1,"index":0,"variant_a":0,"variant_b":0,"b_scale":0.45,"place_b":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.78],[0.0,0.0,0.0,1.0]],"scene":[[-0.2931344185988276,-0.8232243279409148,0.48618300928658165,1.5403781700082835],[-0.5339711906828306,-0.2808505649904419,-0.7974946568255549,-1.8065393026623815],[0.7930617757485777,-0.49338085292274825,-0.357251107508058,1.3861025364359563],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[1.9903805616235992,-1.9966916537889858,0.7434985130268894],"curve":[[1.6737953895368656,-2.573380539726443,1.6000052308353534],[1.5157261670181323,-2.609469129353688,1.4654496986917176],[1.4058201620260427,-2.5961657307609602,1.2861813348012103],[1.3608095675202263,-2.5354956657883934,1.0894921224943959],[1.3875468384921383,-2.4366954018525533,0.9053262113528374],[1.4819614678102702,-2.3148063835110477,0.7617211918796956],[1.6296796840316823,-2.1883851088796686,0.6805396265070126],[1.8082127278056292,-2.0766780710093355,0.6741406726500572],[1.9903805616235992,-1.9966916537889858,0.7434985130268894],[2.1484497841423327,-1.960603064161741,0.8780540451705254],[2.258355789134422,-1.973906462754469,1.0573224090610327],[2.3033663836402383,-2.0345765277270353,1.254011621367847],[2.2766291126683265,-2.133376791662876,1.4381775325094055],[2.1822144833501946,-2.2552658100043814,1.5817825519825475],[2.0344962671287825,-2.38168708463576,1.6629641173552303],[1.855963223354836,-2.4933941225060936,1.6693630712121859]],"curve_closed":true,"third":null}"#,
    );
}

/// Failed in the baseline (W3c, tangency design §4.8): a ball 1e-5 wider than a coaxial cylinder: consolidation's loop dedupe dropped the second of two exact circles 8e-3 apart and B - A served inside out (`intersect_consolidate`, the lateral band).
#[test]
fn a_ball_just_wider_than_a_coaxial_cylinder_keeps_both_circles() {
    passes(
        r#"{"cell":{"a":"Y","b":"S","config":"K"},"seed":1,"index":10090000,"variant_a":0,"variant_b":0,"b_scale":0.8889,"place_b":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,-0.02726275417719121],[0.0,0.0,0.0,1.0]],"scene":[[0.5707106841991072,0.23883114457626994,0.7856519581350102,-0.6292234664017644],[0.20656052395198884,0.8842484134044848,-0.41885259141651643,-2.389757514459214],[-0.7947465412857527,0.4013283291421629,0.455316930656321,2.257379603663707],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-10.0,"rung":{"offset":-10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-1.1072110499645138,-2.5435868583845367,2.8807636431391126],"curve":[[-0.19407395524594234,-2.2130900200613546,1.609169177081908],[-0.15571078878186984,-1.9549590326556954,1.7804317216892511],[-0.19269648918024557,-1.7612837219560222,2.022415461591294],[-0.29940031883162566,-1.6615493983452971,2.298280565976016],[-0.45957758694421214,-1.670939708499358,2.566029073424241],[-0.6488427564584147,-1.7880250628206227,2.7848987008232404],[-0.8383819210087302,-1.9949802774950582,2.9215685314596307],[-0.9993394609691784,-2.260298297142175,2.9552318225341208],[-1.1072110499645138,-2.5435868583845367,2.8807636431391126],[-1.1455742164285865,-2.801717845790196,2.70950109853177],[-1.1085885160302107,-2.9953931564898695,2.4675173586297263],[-1.0018846863788307,-3.0951274801005946,2.1916522542450054],[-0.8417074182662442,-3.0857371699465337,1.9239037467967801],[-0.6524422487520413,-2.968651815625269,1.70503411939778],[-0.462903084201726,-2.7616966009508332,1.5683642887613902],[-0.3019455442412781,-2.4963785813037167,1.5347009976869002]],"curve_closed":true,"third":null}"#,
    );
}

/// Failed in the baseline (W3c, tangency design §4.8): a cone lying on a box along its seam: glued, its apex kept one edge in a ring whose every face holds two fins, and the validator read LowValence (`validation::manifold`).
#[test]
fn a_cone_lying_on_its_seam_glues_into_its_apex() {
    passes(
        r#"{"cell":{"a":"P","b":"C","config":"K"},"seed":1,"index":0,"variant_a":0,"variant_b":0,"b_scale":0.25660011963983365,"place_b":[[-0.4472135954999579,0.0,0.8944271909999159,0.0],[0.0,1.0,0.0,0.04798145298945416],[-0.8944271909999159,0.0,-0.4472135954999579,0.7032795558988644],[0.0,0.0,0.0,1.0]],"scene":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-0.0516397779494322,0.04798145298945433,0.6000000000000003],"curve":[[-0.3098386676965932,0.04798145298945433,0.6000000000000001],[-0.24528894525980283,0.04798145298945433,0.6000000000000001],[-0.1807392228230127,0.04798145298945411,0.6000000000000001],[-0.11618950038622256,0.04798145298945411,0.5999999999999999],[-0.0516397779494322,0.04798145298945433,0.6000000000000003],[0.012909944487357938,0.04798145298945411,0.5999999999999999],[0.07745966692414852,0.04798145298945422,0.6000000000000001],[0.14200938936093865,0.04798145298945422,0.5999999999999999],[0.2065591117977288,0.04798145298945411,0.5999999999999999]],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 12320388 (W3c V cells, tangency design §4.8): a cylinder and a frustum crossing at their tangency: the marcher ran through the crossing and a second branch's march ran back over the first's approach (records parallel within 1e-9); the crossing is now a node every trace ends at, its branch halves seeded from it (`marching_crossing`).
#[test]
fn a_cylinder_and_a_frustum_crossing_at_their_tangency_meet_at_one_node() {
    passes(
        r#"{"cell":{"a":"Y","b":"F","config":"V"},"seed":1,"index":6,"variant_a":2,"variant_b":2,"b_scale":1.0,"place_b":[[0.5205114357408136,-0.850657137908059,-0.07382600482964147,-0.12885194712033887],[0.7553081240407371,0.4183891207174323,0.5044404637054373,0.18183621371549663],[-0.3982178838539269,-0.31832841122237787,0.8602845689581657,0.3003591690226709],[0.0,0.0,0.0,1.0]],"scene":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-1.0578866653597736,0.5664590040364603,0.09155463912512185],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 12320388 (W3c V cells, tangency design §4.8): a cylinder and a NURBS solid crossing at their tangency, the marcher's Newton stalling at the crossing (`newton-exhausted`); the trace ends there and the branch halves leave it from solved departure seeds.
#[test]
fn a_cylinder_and_a_nurbs_solid_crossing_at_their_tangency_meet_at_one_node() {
    passes(
        r#"{"cell":{"a":"Y","b":"W","config":"V"},"seed":1,"index":10104127,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[0.7716510519710897,0.07095019362458457,-0.6320765175329206,-0.1771487009536603],[-0.5986036157509438,-0.254918766219136,-0.7594011679206255,0.6539177212648126],[-0.2150078259082842,0.9643559989204958,-0.15423729167826278,-1.1946969286554787],[0.0,0.0,0.0,1.0]],"scene":[[0.22542527453606853,0.9529585138918373,0.20261667355230029,-0.24964997898114216],[0.9722726875680243,-0.20676832959685942,-0.10923680188002415,-0.8936620824994708],[-0.062203429243003525,0.22162339379403456,-0.9731463429071832,2.2456264058280215],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-1.2063810563367723,-0.9212212600247144,2.9646557587978015],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 2f33ba3d (W3c V cells, tangency design §4.8): a cylinder and a cone crossing at their tangency, a branch leaving the crossing and returning to it: marched as a closed loop, its seam split moved the node off the record's ends (the loop-successor class); one branch half re-traced eight times, each copy within its band `tol / sin θ` of the others.
#[test]
fn a_cylinder_and_a_cone_crossing_at_their_tangency_end_their_loop_at_the_node() {
    passes(
        r#"{"cell":{"a":"Y","b":"C","config":"V"},"seed":1,"index":11100,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[-0.6774916914684359,0.5839565779837337,0.4472132858276518,0.41862172380809504],[-0.23647018818503854,-0.7486614474091546,0.6193447240939098,1.2040840244241662],[0.6964817715424935,0.3138482948699207,0.6453002322301404,1.515499548803094],[0.0,0.0,0.0,1.0]],"scene":[[-0.062116520136444156,-0.7217111713421056,-0.6894015688161323,-2.229024918479282],[0.2340589063829086,0.6609395513598968,-0.7130044444398753,1.9339504045796758],[0.9702360363689733,-0.20564993218618915,0.12786766253764179,-2.02505401404937],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":{"Seam":5e-7},"steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-2.471044253485961,1.9819479103887485,-1.0214504362758043],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 2f33ba3d (W3 near-miss cells, tangency design §4.8): two cylinders overlapping by ten tolerances at an external touch, a thin loop through θ = 0 of the scan: its entry bracket ran across the whole turn and the fitted loop strayed 2 units off both cylinders; the boolean kept nothing and their common was served empty.
#[test]
fn two_cylinders_overlapping_by_ten_tolerances_have_a_common() {
    passes(
        r#"{"cell":{"a":"Y","b":"Y","config":"E"},"seed":1,"index":6065,"variant_a":1,"variant_b":0,"b_scale":1.0,"place_b":[[-0.6029876864903324,-0.6076457476229745,0.5168873139638396,-1.653862486644372],[-0.7923039519824704,0.3805756388297473,-0.4768821980345765,-0.9231837016683294],[0.09306072003802729,-0.6970859549044479,-0.7109225512395551,0.7590800521170626],[0.0,0.0,0.0,1.0]],"scene":[[0.6992644599753052,-0.34510889250743326,0.6260423846096506,2.2859129708763524],[0.13307368087064597,0.9232851934599682,0.36032741638283183,2.0264513398141304],[-0.7023678597964049,-0.16865439173019492,0.691545432835136,2.692388266086395],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-10.0,"rung":{"offset":-10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[2.1675132484196844,1.3163980527123966,3.405419472128168],"curve":[],"curve_closed":false,"third":null}"#,
    );
}
