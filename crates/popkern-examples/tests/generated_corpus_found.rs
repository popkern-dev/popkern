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

/// Failed at ec3e6ea3 (W3 near-miss cells, tangency design §4.8): two cones overlapping by ten tolerances at an external touch met in a loop no starting-point grid held, the pair emitted nothing and their common was served empty; the overlap's deepest point, a near miss of the contact search, now seeds it.
#[test]
fn two_cones_overlapping_by_ten_tolerances_have_a_common() {
    passes(
        r#"{"cell":{"a":"C","b":"C","config":"E"},"seed":1,"index":10097064,"variant_a":1,"variant_b":1,"b_scale":1.0,"place_b":[[-0.8586182736176788,0.45592653857769627,0.23431955025212403,-1.7620309767370106],[0.31211903518213674,0.10237088949307838,0.9445114657119672,1.1031531079760035],[0.4066403424240756,0.8841103960507175,-0.23020086774005494,-0.9180693876928406],[0.0,0.0,0.0,1.0]],"scene":[[0.9067201038272686,0.3244874599408267,-0.2693817767715881,0.6681044108112233],[-0.33871817240709107,0.9408635098093702,-0.006771675597170644,-0.8794986985019548],[0.2512541601579315,0.09738451750839841,0.9630096587018164,-0.30624798324265257],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-10.0,"rung":{"offset":-10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-0.07604110927878394,-0.673664965922524,-1.6464869897562813],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 62c9d0b0 (W3 V cells, tangency design §4.8): a cylinder and a sphere crossing tangentially, the pinch opened half a tolerance, were split at the ruling scan's dip there by the crossing band, a gap of `ε² / 8R`; the two loops passed within the tolerance and both differences were served where the touch declines. A dip splits only past the tolerance.
#[test]
fn a_cylinder_and_sphere_opened_within_the_tolerance_keep_their_crossing() {
    passes(
        r#"{"cell":{"a":"Y","b":"S","config":"V"},"seed":1,"index":10119264,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[-0.7920388210011786,0.053483737753694305,0.6081233393177382,-0.2811922043289632],[-0.445616986255699,-0.7315180123093211,-0.5160493185998876,-0.2844843479435528],[0.4172529299926496,-0.6797211836208388,0.6032239260420922,0.7554964747568667],[0.0,0.0,0.0,1.0]],"scene":[[-0.9560246956801274,0.2431225610310686,0.16403719568260833,0.03340049570317616],[-0.06253135511745594,-0.7154169839151645,0.6958939349878716,-0.816375801926835],[0.28654251146539256,0.6550343192867253,0.6991590875327888,-2.0425493251193556],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.5,"rung":{"offset":0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-0.3418246781010671,-0.843400355681363,-0.8470366042301145],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 62c9d0b0 (W3 near-miss cells, tangency design §4.8): a cylinder overlapping a NURBS torus by two tolerances meets it in a loop 0.018 round, which the overlap rescue seeds; the march took 15 steps a lap, tested closure only after 20, closed on its second lap, and the record covered itself: every assembly declined a loop successor. Closure is a return into the start's ball after leaving it.
#[test]
fn a_cylinder_overlapping_a_torus_by_two_tolerances_traces_its_loop_once() {
    passes(
        r#"{"cell":{"a":"Y","b":"W","config":"I"},"seed":1,"index":10101099,"variant_a":2,"variant_b":2,"b_scale":0.3159427994716938,"place_b":[[-0.09453316723024041,0.7250784197763441,0.6821471728799139,0.23209139101241943],[0.029641106296108055,0.6869620097519152,-0.7260885634515616,-0.7785163680334312],[-0.9950803410333069,-0.048419854731808355,-0.08643282106230454,1.068867147291821],[0.0,0.0,0.0,1.0]],"scene":[[0.8960142112281784,-0.4367778277786548,0.07989782498982734,-1.4782381600317276],[0.438400165884543,0.8416716900423935,-0.3152682361634076,-2.388092460666869],[0.07045443796916176,0.31751203968060926,0.9456332676192758,0.5748002249590707],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-2.0,"rung":{"offset":-2.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-0.6602706139323506,-3.6320235271763677,1.2932892139309258],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 9fae7006 (W3 V cells, tangency design §4.8): two cylinders crossing tangentially with the pinch closed half a tolerance read as two loops (the ruling scan's crossing band was the roots' separation, `(a·tol)²`, a gap of `tol² / 8R`), and both differences were served where the pinch declines. A minimum of the discriminant is a crossing where its double root lies within the reach of the other cylinder.
#[test]
fn two_cylinders_crossing_within_half_a_tolerance_pinch_their_differences() {
    passes(
        r#"{"cell":{"a":"Y","b":"Y","config":"V"},"seed":1,"index":10091009,"variant_a":1,"variant_b":2,"b_scale":1.0,"place_b":[[-0.43852979780804974,-0.8985904469756895,0.015060711751501442,0.18190890973339247],[0.8842535671457894,-0.43440627958195677,-0.1714258243372616,0.2146815140078116],[0.16058407587422194,-0.061857843995677925,0.9850819060930053,-0.528683837681446],[0.0,0.0,0.0,1.0]],"scene":[[0.8630998770542098,0.3365278211711584,-0.3765735357228397,-1.7468327019589616],[-0.4522808092211521,0.8468372362827001,-0.27983703267319354,-1.6517694895800408],[0.2247235453602221,0.4118442919634446,0.883110189805248,1.0686814413185584],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-0.5,"rung":{"offset":-0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-2.86366656888624,-1.4444658293734944,1.3425100758115223],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 9fae7006 (W3 V cells, tangency design §4.8): two cylinders crossing tangentially with the pinch opened half a tolerance read as two loops turning apart near it, and an edge's curve and trims disagreed by 6.1e-4 there; the crossing within the reach is one node.
#[test]
fn two_cylinders_opened_within_half_a_tolerance_meet_at_one_crossing() {
    passes(
        r#"{"cell":{"a":"Y","b":"Y","config":"V"},"seed":1,"index":10119261,"variant_a":2,"variant_b":0,"b_scale":1.0,"place_b":[[0.3298781954401002,0.8625256602204462,0.38370543602412843,-0.5583271378744863],[-0.8329115150465449,0.07461002589598204,0.5483536743185642,-0.10021564316988363],[0.44434084245773775,-0.5004825965976083,0.7430197751250389,-0.23821094599030124],[0.0,0.0,0.0,1.0]],"scene":[[0.7868476372984512,0.010415507712960986,-0.6170594078992125,2.7739678057374153],[0.3255857364918281,-0.8563990719473314,0.40071755359707323,-2.647483226673994],[-0.5242754274910463,-0.5162094020519489,-0.6772467270960045,-1.7763399432279208],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.5,"rung":{"offset":0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[1.8789603298424291,-3.4733094428522207,-1.7571014614401688],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 1da6d059 (W3 X3 cells, design §5): a ball inside the union of two crossing balls, crossing one of them, read N-ary (A ∪ B) − C as A ∪ B whole, validator-clean (12.53 against 11.87, a Monte Carlo 11.85): the shell partition nested its void by vertex boxes (the union's only vertices its ring) and by parity on faces with no trims yet, two merges joined its caps into one ball face bounded by an arc out and back, and the slit sweep retired it. Now nested by face boxes on a trimmed copy, then after the output repairs by parity, and no merge leaves a slit.
#[test]
fn a_ball_inside_two_crossing_balls_cuts_a_void_from_their_union() {
    passes(
        r#"{"cell":{"a":"S","b":"S","config":"X3"},"seed":1,"index":10090,"variant_a":2,"variant_b":0,"b_scale":1.0,"place_b":[[0.3009187793811897,-0.5073458836280031,-0.8074949179911017,-2.301488186386046e-6],[-0.9143246071125679,0.08715969004441138,-0.3954917208486302,7.132133082419954e-6],[0.2710321033460495,0.8573233595307452,-0.43765083817894523,1.0491986158184645],[0.0,0.0,0.0,1.0]],"scene":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":{"Pole":9.999999999999999e-6},"steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":{"variant":0,"scale":0.6,"place":[[-0.9459415701911844,0.039535348750924626,-0.32191847101926174,-0.3376881539653659],[-0.30946718216171354,0.1870675710020001,0.9323281541626319,0.5479496854139418],[0.0970804251593765,0.981551160194187,-0.16472009887175343,0.4201300795358005],[0.0,0.0,0.0,1.0]]}}"#,
    );
}

/// Failed at 1da6d059 (W3 X3 cells, design §5): as the first void scene, a second placement (N-ary 19.62 against the fold's 18.96): the void kept by a merged ball face no trim reader sees.
#[test]
fn a_second_ball_inside_two_crossing_balls_cuts_a_void_from_their_union() {
    passes(
        r#"{"cell":{"a":"S","b":"S","config":"X3"},"seed":1,"index":10098072,"variant_a":2,"variant_b":2,"b_scale":1.0,"place_b":[[-0.4381647209784447,-0.79129007254789,-0.4264641818219165,0.665237623525996],[-0.7453092868148852,0.5850369984255396,-0.3197589364829731,-1.2814197169322021],[0.5025193969165456,0.17774063002136004,-0.8460984128115602,0.11064669958850047],[0.0,0.0,0.0,1.0]],"scene":[[0.7597138673516144,0.31437485429080264,-0.5692128694463,-1.4791763848166941],[0.3985596528642457,0.46653821972207843,0.789615281417024,0.7276742016493403],[0.5137947477956405,-0.826746962817639,0.22913842237700727,-1.3678960497095731],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":{"variant":0,"scale":0.6,"place":[[0.11954461327663143,0.46687580725990674,0.8762054930391335,-0.25479190239023364],[0.2756520095194005,0.832228977323782,-0.48105186721447,-0.38668042286223825],[-0.9537950802370795,0.2990349643403296,-0.02920676321540764,-0.4488059499624154],[0.0,0.0,0.0,1.0]]}}"#,
    );
}

/// Failed at aef1116a (W3 X3 cells, design §5): the fold (A ∩ B) ∩ C of three balls rejected the hit of C's seam meridian on the lens's cap as off the face (the cap's trim box spanned only its ring's three vertices), C's seam stayed unsplit at the record crossing it, its chart walk declined, and the common was one record edge with a valence-1 vertex. A curved face's trim box takes points along its edges.
#[test]
fn three_crossing_balls_share_their_common_through_a_seam() {
    passes(
        r#"{"cell":{"a":"S","b":"S","config":"X3"},"seed":1,"index":10108162,"variant_a":2,"variant_b":1,"b_scale":1.0,"place_b":[[-0.31505366627238607,-0.5196703342156744,0.7941561125524991,0.7794046640100387],[-0.33935531235520694,-0.719766128536323,-0.6056192633892463,0.5433309326561414],[0.8863290355983088,-0.46030366493001407,0.05041206906455542,-1.167797700460516],[0.0,0.0,0.0,1.0]],"scene":[[-0.17696883000700786,0.2114876199291988,0.9612257902400629,2.3658543671525],[-0.6162957080831744,0.7376548388171936,-0.2757624683814658,-2.6590633068689007],[-0.7673732034702603,-0.6412006904131209,-0.00020300129939590517,-1.2526768350854007],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":{"variant":0,"scale":0.6,"place":[[0.9957525969343901,0.08418964054077151,-0.03726754786719661,-0.08064117212719424],[-0.07968068877965083,0.9908092064897399,0.10930738388037807,0.1699392156501659],[0.04612757888745662,-0.10587360747980391,0.9933091289750607,-0.5019149890057615],[0.0,0.0,0.0,1.0]]}}"#,
    );
}

/// Failed at 9fae7006 and cb2f9807 (W3 E cells, tangency design §4.8): a
/// cylinder touching a cone at millimetre scale meets it in a loop 0.0006
/// round, ten steps a lap under the near-tangent step floor; the march
/// asked for twice that floor's ball before closing, never left it, and
/// orbited, and every assembly declined a loop successor.
#[test]
fn a_cylinder_touching_a_cone_at_millimetre_scale_closes_its_small_loop() {
    passes(
        r#"{"cell":{"a":"Y","b":"C","config":"E"},"seed":1,"index":10109173,"variant_a":1,"variant_b":1,"b_scale":1.0,"place_b":[[-0.7086146125048254,-0.6357474628417116,-0.30608903024277545,-1.499473167601841],[-0.6266301406507099,0.3676059242252121,0.6871685028452619,-0.8632706419898311],[-0.32434549137101826,0.6787422544418862,-0.6588694516081608,-0.6941410638184893],[0.0,0.0,0.0,1.0]],"scene":[[0.0008577590420690879,0.000513951390654189,0.000010168273862652906,0.0013218835105378154],[-0.0005072341963484659,0.0008494280603984229,-0.00014555219085434325,-0.0007194902396543323],[-0.00008344396804710688,0.0001196910115773047,0.0009892983199946087,0.00194726862842652],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":-2.0,"rung":{"offset":-2.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[0.00045475071509627144,-0.0012472428458575376,0.0021453778636017938],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at cb2f9807 (W3 V cells, scale; tangency design §4.8): two
/// millimetre cylinders crossing tangentially, overlapping by two
/// tolerances. The clip's uv boundary band (a thousandth of the 2π of
/// u) was wider than the face's height in v, and residues within 400
/// tolerances of any record were dropped as covered.
#[test]
fn two_millimetre_cylinders_crossing_tangentially_clip_in_space() {
    passes(
        r#"{"cell":{"a":"Y","b":"Y","config":"V"},"seed":1,"index":10093027,"variant_a":0,"variant_b":1,"b_scale":1.0,"place_b":[[-0.4239798747602541,0.6018288210599926,0.6767888414711184,0.3117753991173719],[0.5705842468748785,0.75782456255649,-0.3164420161803619,0.07279914771791562],[-0.7033311332626382,0.2520000050150875,-0.6646964077350568,-0.11534992859605891],[0.0,0.0,0.0,1.0]],"scene":[[0.0008225129071654389,-0.0004999997999831893,-0.00027105851317202486,-0.000336461028824381],[0.00048114745599601717,0.0008658462602346747,-0.0001371399986371176,0.0021308187979465877],[0.00030326497182302264,-0.000017619695071102788,0.0009527433564243745,0.0019098085426906863],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":-2.0,"rung":{"offset":-2.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-0.0002822217848290445,0.0013254311188369706,0.0019231814839848979],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 1a7f4c16 (W3 scale; tangency design §4.8): two balls a
/// thousand across crossing in a circle; the degree-1 sphere pcurve
/// reached its 4097-sample cap 3.4e-4 off the exact circle and every
/// boolean declined. The cubic chart fit holds it.
#[test]
fn two_balls_a_thousand_across_keep_their_circle_pcurves_on_it() {
    passes(
        r#"{"cell":{"a":"S","b":"S","config":"X"},"seed":1,"index":10109171,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[0.6270455205369367,0.5347522557868213,0.5664397056222913,0.608112211941213],[0.7775571631809831,-0.38569371614098746,-0.4966338845822922,-0.3198232634796464],[-0.04710385504928405,0.7518513032909244,-0.6576479640196832,0.7906882117273409],[0.0,0.0,0.0,1.0]],"scene":[[620.6267901632712,484.0887167101862,616.8310155022365,2392.279206535422],[-728.2419017258311,64.24731495839767,682.3019969862183,-1954.0848387976412],[290.6649616007753,-872.6570900843709,392.4072925202822,2920.053110608401],[0.0,0.0,0.0,1.0]],"scale":1000.0,"offset":0.0,"rung":{"offset":0.0,"scale":1000.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 1a7f4c16 (W3 scale; tangency design §4.8): a ball and a
/// NURBS torus scaled a thousandfold apart; the torus's edges, 9e-7 off
/// its surfaces at unit size, sat 9e-4 off and the union failed
/// validation. The scale records the imprecision it carries.
#[test]
fn a_nurbs_torus_scaled_a_thousandfold_states_its_imprecision() {
    passes(
        r#"{"cell":{"a":"S","b":"W","config":"D"},"seed":1,"index":10095045,"variant_a":2,"variant_b":0,"b_scale":1.0,"place_b":[[1.0,0.0,0.0,-5.540054093434559],[0.0,1.0,0.0,-10.719288708316363],[0.0,0.0,1.0,2.6408630490731277],[0.0,0.0,0.0,1.0]],"scene":[[889.8115053996416,452.7433613328502,-57.08707058043817,-2467.15527083425],[-415.9657204531726,753.2947118438515,-509.43065918332985,-1315.103245031037],[-187.63796062176957,477.04352618720355,858.6160200325768,2135.69367784366],[0.0,0.0,0.0,1.0]],"scale":1000.0,"offset":0.0,"rung":{"offset":0.0,"scale":1000.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 4412919e (W3 X3 cells, scale; tangency design §4.8): three
/// millimetre bodies, a cylinder crossing a ball in two loops, one
/// winding round its pole far off, the other 18 tolerances beside it.
/// The ruling scan declined for the first, left to the marching path,
/// which stalled at the pole and lost the second.
#[test]
fn a_cylinder_through_a_millimetre_ball_keeps_its_loop_beside_the_pole() {
    passes(
        r#"{"cell":{"a":"Y","b":"S","config":"X3"},"seed":1,"index":10115225,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[0.8263006198915203,-0.5620588984327377,0.03629160039825638,-0.1598919158016851],[0.5507629727185859,0.7928506325101672,-0.2608601587257318,0.4667096266630135],[0.11784495512781021,0.23553698058013126,0.9646941988682836,0.0007895526368704164],[0.0,0.0,0.0,1.0]],"scene":[[0.0009581145618857209,-0.00016111759101063635,0.0002367648795101644,0.0027209520316404563],[0.00012229992053139577,0.0009777482486602463,0.000170443807983263,0.0009435125769824135],[-0.00025895794203026026,-0.00013434836846329943,0.0009564995034764534,0.000988084499735152],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.0,"rung":{"offset":0.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":{"variant":1,"scale":0.6,"place":[[0.9999901839294708,-0.001283667362401549,0.004240783277405582,0.0857338838260544],[0.0012775212912498686,0.9999981302972736,0.001451668456439928,-0.6267138375234889],[-0.004242638807820086,-0.00144623651583151,0.9999899541574838,-0.838313000406333],[0.0,0.0,0.0,1.0]]}}"#,
    );
}

/// Failed at e157a3c5 (W3 X cells; PK-1293): a millimetre ball poking out
/// through a cylinder's cap, minus the cylinder, is two solids. The
/// mid-assembly shell partition found no interior point for the nested
/// cap piece, declined, and served both as one shell; every component is
/// now its own solid and the post-repair cavity nesting alone nests.
#[test]
fn a_ball_through_a_cylinder_cap_less_the_cylinder_is_two_solids() {
    passes(
        r#"{"cell":{"a":"Y","b":"S","config":"X"},"seed":1,"index":10090000,"variant_a":1,"variant_b":1,"b_scale":1.0,"place_b":[[-0.3790964582313595,-0.4904949947156083,0.784665875080199,-0.09704597115379465],[0.8962418992036536,-0.4056556258107373,0.17942678551411195,0.05429041130689899],[0.23029618639540506,0.771270493022412,0.5933847766202791,-0.252440629771427],[0.0,0.0,0.0,1.0]],"scene":[[0.0009454821713099151,-0.0000934711357344402,-0.000311972130998271,0.001898227621617509],[0.00013461767319900095,0.0009844297024726083,0.00011303204391755003,-0.0024896809632662677],[0.0002965493985790149,-0.00014886674468868403,0.0009433435993996136,-0.0003353912604804439],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":0.0,"rung":{"offset":0.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at e157a3c5 (W3 K cells; tangency design §4.8): a box face
/// cuts a cone frustum two tolerances deep along a ruling in a long thin
/// loop. Its march, rounding the loop's tip, passed its start 8.6e-4 off
/// on the far branch and closed there. A march closes only arriving the
/// start's way.
#[test]
fn a_box_overlapping_a_frustum_along_a_ruling_marches_its_thin_loop_whole() {
    passes(
        r#"{"cell":{"a":"P","b":"F","config":"K"},"seed":1,"index":10108162,"variant_a":1,"variant_b":0,"b_scale":0.29234229107591525,"place_b":[[-0.29827499313594674,0.0,0.9544799780350298,0.0],[0.0,1.0,0.0,0.04489383027410398],[-0.9544799780350298,0.0,-0.29827499313594674,1.0092741476736373],[0.0,0.0,0.0,1.0]],"scene":[[0.9615849751514804,-0.18732219096602867,0.2006607393951627,1.4459786484013355],[0.035843269238013675,0.8104058845599191,0.5847713761129099,1.9193736776992338],[-0.27215729939367417,-0.5551150322611897,0.7859883620922133,2.5948823753831425],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-2.0,"rung":{"offset":-2.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[1.535211129500741,2.421228896165955,3.216550584759216],"curve":[[1.2995963770331778,2.4124463098866586,3.2832366045972226],[1.3585000651500685,2.4146419564564825,3.266565099637721],[1.4174037532669594,2.4168376030263063,3.2498935946782193],[1.4763074413838502,2.4190332495961306,3.2332220897187174],[1.535211129500741,2.421228896165955,3.216550584759216],[1.5941148176176319,2.423424542735779,3.199879079799714],[1.6530185057345226,2.4256201893056026,3.183207574840212],[1.7119221938514133,2.427815835875427,3.166536069880711],[1.770825881968304,2.4300114824452512,3.149864564921209]],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at e157a3c5 (W3 E cells, scale; tangency design §4.8): two
/// millimetre cylinders overlapping ten tolerances by one's end. The
/// loop leaves that cylinder 70 tolerances past its rim and closes across
/// its cap in a lens 7.8 tolerances deep, which the cap's clip read
/// outside (its tuned polygon's sag, the exact trim abstaining in its
/// band). The lens takes its side of the rim against its neighbour.
#[test]
fn two_millimetre_cylinders_overlapping_by_an_end_close_their_loop_across_the_cap() {
    passes(
        r#"{"cell":{"a":"Y","b":"Y","config":"E"},"seed":1,"index":10115229,"variant_a":0,"variant_b":0,"b_scale":1.0,"place_b":[[-0.29861548299557766,-0.1975279663875598,0.93370846403474,-1.267256637849513],[-0.10614277797987613,-0.9654036940393652,-0.23817938243656372,-1.3176519293883506],[0.9484526893830988,-0.17023046152185342,0.2673183233001886,0.2576552769766601],[0.0,0.0,0.0,1.0]],"scene":[[0.0009054186755778824,-0.00023303350570886484,0.0003548413830598542,0.0007961741928522246],[0.00014976777942099512,-0.0006067765106883986,-0.0007806355605044623,0.002082350788421522],[0.00039722365760626327,0.000759945821288497,-0.0005144858739985522,0.0019071119665181229],[0.0,0.0,0.0,1.0]],"scale":0.001,"offset":-10.0,"rung":{"offset":-10.0,"scale":0.001,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[0.0009780211216122636,0.0021265705607053658,0.0009781393700735227],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at cf0e39b1 (W3 V cells; tangency design §4.8): a cylinder
/// crossing a ball tangentially, opened by ten tolerances. The crossing's
/// branch turns within 3e-3, a twentieth of the sample spacing, and the
/// projected pcurve's refinement stalled before reaching the turn: the
/// edge's curve and trims disagreed by 6.5e-4. Refinement now waits
/// while a straying interval turns sharply.
#[test]
fn a_cylinder_crossing_a_ball_tangentially_keeps_curve_and_trims_together() {
    passes(
        r#"{"cell":{"a":"Y","b":"S","config":"V"},"seed":1,"index":10094038,"variant_a":1,"variant_b":1,"b_scale":1.0,"place_b":[[0.4289855090521693,0.7895695193048712,0.43880679940940515,-0.027414691982398653],[-0.7195191374577363,0.5923723978409918,-0.36247365850250296,0.1474836762611722],[-0.5461351882595797,-0.16023394291094673,0.8222295541294324,-0.7624063142366733],[0.0,0.0,0.0,1.0]],"scene":[[0.5886599444829197,0.33446426816981223,0.7359436955902251,-2.4530460272721872],[-0.6735579643726737,-0.3004605036439116,0.6753096729502599,1.6199059731679082],[0.44698896898190454,-0.8932284920823332,0.04841200822951974,2.103989529923137],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":10.0,"rung":{"offset":10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-3.235386656127694,1.2773516975274484,3.0269538286207402],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 877ffe1f (W3 V cells; tangency design §4.8): two cones
/// crossing tangentially, opened by ten tolerances. The projected
/// pcurve's local refinement worsened once as it split the crossing's
/// overshoot (4.34e-3 to 4.45e-3) and the stall rule stopped it there;
/// it now waits a round and converges to 8e-7.
#[test]
fn two_cones_crossing_tangentially_refine_their_pcurve_past_a_worse_round() {
    passes(
        r#"{"cell":{"a":"C","b":"F","config":"V"},"seed":1,"index":10113207,"variant_a":1,"variant_b":0,"b_scale":1.0,"place_b":[[0.45650410872202984,0.8770992715893183,0.14933474645035644,0.40517823760558547],[-0.8504698077916277,0.4794778977308226,-0.21633781829005894,-0.06297950899112716],[-0.2613524531256154,-0.028245590188874348,0.9648300792782659,0.30943236481905806],[0.0,0.0,0.0,1.0]],"scene":[[0.8712705624921834,0.47307263285267864,0.1307283097893841,1.8513403634388954],[-0.3835731762923137,0.4901448234576508,0.7827066311631032,-1.8640863881920677],[0.30620128243294487,-0.7320931198170784,0.6085067284368488,2.012415620374168],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-10.0,"rung":{"offset":-10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[1.5898510433650412,-1.4259706956654234,2.137672157307713],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at 877ffe1f (W3 V cells; tangency design §4.8): a cylinder
/// crossing a filleted box's round tangentially, opened by ten
/// tolerances. The round's projected pcurve was served 1.76e-3 off its
/// edge, the curve itself on both surfaces: local refinement now runs
/// wherever the gap is interpolation error rather than misfit-sized, and
/// waits out a round that worsens.
#[test]
fn a_cylinder_crossing_a_round_tangentially_refines_its_pcurve_locally() {
    passes(
        r#"{"cell":{"a":"Y","b":"L","config":"V"},"seed":1,"index":10098072,"variant_a":2,"variant_b":0,"b_scale":1.0,"place_b":[[-0.22715615569328007,0.7873532272409386,-0.573127364975648,0.28185326536956234],[0.7580724011919119,0.5123931762972456,0.40345937519798536,0.03833494747932231],[0.6113315920854894,-0.34282375710749563,-0.7132640156918906,-0.7502189020060916],[0.0,0.0,0.0,1.0]],"scene":[[0.9296426731462757,-0.3160617374835453,0.18939239257248244,0.7476014960129445],[0.3289169522034055,0.9435133309888379,-0.0399528822435901,-0.5786155343052997],[-0.16606666980062398,0.09943627278428134,0.9810883185708107,2.7990635241584707],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":10.0,"rung":{"offset":10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[1.6009983299743071,0.14690408065063687,2.091441095800158],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline (W3 V cells, tangency design §4.8): the cylinder's seam pierced the cone 7.5e-3 rad off tangency, and a record 6.9e-9 off both surfaces passed the seam hit 1.06e-6 away, past its linear deviation but within it read laterally (`tolerance_zone::incidence`, `lateral_deviation`): the record was not split there and every operation's loop broke at the hit. B - A still declines on two triple points of one face triple left apart within the tolerance (PK-1292), so only the other operations are judged here.
#[test]
fn a_cylinder_seam_piercing_a_cone_shallowly_splits_the_record_it_passes() {
    let scene: Scene = serde_json::from_str(r#"{"cell":{"a":"Y","b":"C","config":"V"},"seed":1,"index":4039,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[-0.8640105997777957,0.2407390928134664,0.4421881643180461,0.2089588974777007],[-0.10927420586284498,-0.9470053825072731,0.3020595196899092,0.5755486504238931],[0.49147210643602,0.21266286627806746,0.8445292617202427,2.2148542596392646],[0.0,0.0,0.0,1.0]],"scene":[[-0.9555387502212048,0.13487484935135885,0.2622107393646307,0.5375410677548222],[0.26170775970303906,-0.021765974764099916,0.9649017000989196,-1.0850906767670334],[0.13584824377559843,0.9906235497680324,-0.014499562340594241,-2.9047898602897186],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-2.0,"rung":{"offset":-2.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":true,"shift":0.0,"contact":[-0.07417262089353538,-0.3211550179129533,-2.7911154580885604],"curve":[],"curve_closed":false,"third":null}"#).expect("a scene descriptor");
    let failures: Vec<String> = judge(&scene)
        .failures
        .iter()
        .filter(|f| f.op != Some(Op::BminusA))
        .map(|f| format!("{} {:?}: {}", f.oracle, f.op, f.detail))
        .collect();
    assert!(failures.is_empty(), "{failures:#?}");
}

/// Failed in the baseline (W3 V cells, tangency design §4.8, PK-1292): a box face opened half a tolerance from a tangent crossing on a torus's saddle cut one loop with a waist, the material between its branches thinner than the tolerance, and the union and common served a neck where they pinch (`plane_torus_spiric_crossing`: the waist is spelled as the crossing).
#[test]
fn a_box_face_opened_from_a_torus_saddle_within_the_tolerance_pinches() {
    passes(
        r#"{"cell":{"a":"P","b":"T","config":"V"},"seed":1,"index":10097063,"variant_a":0,"variant_b":0,"b_scale":1.0,"place_b":[[-0.36588415590292733,0.9016654093657778,0.2304957136530889,0.6741711452802353],[0.9078248712434012,0.29125941045829773,0.3016984570252007,-0.05591542132343574],[0.20489701708688784,0.3196364268551525,-0.9251214877064343,0.4147368288868008],[0.0,0.0,0.0,1.0]],"scene":[[0.3281449489708592,0.898428668516936,-0.2917992804168605,-0.5261848644142586],[-0.7337304677411749,0.4369721470635508,0.5202835221290284,2.9870594154714496],[0.5949457901130395,0.0433736126873924,0.8025946900828682,-2.8276664908851004],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-0.5,"rung":{"offset":-0.5,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-1.013491734975748,2.1116066606273365,-2.358281269013464],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline (W3 E cells, tangency design §4.8): a cylinder and a cone touching at scale 1000 read as a near miss, contact discovery's Newton differencing the angle over a sixth of a radian, measuring its gradient in parameter units (the angle's component a thousandfold the length's) and stopping 4e-5 from the touch (`contact_discovery_newton`): the union served where it is non-manifold.
#[test]
fn a_cone_touching_a_cylinder_at_scale_1000_is_found() {
    passes(
        r#"{"cell":{"a":"Y","b":"C","config":"E"},"seed":1,"index":8077,"variant_a":0,"variant_b":2,"b_scale":1.0,"place_b":[[-0.9899772733113584,-0.13559395155459694,-0.03948770224788997,-1.978310464075825],[0.06701819594068406,-0.2049322762273492,-0.9764790441034223,-0.6886645932508695],[0.1243123474955907,-0.9693384560936842,0.21196555804669875,1.2575226242520232],[0.0,0.0,0.0,1.0]],"scene":[[726.0972228303457,525.812542835446,443.0620642697607,601.7868999586549],[-676.0817276354605,663.3721315268426,320.7037771386307,2546.0185732499563],[-125.2849574191185,-532.408287792795,837.1649147772702,-833.7283970866407],[0.0,0.0,0.0,1.0]],"scale":1000.0,"offset":0.0,"rung":{"offset":0.0,"scale":1000.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[690.878097843534,3558.358953643982,-225.41403754263445],"curve":[],"curve_closed":false,"third":null}"#,
    );
}

/// Failed in the baseline (W3 X3 folds, tangency design §4.8): a box face's ellipse on a frustum carried a cone pcurve made of its 3D control points projected to the chart, 1.3e-2 off the section; the clip's crossing read on it landed 7.6e-3 from the junction the record had to reach, the record was dropped from the cone's walk, and the union fold's loop broke (`plane_cone::cone_uv_curve`: the chart fit co-parameterized with the section).
#[test]
fn a_box_face_cutting_a_frustum_keeps_its_record_on_the_cone() {
    passes(
        r#"{"cell":{"a":"P","b":"F","config":"X3"},"seed":1,"index":0,"variant_a":2,"variant_b":1,"b_scale":1.0,"place_b":[[0.02503724152751606,0.9994881747385365,-0.01991293786767423,-1.8361137960772658],[-0.8128230691605336,0.03194923811427587,0.5816338233152951,-0.7410278314740214],[0.5819723316250898,0.0016231887586883875,0.8132069665719243,0.16091725895259557],[0.0,0.0,0.0,1.0]],"scene":[[0.2523748592291189,-0.9161291003879675,0.3114713499688466,1.8525168060287074],[0.9214393829799667,0.32580153736156586,0.21166677052938251,0.8073787184941885],[-0.29539193273104214,0.23358259711539622,0.9263815501197366,2.5834397895213312],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":0.0,"rung":{"offset":0.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":null,"curve":[],"curve_closed":false,"third":{"variant":0,"scale":0.6,"place":[[0.050700391748289,-0.2405566661012032,-0.9693100436242494,-0.6610550492100404],[-0.9321029901226184,0.33711932202689837,-0.1324181955797552,-1.01938997755065],[0.3586272244002733,0.9102104444085437,-0.2071315060741785,-0.7794291969042022],[0.0,0.0,0.0,1.0]]}}"#,
    );
}
