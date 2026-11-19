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
