//! A cone face whose boundary turns round the apex a few tolerances from
//! it, meshed in the cone's development (generated corpus PC-K: a plane
//! missing a cone's apex by about ten tolerances).

use super::passes;

/// Failed at bb7b0667 (generated corpus PC-K 292): the common's cone face
/// is a ruling sliver whose conic turns about pi round the apex 1.6e-5
/// from it; in the face chart the turn runs along the apex row, the
/// constrained recovery declined and the fallback filled the sliver's
/// complement (mesh 0.082 for 1.4e-8). The development
/// (`boundary_first_cone_development`) meshes it at its true size.
#[test]
fn a_cone_sliver_turning_round_its_apex_meshes_its_own_area() {
    passes(
        r#"{"cell":{"a":"P","b":"C","config":"K"},"seed":1,"index":10092018,"variant_a":1,"variant_b":0,"b_scale":0.3047126420723025,"place_b":[[-0.4472135954999579,0.0,0.8944271909999159,0.0],[0.0,1.0,0.0,-0.01269616920685448],[-0.8944271909999159,0.0,-0.4472135954999579,0.9226344726299016],[0.0,0.0,0.0,1.0]],"scene":[[-0.037405012580212005,-0.4421311358920387,-0.8961701421653647,-0.09905730640635113],[0.9992885319603254,-0.02088071831153432,-0.03140741147200993,1.6761786644000622],[-0.004826481788214265,-0.8967073403723258,0.44259761725063335,1.8703286254668114],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-10.0,"rung":{"offset":-10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-0.8080862894049314,1.5900392328516069,2.2360874380465003],"curve":[[-0.796617494300894,1.2836461953331504,2.237567291330434],[-0.7994846930769034,1.3602444547127646,2.2371973280094504],[-0.8023518918529127,1.4368427140923785,2.236827364688467],[-0.805219090628922,1.5134409734719927,2.236457401367484],[-0.8080862894049314,1.5900392328516069,2.2360874380465003],[-0.8109534881809407,1.666637492231221,2.235717474725517],[-0.8138206869569501,1.7432357516108352,2.235347511404534],[-0.8166878857329594,1.8198340109904494,2.2349775480835503],[-0.8195550845089687,1.8964322703700636,2.234607584762567]],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at bb7b0667 (generated corpus PC-K 306): the cone face beside
/// the sliver spans a whole period between its seam's two uses, which
/// meet at the apex, its conic passing the apex 1.2e-6 off; the
/// development takes the ring's period where it passes the apex (the
/// seam's uses two rays from the origin), and samples the face chart
/// spread along the apex row at one position merge.
#[test]
fn a_cone_face_closing_round_its_apex_through_it_meshes_one_sheet() {
    passes(
        r#"{"cell":{"a":"P","b":"C","config":"K"},"seed":1,"index":10109171,"variant_a":1,"variant_b":1,"b_scale":0.24931034351370202,"place_b":[[-0.4472135954999579,0.0,0.8944271909999159,0.0],[0.0,1.0,0.0,-0.03675292335594439],[-0.8944271909999159,0.0,-0.4472135954999579,0.9226424726299015],[0.0,0.0,0.0,1.0]],"scene":[[-0.16738405741615373,-0.7881750097786117,-0.5922522530842675,-1.6740504850700526],[0.9467851462878439,0.03900687175319534,-0.319493897789513,-0.6134126632295094],[0.2749190136960552,-0.614213821007404,0.7397032634725043,-0.09159746448615635],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-2.0,"rung":{"offset":-2.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-2.108620187087774,-0.9285003805091685,0.5058806510541513],"curve":[[-2.0572983634666304,-1.2187952929099106,0.42158740742743783],[-2.0701288193719165,-1.146221564809725,0.44266071833411613],[-2.082959275277202,-1.0736478367095397,0.46373402924079454],[-2.0957897311824882,-1.001074108609354,0.48480734014747284],[-2.108620187087774,-0.9285003805091685,0.5058806510541513],[-2.12145064299306,-0.855926652408983,0.5269539619608296],[-2.1342810988983456,-0.7833529243087975,0.548027272867508],[-2.1471115548036313,-0.710779196208612,0.5691005837741863],[-2.1599420107089173,-0.6382054681084265,0.5901738946808647]],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ec1c1568 (generated corpus PC-K 311): the union's cone face
/// holds the apex, its outer ring an ellipse's near end sweeping round it
/// 9e-4 off, turning once round the axis without passing the apex; the
/// face chart's apex cap met the sweep along its apex row and the mesh's
/// winding votes turned the face inside out (the union's B-rep volume
/// 3.8929 over V(A) + V(B) = 3.8658). The development cuts it along a
/// ruling to the apex (`boundary_first_cone_development_cut`).
#[test]
fn a_cone_face_holding_its_apex_beside_a_sweep_is_cut_to_it() {
    passes(
        r#"{"cell":{"a":"P","b":"C","config":"K"},"seed":1,"index":10115225,"variant_a":0,"variant_b":0,"b_scale":0.25660011963983365,"place_b":[[-0.4472135954999579,0.0,0.8944271909999159,0.0],[0.0,1.0,0.0,-0.035940862130413136],[-0.8944271909999159,0.0,-0.4472135954999579,0.7032695558988644],[0.0,0.0,0.0,1.0]],"scene":[[-0.6057148745815699,-0.7825017206507303,0.1442239504703686,-0.9956990011146551],[-0.6168456675064917,0.5762916337709043,0.536087096770788,1.1797180333541997],[-0.5026041316872589,0.235752009586998,-0.8317512108717297,-1.8836201478850714],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-10.0,"rung":{"offset":-10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-0.8497618627496546,1.5125116465595063,-2.365189639124846],"curve":[[-0.6933669546293523,1.6717805130550292,-2.2354178103408597],[-0.7324656816594278,1.6319632964311483,-2.2678607675368565],[-0.7715644086895035,1.5921460798072677,-2.300303724732853],[-0.810663135719579,1.5523288631833871,-2.3327466819288496],[-0.8497618627496546,1.5125116465595063,-2.365189639124846],[-0.8888605897797301,1.4726944299356257,-2.3976325963208427],[-0.9279593168098057,1.4328772133117451,-2.430075553516839],[-0.9670580438398813,1.3930599966878645,-2.462518510712836],[-1.0061567708699568,1.3532427800639837,-2.494961467908832]],"curve_closed":false,"third":null}"#,
    );
}

/// Failed at ec1c1568 (generated corpus PC-K 300): the common's cone face
/// is a thin cap round the apex split along the seam, a hyperbola's arms
/// sweeping round the apex within 1e-4 of it, its split vertex 5.5e-6
/// off; one ring step from there turns -4.77 along the curve but +1.51
/// between its positions, so the development took the short way and the
/// face chart met the sweep (mesh 0.036 for 1.5e-9). The step is read
/// from the edge's own curve (`boundary_first_cone_development_azimuth`).
#[test]
fn a_cone_cap_whose_arms_sweep_round_its_apex_reads_its_curves() {
    passes(
        r#"{"cell":{"a":"P","b":"C","config":"K"},"seed":1,"index":10100090,"variant_a":2,"variant_b":0,"b_scale":0.3528251645047713,"place_b":[[-0.4472135954999579,0.0,0.8944271909999159,0.0],[0.0,1.0,0.0,-0.032374467339825186],[-0.8944271909999159,0.0,-0.4472135954999579,1.1420073893609386],[0.0,0.0,0.0,1.0]],"scene":[[1.0,0.0,0.0,0.0],[0.0,1.0,0.0,0.0],[0.0,0.0,1.0,0.0],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":-2.0,"rung":{"offset":-2.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[-0.07100469468046988,-0.032374467339824964,1.0000000000000004],"curve":[[-0.4260281680828166,-0.032374467339825186,0.9999999999999996],[-0.3372722997322297,-0.03237446733982541,0.9999999999999996],[-0.2485164313816428,-0.032374467339824964,0.9999999999999996],[-0.15976056303105635,-0.03237446733982474,0.9999999999999996],[-0.07100469468046988,-0.032374467339824964,1.0000000000000004],[0.017751173670117026,-0.032374467339824964,1.0],[0.10650704202070393,-0.032374467339825186,1.0],[0.1952629103712904,-0.032374467339824964,1.0],[0.2840187787218773,-0.032374467339824964,1.0]],"curve_closed":false,"third":null}"#,
    );
}
