//! K cells: operands tangent along a curve (generated corpus §3,
//! W3 K cells).

use super::passes;

/// Failed at 4cfacba0 (W3 K cells): a cone lying on a box face along a
/// ruling, sunk ten tolerances, its apex on the face. The face-pair
/// solver's window on the cone padded its v range past the apex row into
/// the mirrored nappe and traced the plane's two lines there, none on
/// the face's own nappe: the cone's face walked in 3D with no records,
/// was never selected, and the union left the box's edge on four faces.
/// A finite end of a half-infinite chart range is never padded past
/// (`uv_projection_padding::padded_within`).
#[test]
fn a_cone_on_a_box_along_a_ruling_meets_it_on_its_own_nappe() {
    passes(
        r#"{"cell":{"a":"P","b":"C","config":"K"},"seed":1,"index":10119261,"variant_a":0,"variant_b":0,"b_scale":0.25660011963983365,"place_b":[[-0.4472135954999579,0.0,0.8944271909999159,0.0],[0.0,1.0,0.0,-0.007964401950836275],[-0.8944271909999159,0.0,-0.4472135954999579,0.7032895558988643],[0.0,0.0,0.0,1.0]],"scene":[[0.3221974391280571,0.39728463182107643,0.8592751198179277,-0.34922751467367075],[0.027141237999861206,-0.9111829028639623,0.41110712804357685,1.3517640426778161],[0.9462833420382909,-0.10913587332859301,-0.3043635946259807,-2.0309546809258476],[0.0,0.0,0.0,1.0]],"scale":1.0,"offset":10.0,"rung":{"offset":10.0,"scale":1.0,"steer_a":"Free","steer_b":"Free"},"aligned":false,"shift":0.0,"contact":[0.14653521850792423,1.6042837788895152,-2.2615694973990954],"curve":[[0.06334419744568137,1.5972759413715876,-2.5058988056996157],[0.08414195271124209,1.5990279007510695,-2.4448164786244853],[0.1049397079768028,1.6007798601305514,-2.3837341515493558],[0.1257374632423635,1.602531819510033,-2.3226518244742254],[0.14653521850792423,1.6042837788895152,-2.2615694973990954],[0.167332973773485,1.606035738268997,-2.2004871703239655],[0.1881307290390457,1.607787697648479,-2.1394048432488355],[0.20892848430460642,1.6095396570279608,-2.0783225161737056],[0.22972623957016713,1.6112916164074427,-2.0172401890985756]],"curve_closed":false,"third":null}"#,
    );
}
