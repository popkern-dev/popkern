//! The two anti-slide lip pieces: the slide-in dovetail lip (mounts
//! in the arm tips' through-slots) and the peg lip (the holder-v1 W3
//! design — an L-bracket on a hexagram peg shank for any spare
//! stacking socket).

use popkern_part::prelude::*;

use super::{StandParams, hexagram};

/// The slide-in lip: dovetail tab at fit clearance + blade, one
/// 10-gon profile × one extrude (no booleans). The slot is a
/// through-slot so identical arms serve both sides; the lip inserts
/// from the outboard edge.
pub(crate) fn build_lip(p: &StandParams) -> Result<Part, PartError> {
    let c = p.fit_clearance;
    // The blade's two long outer corners — the edges the laptop
    // leans on — are rounded IN THE SKETCH (r = 2 arcs), so the
    // extrude carries true cylinder fillets by construction: no
    // dress-up pass, no blend-corner patches on the non-convex end
    // caps (a post-extrude fillet left the caps' corner trims open —
    // the flagship's blend-corner frontier). The dovetail tab stays
    // sharp for the joint.
    let r = 2.0;
    let d = r - r / 2.0_f64.sqrt(); // arc mid-point sagitta offset
    // Tab = the arm's slot dovetail inset by a UNIFORM perpendicular
    // clearance c (slot flank (1,4.5)->(6,7), slope 0.5), so the
    // flanks are exactly parallel; nose keeps the 0.4 mm single-stop
    // setback (the blade shoulder at x = 0 is the only datum).
    let slope = 0.5_f64;
    let flank_shift = c * (1.0 + slope * slope).sqrt();
    let x_throat = 1.0 + (flank_shift - c) / slope;
    let y_nose = 4.5 + slope * (5.6 - 1.0) - flank_shift;
    // The throat land is not a working surface (flanks carry the
    // fit), so it tapers straight to the blade edge — no sharp
    // half-width step at the top of the tab.
    let lip = Profile::on(Plane::XY)
        .path(x_throat, -(4.5 - c))
        .line_to(5.6, -y_nose)
        .line_to(5.6, y_nose)
        .line_to(x_throat, 4.5 - c)
        .line_to(0.0, 4.0)
        .line_to(-17.0 + r, 4.0)
        .arc_to((-17.0 + d, 4.0 - d), (-17.0, 4.0 - r))
        .line_to(-17.0, -4.0 + r)
        .arc_to((-17.0 + d, -4.0 + d), (-17.0 + r, -4.0))
        .line_to(0.0, -4.0)
        .close()
        .extrude(40.0)?
        .rotate_y(90.0_f64.to_radians())?;
    println!("  lip: print on a side face; two per holder");
    Ok(lip)
}

/// Anti-slide PEG LIP — the holder-v1 W3 design, shippable since the
/// star-boss union closed: an L-bracket (socket-covering base pad +
/// rounded-top blade wall) with a hexagram peg shank underneath, so
/// it plugs into any spare stacking socket as a rear stop. The shank
/// embeds 1 mm into the pad (transversal union — the star walls
/// cross the pad underside; the proven boss-onto-plate class) and
/// carries a chamfer lead-in on its free rim, dressed on the tool
/// BEFORE the union (dress-first).
pub(crate) fn build_peg_lip(p: &StandParams) -> Result<Part, PartError> {
    // L-profile in (depth x, height y): pad 34 deep × 4 tall, blade
    // wall 5 thick rising to 30 at the rear with a HALF-ROUND top —
    // one 180° sketch arc → one true cylinder face by construction.
    // (Two corner arcs flanking a 1 mm flat exported a merged top
    // surface OCC's reader rejects — "Entity has illegal type", the
    // face dropped and its quad left as free bounds; single-arc
    // crowns are the proven class.)
    let bracket = Profile::on(Plane::XY)
        .path(-17.0, 0.0)
        .line_to(17.0, 0.0)
        .line_to(17.0, 27.5)
        .arc_to((14.5, 30.0), (12.0, 27.5))
        .line_to(12.0, 4.0)
        .line_to(-17.0, 4.0)
        .close()
        .extrude(34.0)?
        .rotate_x(90.0_f64.to_radians())?
        .translate(0.0, 17.0, 0.0)?;
    let tip_z = -p.socket_depth;
    let shank = Profile::on(Plane::xy_at(tip_z))
        .polygon(&hexagram(p.star_r, p.fit_clearance))
        .extrude(p.socket_depth + 1.0)?
        .chamfer(p.lead_in, move |e: EdgeQuery| {
            e.straight()
                .at_height(tip_z)
                .between(SurfaceKind::Plane, SurfaceKind::Plane)
        })?;
    let lip = bracket.union(&shank)?;
    println!("  peg lip: print blade-down; plugs any spare stacking socket");
    Ok(lip)
}
