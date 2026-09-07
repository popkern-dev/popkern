//! Round-trip emit for the fastener features (PK-460): the
//! `Feature::Fastener` model-level statement and the
//! `TappedHole` / `ClearanceHole` body-chained calls. Threads emit
//! in the `#{ major, pitch }` spec-map form — value-faithful to the
//! recorded parameters (a designation string is a record-time
//! convenience, not stored).

use popkern_history::FeatureTree;
use popkern_history::feature::Feature;
use popkern_history::feature_fastener::{ClearanceFitDesc, FastenerKind, ThreadDesc};

use super::{Statement, axis, label_option, options};
use crate::emit::EmitError;
use crate::emit_spelling::{Speller, quoted, scalar};

fn spec_map(tree: &FeatureTree, thread: &ThreadDesc, context: &str) -> Result<String, EmitError> {
    Ok(format!(
        "#{{ major: {}, pitch: {} }}",
        scalar(tree, thread.major_diameter, context)?,
        scalar(tree, thread.pitch, context)?
    ))
}

fn thread_pairs(thread: &ThreadDesc, pairs: &mut Vec<(&'static str, String)>) {
    if thread.left_handed {
        pairs.push(("left_hand", "true".to_string()));
    }
    if thread.starts > 1 {
        pairs.push(("starts", thread.starts.to_string()));
    }
}

/// The model-level `m.bolt(..)` / `m.threaded_rod(..)` / `m.nut(..)`
/// statement.
pub(super) fn fastener_statement(
    tree: &FeatureTree,
    context: &str,
    label: &str,
    feature: &Feature,
) -> Result<Statement, EmitError> {
    let Feature::Fastener {
        kind,
        thread,
        length,
        thread_length,
        clearance,
        tip_chamfer,
    } = feature
    else {
        return Err(EmitError::not_applicable(
            label,
            "fastener",
            "not a fastener feature",
        ));
    };
    let spec = spec_map(tree, thread, context)?;
    let mut pairs = Vec::new();
    match kind {
        FastenerKind::SocketCapScrew => pairs.push(("head", quoted("socket"))),
        FastenerKind::CountersunkScrew => pairs.push(("head", quoted("countersunk"))),
        FastenerKind::HexBolt | FastenerKind::ThreadedRod | FastenerKind::HexNut => {}
    }
    if let Some(id) = thread_length {
        pairs.push(("thread_length", scalar(tree, *id, context)?));
    }
    if let Some(id) = clearance {
        pairs.push(("clearance", scalar(tree, *id, context)?));
    }
    thread_pairs(thread, &mut pairs);
    if !tip_chamfer && !matches!(kind, FastenerKind::HexNut) {
        pairs.push(("tip_chamfer", "false".to_string()));
    }
    let call = match kind {
        FastenerKind::HexNut => {
            let thickness = scalar(tree, *length, context)?;
            if thickness != "0" {
                pairs.push(("thickness", thickness));
            }
            pairs.push(label_option(label));
            format!("m.nut({spec}{})", options(pairs))
        }
        FastenerKind::ThreadedRod => {
            pairs.push(label_option(label));
            format!(
                "m.threaded_rod({spec}, {}{})",
                scalar(tree, *length, context)?,
                options(pairs)
            )
        }
        _ => {
            pairs.push(label_option(label));
            format!(
                "m.bolt({spec}, {}{})",
                scalar(tree, *length, context)?,
                options(pairs)
            )
        }
    };
    Ok(Statement {
        call,
        starts_body: true,
    })
}

/// The body-chained `.tapped_hole(..)` / `.clearance_hole(..)` call.
pub(super) fn hole_call(
    tree: &FeatureTree,
    context: &str,
    label: &str,
    spell: Speller<'_>,
    feature: &Feature,
) -> Result<String, EmitError> {
    match feature {
        Feature::TappedHole {
            axis: a,
            thread,
            depth,
            clearance,
            entry_chamfer,
        } => {
            let spec = spec_map(tree, thread, context)?;
            let depth = match depth {
                Some(id) => format!("blind({})", scalar(tree, *id, context)?),
                None => "through_all()".to_string(),
            };
            let mut pairs = Vec::new();
            let clearance_text = scalar(tree, *clearance, context)?;
            if clearance_text != "0" {
                pairs.push(("clearance", clearance_text));
            }
            thread_pairs(thread, &mut pairs);
            if !entry_chamfer {
                pairs.push(("entry_chamfer", "false".to_string()));
            }
            pairs.push(label_option(label));
            Ok(format!(
                ".tapped_hole({}, {spec}, {depth}{})",
                axis(spell, a)?,
                options(pairs)
            ))
        }
        Feature::ClearanceHole {
            axis: a,
            thread,
            fit,
            depth,
            countersink,
            counterbore,
            entry_chamfer,
        } => {
            let spec = spec_map(tree, thread, context)?;
            let depth = match depth {
                Some(id) => format!("blind({})", scalar(tree, *id, context)?),
                None => "through_all()".to_string(),
            };
            let mut pairs = Vec::new();
            match fit {
                ClearanceFitDesc::Normal => {}
                ClearanceFitDesc::Close => pairs.push(("fit", quoted("close"))),
                ClearanceFitDesc::Loose => pairs.push(("fit", quoted("loose"))),
                ClearanceFitDesc::Custom(id) => {
                    pairs.push(("diameter", scalar(tree, *id, context)?));
                }
            }
            if let Some(id) = countersink {
                pairs.push(("countersink", scalar(tree, *id, context)?));
            }
            if let Some(id) = counterbore {
                pairs.push(("counterbore", scalar(tree, *id, context)?));
            }
            if !entry_chamfer {
                pairs.push(("entry_chamfer", "false".to_string()));
            }
            pairs.push(label_option(label));
            Ok(format!(
                ".clearance_hole({}, {spec}, {depth}{})",
                axis(spell, a)?,
                options(pairs)
            ))
        }
        _ => Err(EmitError::not_applicable(
            label,
            "fastener_hole",
            "not a fastener hole feature",
        )),
    }
}
