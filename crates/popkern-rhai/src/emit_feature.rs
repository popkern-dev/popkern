//! One statement per feature: the canonical Rhai spelling of every
//! feature the skin can record. Features the skin has no spelling for
//! decline with `NotApplicable` naming the feature — never silently
//! dropped.

use popkern_history::feature::{
    AxisDesc, BoolOpKind, ExtrudeEnd, Feature, JoinStyleKind, OperationMode,
    RevolveEnd as TreeRevolveEnd, ShellDirectionKind, TrimPolicyKind,
};
use popkern_history::feature_blend::{
    BlendProfileKind, CapKindDesc, ChamferKindDesc, ContinuityKind, FilletRadiusLaw,
    FilletRecipeDesc, InterferencePolicyKind,
};
use popkern_history::feature_sketch::SketchPayload;
use popkern_history::feature_transform::RigidMotionDesc;
use popkern_history::node::FeatureNode;
use popkern_script::Model;

use crate::emit::EmitError;
use crate::emit_spelling::{Speller, float, plane, quoted, scalar, shape, vec3};

/// The receiver expression and the statement's right-hand side for a
/// feature, plus whether the served body is a new one.
pub struct Statement {
    /// Right-hand side after the receiver (`.extrude(blind(0.4), #{..})`).
    pub call: String,
    /// `true` when the call starts a body (`m.profile(..)...sketch(..)`).
    pub starts_body: bool,
}

fn options(pairs: Vec<(&str, String)>) -> String {
    if pairs.is_empty() {
        String::new()
    } else {
        let inner: Vec<String> = pairs.iter().map(|(k, v)| format!("{k}: {v}")).collect();
        format!(", #{{ {} }}", inner.join(", "))
    }
}

fn label_option(label: &str) -> (&'static str, String) {
    ("label", quoted(label))
}

fn fillet_options(
    model: &Model,
    recipe: &FilletRecipeDesc,
    label: &str,
    context: &str,
) -> Result<String, EmitError> {
    let tree = model.tree();
    let mut pairs = vec![label_option(label)];
    match recipe.continuity {
        ContinuityKind::G1 => {}
        ContinuityKind::G0 => pairs.push(("continuity", quoted("g0"))),
        ContinuityKind::G2 => pairs.push(("continuity", quoted("g2"))),
    }
    match recipe.cap {
        CapKindDesc::Open => {}
        CapKindDesc::Flat => pairs.push(("cap", quoted("flat"))),
        CapKindDesc::Rounded => pairs.push(("cap", quoted("rounded"))),
        CapKindDesc::TangentMatch => pairs.push(("cap", quoted("tangent_match"))),
    }
    match recipe.profile {
        BlendProfileKind::Circular => {}
        BlendProfileKind::Conic { rho } => {
            pairs.push(("profile", format!("conic({})", float(rho))))
        }
        BlendProfileKind::CurvatureContinuous => {
            pairs.push(("profile", quoted("curvature_continuous")))
        }
    }
    if recipe.interference == InterferencePolicyKind::ServeNonInterfering {
        pairs.push(("interference", quoted("serve_non_interfering")));
    }
    match &recipe.radius_law {
        FilletRadiusLaw::Constant => {}
        FilletRadiusLaw::Linear { end } => {
            pairs.push(("law", format!("linear({})", scalar(tree, *end, context)?)))
        }
        FilletRadiusLaw::ByStations(stations) => {
            let items: Vec<String> = stations
                .iter()
                .map(|s| {
                    Ok(format!(
                        "[{}, {}]",
                        float(s.t),
                        scalar(tree, s.radius, context)?
                    ))
                })
                .collect::<Result<_, EmitError>>()?;
            pairs.push(("law", format!("by_stations([{}])", items.join(", "))));
        }
    }
    Ok(options(pairs))
}

fn end_condition(
    model: &Model,
    spell: Speller<'_>,
    end: &ExtrudeEnd,
    context: &str,
) -> Result<String, EmitError> {
    let tree = model.tree();
    Ok(match end {
        ExtrudeEnd::Blind { distance } => format!("blind({})", scalar(tree, *distance, context)?),
        ExtrudeEnd::ThroughAll => "through_all()".into(),
        ExtrudeEnd::UpToSurface { face } => format!(
            "up_to_surface({})",
            spell("end_face", std::slice::from_ref(face), true)?
        ),
        ExtrudeEnd::OffsetFromSurface { face, offset } => format!(
            "offset_from_surface({}, {})",
            spell("end_face", std::slice::from_ref(face), true)?,
            scalar(tree, *offset, context)?
        ),
        ExtrudeEnd::MidPlane { half_distance } => {
            format!("mid_plane({})", scalar(tree, *half_distance, context)?)
        }
        ExtrudeEnd::UpToVertex { .. } => {
            return Err(EmitError::not_applicable(
                context,
                "extrude",
                "up_to_vertex has no spelling in the skin",
            ));
        }
    })
}

/// The body owning the feature that minted `operation` (loft sections
/// name their sketch's operation; that sketch's node carries the body).
fn body_of_operation(
    tree: &popkern_history::FeatureTree,
    operation: popkern_core::naming::OperationId,
) -> Option<popkern_history::node::BodyId> {
    tree.timeline()
        .iter()
        .filter_map(|id| tree.feature(*id))
        .find(|n| n.operation_id == operation)
        .map(|n| n.body)
}

fn axis(spell: Speller<'_>, axis: &AxisDesc) -> Result<String, EmitError> {
    Ok(match axis {
        AxisDesc::Literal { origin, direction } => {
            format!("axis_at({}, {})", vec3(*origin), vec3(*direction))
        }
        AxisDesc::Edge(edge) => {
            format!(
                "axis({})",
                spell("axis", std::slice::from_ref(edge), false)?
            )
        }
    })
}

/// The statement for `node`. `tool_var` names the tool body of a
/// boolean.
pub fn statement(
    model: &Model,
    node: &FeatureNode,
    tool_var: &dyn Fn(popkern_history::node::BodyId) -> String,
    spell: Speller<'_>,
) -> Result<Statement, EmitError> {
    let tree = model.tree();
    let label = node.label.as_str();
    let context = label;
    let op = node.feature.operation_label();
    let decline = |detail: &str| EmitError::not_applicable(label, op, detail);
    if node.suppress_if.is_some() {
        return Err(decline("suppress_if has no spelling in the skin"));
    }
    let call = match &node.feature {
        Feature::Sketch {
            plane: p,
            sketch_data,
        } => {
            let profile = match sketch_data {
                SketchPayload::Analytic(profile) => profile,
                SketchPayload::Constraint(_) => {
                    return Err(decline(
                        "constraint sketches are not spelled by the skin yet",
                    ));
                }
            };
            return Ok(Statement {
                call: format!(
                    "m.profile({}){}.sketch({})",
                    plane(model, p, spell, context)?,
                    shape(tree, profile, context)?,
                    quoted(label)
                ),
                starts_body: true,
            });
        }
        Feature::Extrude {
            direction,
            end_condition: end,
            mode,
            draft_angle,
            second_direction,
            ..
        } => {
            if direction.is_some() || draft_angle.is_some() || second_direction.is_some() {
                return Err(decline(
                    "extrude direction / draft / second direction have no spelling in the skin",
                ));
            }
            let mut pairs = vec![label_option(label)];
            match mode {
                OperationMode::NewBody => {}
                OperationMode::Boss => pairs.push(("mode", quoted("boss"))),
                OperationMode::Cut => pairs.push(("mode", quoted("cut"))),
            }
            format!(
                ".extrude({}{})",
                end_condition(model, spell, end, context)?,
                options(pairs)
            )
        }
        Feature::Revolve {
            axis: a,
            end_condition: end,
            mode,
            ..
        } => {
            if !matches!(mode, OperationMode::NewBody) {
                return Err(decline(
                    "revolve boss/cut modes have no spelling in the skin",
                ));
            }
            let end = match end {
                TreeRevolveEnd::Blind { angle } => {
                    format!("blind({})", scalar(tree, *angle, context)?)
                }
                TreeRevolveEnd::FullCircle => "full_circle()".into(),
                TreeRevolveEnd::MidPlane { half_angle } => {
                    format!("mid_plane({})", scalar(tree, *half_angle, context)?)
                }
                TreeRevolveEnd::UpToSurface { .. } => {
                    return Err(decline("revolve up_to_surface has no spelling in the skin"));
                }
            };
            format!(
                ".revolve({}, {}{})",
                axis(spell, a)?,
                end,
                options(vec![label_option(label)])
            )
        }
        Feature::Sweep { mode, .. } => {
            if !matches!(mode, OperationMode::NewBody) {
                return Err(decline("sweep boss/cut modes have no spelling in the skin"));
            }
            // The path edge lives on another body; spelling it needs
            // cross-body selection verification the emitter does not
            // have yet (the recorded query would be re-verified against
            // the wrong body's state).
            return Err(decline(
                "sweep paths cross bodies; the skin has no verified spelling for them yet",
            ));
        }
        Feature::Loft {
            sections,
            closed,
            mode,
        } => {
            if !matches!(mode, OperationMode::NewBody) {
                return Err(decline("loft boss/cut modes have no spelling in the skin"));
            }
            let vars = sections
                .iter()
                .map(|s| {
                    body_of_operation(tree, s.operation)
                        .map(&tool_var)
                        .ok_or_else(|| {
                            EmitError::not_applicable(
                                label,
                                op,
                                "a loft section's sketch is not in the timeline",
                            )
                        })
                })
                .collect::<Result<Vec<_>, _>>()?;
            let mut pairs = Vec::new();
            if *closed {
                pairs.push(("closed", "true".to_string()));
            }
            pairs.push(label_option(label));
            // Model-level call: sections are sketch bodies, so the
            // statement starts (re-binds) the first section's body.
            return Ok(Statement {
                call: format!("m.loft([{}]{})", vars.join(", "), options(pairs)),
                starts_body: true,
            });
        }
        Feature::Pipe {
            path_edge,
            radius,
            mode,
        } => {
            if !matches!(mode, OperationMode::NewBody) {
                return Err(decline("pipe boss/cut modes have no spelling in the skin"));
            }
            format!(
                ".pipe({}, {}{})",
                spell("path_edge", std::slice::from_ref(path_edge), false)?,
                scalar(tree, *radius, context)?,
                options(vec![label_option(label)])
            )
        }
        Feature::HelicalSweep {
            axis: a,
            start_radius,
            end_radius,
            pitch,
            turns,
            start_angle,
            left_handed,
            mode,
            ..
        } => {
            if !matches!(mode, OperationMode::NewBody) {
                return Err(decline(
                    "helical_sweep boss/cut modes have no spelling in the skin",
                ));
            }
            let r = scalar(tree, *start_radius, context)?;
            let taper = scalar(tree, *end_radius, context)?;
            let phase = scalar(tree, *start_angle, context)?;
            let mut pairs = Vec::new();
            if taper != r {
                pairs.push(("taper_to", taper));
            }
            if phase != "0" {
                pairs.push(("start_angle", phase));
            }
            if *left_handed {
                pairs.push(("left_hand", "true".to_string()));
            }
            pairs.push(label_option(label));
            format!(
                ".helical_sweep({}, {}, {}, {}{})",
                axis(spell, a)?,
                r,
                scalar(tree, *pitch, context)?,
                scalar(tree, *turns, context)?,
                options(pairs)
            )
        }
        Feature::Fillet {
            edges,
            radius,
            recipe,
        } => format!(
            ".fillet({}, {}{})",
            spell("edges", edges, false)?,
            scalar(tree, *radius, context)?,
            fillet_options(model, recipe, label, context)?
        ),
        Feature::Chamfer {
            edges,
            distance,
            kind,
        } => {
            let mut pairs = vec![label_option(label)];
            match kind {
                ChamferKindDesc::Equal => {}
                ChamferKindDesc::TwoDistance { face, distance_2 } => pairs.push((
                    "kind",
                    format!(
                        "two_distance({}, {})",
                        spell("face", std::slice::from_ref(face), true)?,
                        scalar(tree, *distance_2, context)?
                    ),
                )),
                ChamferKindDesc::DistanceAngle { face, angle } => pairs.push((
                    "kind",
                    format!(
                        "distance_angle({}, {})",
                        spell("face", std::slice::from_ref(face), true)?,
                        scalar(tree, *angle, context)?
                    ),
                )),
            }
            format!(
                ".chamfer({}, {}{})",
                spell("edges", edges, false)?,
                scalar(tree, *distance, context)?,
                options(pairs)
            )
        }
        Feature::Shell {
            faces_to_remove,
            thickness,
            direction,
            join_style,
            trim_policy,
        } => {
            let mut pairs = vec![label_option(label)];
            match direction {
                ShellDirectionKind::Inward => {}
                ShellDirectionKind::Outward => pairs.push(("direction", quoted("outward"))),
                ShellDirectionKind::Midplane => pairs.push(("direction", quoted("midplane"))),
            }
            if matches!(join_style, JoinStyleKind::Intersection) {
                pairs.push(("join", quoted("intersection")));
            }
            if matches!(trim_policy, TrimPolicyKind::Decline) {
                pairs.push(("trim", quoted("decline")));
            }
            format!(
                ".shell({}, {}{})",
                spell("faces_to_remove", faces_to_remove, true)?,
                scalar(tree, *thickness, context)?,
                options(pairs)
            )
        }
        Feature::Boolean {
            operation,
            tool_body,
        } => {
            let verb = match operation {
                BoolOpKind::Union => "union",
                BoolOpKind::Subtraction => "subtract",
                BoolOpKind::Intersection => "intersect",
            };
            format!(
                ".{verb}({}{})",
                tool_var(*tool_body),
                options(vec![label_option(label)])
            )
        }
        Feature::Transform { motion } => match motion {
            RigidMotionDesc::Translate { dx, dy, dz } => format!(
                ".translate({}, {}, {}{})",
                scalar(tree, *dx, context)?,
                scalar(tree, *dy, context)?,
                scalar(tree, *dz, context)?,
                options(vec![label_option(label)])
            ),
            RigidMotionDesc::Rotate { axis: a, angle } => format!(
                ".rotate({}, {}{})",
                axis(spell, a)?,
                scalar(tree, *angle, context)?,
                options(vec![label_option(label)])
            ),
            _ => return Err(decline("general transforms have no spelling in the skin")),
        },
        _ => return Err(decline("this operation has no spelling in the Rhai skin")),
    };
    Ok(Statement {
        call,
        starts_body: false,
    })
}
