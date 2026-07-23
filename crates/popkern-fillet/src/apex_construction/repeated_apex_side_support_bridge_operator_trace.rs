//! PK-167 stage 1b observability: dump the face-family bridge cycle
//! structure at execution time (post-exposure) under
//! `PK_SLIR_TRACE=1`. The in-place support-loop rewrite derivation
//! needs the measured near/far classification of the source loop's
//! sides — this trace is how that structure was captured on the
//! V-free mixed-cube fixture (each shared-start source is a 3-sided
//! corner fragment: two dead sides where the corner material was
//! removed plus one rail side already twinning its strip).

use popkern_core::topology::BRep;

use super::super::repeated_apex_side_support_bridge::FaceFamilyBridgePrimitivePlan;

pub(super) fn trace_cycle_structure(primitive: &FaceFamilyBridgePrimitivePlan, candidate: &BRep) {
    if std::env::var_os("PK_SLIR_TRACE").is_none() {
        return;
    }
    for region in primitive.regions() {
        for (ci, cycle) in region.boundary_cycles().iter().enumerate() {
            for (label, sides) in [
                ("support", cycle.support_interval()),
                ("outside", cycle.outside_interval()),
            ] {
                for &h in sides {
                    let Some(data) = candidate.half_edges.get(h) else {
                        eprintln!("[slir] c{ci} {label} {h:?} STALE");
                        continue;
                    };
                    let Some(v) = candidate.vertices.get(data.origin) else {
                        eprintln!("[slir] c{ci} {label} {h:?} stale-origin");
                        continue;
                    };
                    let p = v.position;
                    let twin = candidate.half_edges.get(data.twin).map(|t| {
                        let f = candidate.loops.get(t.loop_).map(|l| l.face);
                        let kind = f.and_then(|f| candidate.faces.get(f)).map(|fd| {
                            (
                                fd.role,
                                fd.surface
                                    .as_ref()
                                    .map(|s| std::mem::discriminant(&s.geometry)),
                            )
                        });
                        format!("{f:?} {kind:?}")
                    });
                    eprintln!(
                        "[slir] c{ci} {label} {h:?} o=({:.3},{:.3},{:.3}) twin {}",
                        p.x,
                        p.y,
                        p.z,
                        twin.unwrap_or_else(|| "DEAD".into())
                    );
                }
            }
        }
    }
}
