//! Global self-intersection trim stage (PK-327, design §5 Stage 4).
//!
//! Runs post-assembly on every offset/shell pipeline, before the
//! serving gate: a BVH + SSI scan over the assembled result
//! (`popkern_bool::ops::self_intersect`) detects transversal
//! collisions between non-adjacent faces; under
//! [`TrimPolicy::Trim`] the unrealized region is removed and the skin
//! re-closed along the exact intersection curves
//! (`popkern_bool::ops::self_trim`); under [`TrimPolicy::Decline`]
//! the operation declines with reproducible
//! [`TrimRegionReport`]s. The realization predicate is measured:
//! a face survives iff some sample of it achieves its expected offset
//! distance against the SOURCE body's trimmed boundary
//! ([`BodyDistanceIndex`] — the PK-326 substrate).

use std::collections::{BTreeMap, BTreeSet};

use popkern_bool::ops::self_intersect::{SelfIntersectionScan, find_self_intersections};
use popkern_bool::ops::self_trim::trim_self_intersections;
use popkern_core::geometry::curve::CurveEval;
use popkern_core::geometry::surface::SurfaceEval;
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::{BRep, FaceHandle};

use crate::body_distance::{BodyDistanceIndex, FaceTrimEntry};
use crate::error::{ShellError, ShellWarning, TrimRegionReport};
use crate::types::TrimPolicy;

/// One free-form source/offset surface pair for the tier-3 distance
/// map (same parameterization; analytic classes never appear here).
#[derive(Debug)]
pub(crate) struct FreeFormOffsetPair {
    /// The result face carrying the offset surface.
    pub face: FaceHandle,
    /// Source surface.
    pub source: popkern_core::geometry::surface::SurfaceGeometry,
    /// Offset surface (fitted).
    pub offset: popkern_core::geometry::surface::SurfaceGeometry,
    /// |offset distance|.
    pub distance: f64,
    /// Shared u range.
    pub u_range: (f64, f64),
    /// Shared v range.
    pub v_range: (f64, f64),
}

/// Per-face expected offset distances for the realization predicate.
#[derive(Debug, Default)]
pub(crate) struct TrimProvenance {
    /// Result-face → expected |offset distance|. `0.0` means the face
    /// lies ON the source boundary (trivially realized, never
    /// sampled).
    pub expected: BTreeMap<FaceHandle, f64>,
    /// Faces exempt from the realization test (side walls — their
    /// points legitimately span every distance in `[0, d]`).
    pub always_keep: BTreeSet<FaceHandle>,
    /// Expected |distance| for faces absent from `expected` (join
    /// windows, corner spheres, connectors — all rolling-ball surfaces
    /// at the offset radius).
    pub default_distance: f64,
    /// Free-form offset pairs for the tier-3 distance map.
    pub free_form_offsets: Vec<FreeFormOffsetPair>,
}

/// Grid resolution for per-face realization sampling.
const GRID: usize = 5;

/// Run the trim stage on the assembled `result`.
///
/// `deferred` carries the chord-gate decline that was routed here
/// (inverted skin edges): if the scan is clean or the trim cannot
/// serve, that decline is returned so behavior degrades to the honest
/// pre-trim decline. Returns the warnings to append on success.
///
/// # Errors
///
/// [`ShellError::SelfIntersecting`] under the Decline policy or when
/// the trim mechanism fails without a deferred decline to fall back
/// to; the deferred [`ShellError::OffsetTooLarge`] otherwise.
pub(crate) fn run_trim_stage(
    result: &mut BRep,
    source: &BRep,
    provenance: &TrimProvenance,
    policy: TrimPolicy,
    deferred: Option<ShellError>,
    inverted_edges: usize,
    tolerance: &Tolerance,
) -> Result<Vec<ShellWarning>, ShellError> {
    let mut warnings = Vec::new();

    // Full eversion bail: when the majority of the source's edges
    // inverted, the skin passed through itself wholesale — there is
    // no local collision to trim, and scanning the everted body burns
    // minutes of grazing-pair marching before declining anyway
    // (measured 477 s on the over-thick cube). The deferred chord
    // decline is the honest verdict.
    if inverted_edges > 0
        && inverted_edges * 2 >= source.edges.len()
        && let Some(err) = deferred
    {
        return Err(err);
    }

    // Spindle arm first: the pinch class has DEGENERATE trim curves
    // (zero-radius circles) invisible to the transversal scan.
    let spindle_groups = crate::spindle_pinch::detect_spindle_groups(result, tolerance);
    if !spindle_groups.is_empty() {
        let regions = crate::spindle_pinch::spindle_region_reports(&spindle_groups);
        if policy == TrimPolicy::Decline {
            return Err(ShellError::SelfIntersecting {
                regions,
                policy_note: "TrimPolicy::Decline".into(),
            });
        }
        let mut removed = 0usize;
        for group in &spindle_groups {
            match crate::spindle_pinch::serve_spindle_group(result, group, tolerance) {
                Ok(n) => removed += n,
                Err(e) => {
                    return Err(match deferred {
                        Some(err) => err,
                        None => ShellError::SelfIntersecting {
                            regions,
                            policy_note: format!("spindle pinch declined: {e}"),
                        },
                    });
                }
            }
        }
        warnings.push(ShellWarning::RegionTrimmed {
            removed_faces: removed,
            welded_seams: 0,
            regions,
        });
    }

    // Tier 3: numeric distance map on free-form offsets only (the
    // analytic corpus never enters it — pinned). Fold regions decline
    // with reports under BOTH policies: serving a free-form fold needs
    // marched self-intersection curves (decline-class registry row).
    if !provenance.free_form_offsets.is_empty() {
        let mut fold_regions: Vec<TrimRegionReport> = Vec::new();
        for pair in &provenance.free_form_offsets {
            let report = popkern_fillet::offset::distance_map::scan_offset_self_fold(
                &pair.source,
                &pair.offset,
                pair.distance,
                pair.u_range,
                pair.v_range,
            );
            for region in &report.regions {
                use popkern_core::geometry::surface::SurfaceEval;
                let witness = pair.offset.point_at(
                    (region.r_range.0 + region.r_range.1) * 0.5,
                    (region.t_range.0 + region.t_range.1) * 0.5,
                );
                fold_regions.push(TrimRegionReport {
                    face_a: pair.face,
                    face_b: pair.face,
                    witness,
                    curve_length: 0.0,
                });
            }
        }
        if !fold_regions.is_empty() {
            return Err(ShellError::SelfIntersecting {
                regions: fold_regions,
                policy_note: if policy == TrimPolicy::Decline {
                    "TrimPolicy::Decline".into()
                } else {
                    "free-form self-fold detected by the distance map; serving this                      class needs marched self-intersection curves (registry row)"
                        .into()
                },
            });
        }
        warnings.push(ShellWarning::FreeFormFoldScan {
            faces: provenance.free_form_offsets.len(),
        });
    }

    let scan = find_self_intersections(result, tolerance)
        .map_err(|e| ShellError::Topology(format!("self-intersection scan failed: {e}")))?;
    if scan.is_clean() {
        return match deferred {
            // Inverted edges with no transversal curve to trim along:
            // the honest decline stands. (A served spindle already
            // consumed its inverted edges — the spindle warning set is
            // non-empty exactly then.)
            Some(err) if warnings.is_empty() => Err(err),
            _ => Ok(warnings),
        };
    }

    let regions = region_reports(&scan);
    if policy == TrimPolicy::Decline {
        return Err(ShellError::SelfIntersecting {
            regions,
            policy_note: "TrimPolicy::Decline".into(),
        });
    }

    let index = BodyDistanceIndex::build(source);
    let band_tol = *tolerance;
    let keep = |brep: &BRep, fh: FaceHandle, origin: FaceHandle| -> bool {
        if provenance.always_keep.contains(&origin) {
            return true;
        }
        let expected = provenance
            .expected
            .get(&origin)
            .copied()
            .unwrap_or(provenance.default_distance);
        if expected <= band_tol.linear() * 10.0 {
            return true; // on the source boundary: trivially realized
        }
        face_realizes_distance(brep, fh, expected, &index, &band_tol)
    };

    match trim_self_intersections(result, &scan, &keep, tolerance) {
        Ok(outcome) => {
            warnings.push(ShellWarning::RegionTrimmed {
                removed_faces: outcome.removed_faces,
                welded_seams: outcome.welded_seams,
                regions,
            });
            *result = outcome.body;
            Ok(warnings)
        }
        Err(e) => match deferred {
            // Degrade to the honest chord decline, carrying WHY the
            // trim didn't apply (explicit fallback boundary, §7).
            Some(ShellError::OffsetTooLarge {
                requested,
                max_safe,
                binding,
            }) => Err(ShellError::OffsetTooLarge {
                requested,
                max_safe,
                binding: format!("{binding}; trim declined: {e}"),
            }),
            Some(err) => Err(err),
            None => Err(ShellError::SelfIntersecting {
                regions,
                policy_note: format!("trim mechanism declined: {e}"),
            }),
        },
    }
}

/// Reproducible reports for the scan's curves.
fn region_reports(scan: &SelfIntersectionScan) -> Vec<TrimRegionReport> {
    scan.curves
        .iter()
        .map(|c| {
            let (t0, t1) = c.curve.curve_3d.parameter_range();
            let witness = c.curve.curve_3d.point_at((t0 + t1) * 0.5);
            let mut length = 0.0;
            let mut prev = c.curve.curve_3d.point_at(t0);
            for i in 1..=16 {
                let t = t0 + (t1 - t0) * (i as f64) / 16.0;
                let p = c.curve.curve_3d.point_at(t);
                length += (p - prev).length();
                prev = p;
            }
            TrimRegionReport {
                face_a: c.face_a,
                face_b: c.face_b,
                witness,
                curve_length: length,
            }
        })
        .collect()
}

/// Measured realization: does any in-trim sample of `fh` sit at its
/// expected offset distance from the source boundary? Mirrors the
/// PK-326 face-death test, applied to RESULT faces (fragments carry
/// their own trims after the split).
fn face_realizes_distance(
    brep: &BRep,
    fh: FaceHandle,
    expected: f64,
    index: &BodyDistanceIndex,
    tolerance: &Tolerance,
) -> bool {
    let Some(entry) = FaceTrimEntry::build(brep, fh) else {
        return true; // no sampleable trim: keep (conservative)
    };
    let Some((u0, v0, u1, v1)) = entry.uv_bbox() else {
        return true;
    };
    // The split curves ARE the equidistance loci, so every fragment is
    // homogeneously realized or unrealized — but samples NEAR a
    // crossing sit within any fixed band of the expected distance on
    // both sides. Classify by the MEAN realization deficit instead of
    // any-sample-realizes: an unrealized fragment's deficit grows away
    // from the crossing, a realized fragment's stays at numeric zero.
    let band = (tolerance.linear() * 100.0).max(expected * 1.0e-3);
    let mut valid = 0usize;
    let mut deficit_sum = 0.0;
    for i in 0..GRID {
        for j in 0..GRID {
            let u = u0 + (u1 - u0) * (i as f64 + 0.5) / GRID as f64;
            let v = v0 + (v1 - v0) * (j as f64 + 0.5) / GRID as f64;
            if !entry.contains(u, v) {
                continue;
            }
            let p = entry.surface.point_at(u, v);
            let Some(dist) = index.distance(&p) else {
                continue;
            };
            valid += 1;
            deficit_sum += (expected - dist).max(0.0);
        }
    }
    if valid == 0 {
        // No sampleable interior: keep (conservative — the boundary
        // postcondition and the serving validator still guard).
        return true;
    }
    deficit_sum / valid as f64 <= band
}

/// Collect the tier-3 pairs: source faces whose offset surface is
/// free-form (NURBS/Gregory), mapped to their result-face images.
pub(crate) fn collect_free_form_pairs(
    source: &BRep,
    offset_surfaces: &std::collections::HashMap<
        FaceHandle,
        popkern_core::geometry::surface::SurfaceGeometry,
    >,
    image_of: impl Fn(FaceHandle) -> Option<FaceHandle>,
    distance_of: impl Fn(FaceHandle) -> f64,
) -> Vec<FreeFormOffsetPair> {
    use popkern_core::geometry::surface::SurfaceGeometry as SG;
    let mut pairs = Vec::new();
    let mut handles: Vec<FaceHandle> = offset_surfaces.keys().copied().collect();
    handles.sort();
    for fh in handles {
        let offset = &offset_surfaces[&fh];
        let free_form = matches!(
            offset,
            SG::Nurbs(_) | SG::Gregory(_) | SG::PolyhedralSpline(_)
        );
        if !free_form {
            continue;
        }
        let Some(face) = image_of(fh) else { continue };
        let Some(sof) = source.faces.get(fh).and_then(|f| f.surface.as_ref()) else {
            continue;
        };
        let (u_range, v_range) = offset.parameter_range();
        if !u_range.0.is_finite()
            || !u_range.1.is_finite()
            || !v_range.0.is_finite()
            || !v_range.1.is_finite()
        {
            continue;
        }
        pairs.push(FreeFormOffsetPair {
            face,
            source: sof.geometry.clone(),
            offset: offset.clone(),
            distance: distance_of(fh).abs(),
            u_range,
            v_range,
        });
    }
    pairs
}
