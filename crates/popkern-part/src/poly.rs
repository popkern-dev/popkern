//! Polyhedral modeling: convex hulls, vertex truncation, and the
//! generic planar-face BRep builder — promoted from the dice-set
//! example with its conventions intact.
//!
//! Representation: a polyhedron is `(Vec<Point3>, Vec<PolyFace>)`
//! where each face is an outer vertex-index cycle wound CCW viewed
//! from OUTSIDE plus optional inner cycles wound CW (engraved
//! pockets). The builder wires half-edge topology, attaches analytic
//! line curves and right-handed plane charts (v_axis = normal ×
//! u_axis, so STEP export needs no sense compensation), recomputes
//! pcurves, and validates.

use std::collections::HashMap;

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, SurfaceOnFace, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::*;

use crate::error::PartError;

/// Hull-membership epsilon on unit-scale canonical coordinates.
const HULL_EPS: f64 = 1e-9;

/// A polyhedron face: outer boundary plus engraved inner loops, as
/// vertex-index cycles (outer CCW from outside, inners CW).
#[derive(Debug, Clone)]
pub struct PolyFace {
    /// Outer cycle, CCW viewed from outside the solid.
    pub outer: Vec<usize>,
    /// Inner (hole) cycles, CW viewed from outside.
    pub inners: Vec<Vec<usize>>,
}

impl PolyFace {
    /// A face with no inner loops.
    pub fn simple(outer: Vec<usize>) -> PolyFace {
        PolyFace {
            outer,
            inners: Vec::new(),
        }
    }
}

fn origin() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}

/// Faces of the convex hull of `points` (assumed in convex
/// position), coplanar triangles merged into polygons, each wound
/// CCW viewed from outside. O(n⁴) plane enumeration — intended for
/// canonical solids (n ≤ a few dozen), not point clouds.
pub fn convex_hull_faces(points: &[Point3]) -> Vec<Vec<usize>> {
    let n = points.len();
    let mut faces: Vec<(Vec3, f64, Vec<usize>)> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let normal = (points[j] - points[i]).cross(&(points[k] - points[i]));
                let Some(mut normal) = normal.normalized() else {
                    continue;
                };
                let mut offset = normal.dot(&(points[i] - origin()));
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                for p in points {
                    let d = normal.dot(&(*p - origin())) - offset;
                    lo = lo.min(d);
                    hi = hi.max(d);
                }
                if hi > HULL_EPS && lo < -HULL_EPS {
                    continue;
                }
                if hi > HULL_EPS {
                    normal = normal * -1.0;
                    offset = -offset;
                }
                if faces.iter().any(|(fn_, fo, _)| {
                    fn_.dot(&normal) > 1.0 - 1e-7 && (fo - offset).abs() < HULL_EPS
                }) {
                    continue;
                }
                let members: Vec<usize> = (0..n)
                    .filter(|&m| (normal.dot(&(points[m] - origin())) - offset).abs() < HULL_EPS)
                    .collect();
                faces.push((normal, offset, members));
            }
        }
    }
    faces
        .into_iter()
        .map(|(normal, _, members)| {
            let centroid = members.iter().fold(Vec3::new(0.0, 0.0, 0.0), |acc, &m| {
                acc + (points[m] - origin())
            }) * (1.0 / members.len() as f64);
            let e1 = (points[members[0]] - origin() - centroid)
                .normalized()
                .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
            let e2 = normal.cross(&e1);
            let mut ordered = members;
            ordered.sort_by(|&a, &b| {
                let pa = points[a] - origin() - centroid;
                let pb = points[b] - origin() - centroid;
                let ta = pa.dot(&e2).atan2(pa.dot(&e1));
                let tb = pb.dot(&e2).atan2(pb.dot(&e1));
                ta.partial_cmp(&tb).unwrap_or(std::cmp::Ordering::Equal)
            });
            ordered
        })
        .collect()
}

/// Uniform vertex truncation: every vertex is replaced by a flat
/// corner facet through points a fraction `t` along each incident
/// edge (`0 < t < 0.5`). Original n-gon faces become 2n-gons; a
/// valence-k vertex yields a k-gon facet.
pub fn truncate_vertices(
    points: &[Point3],
    faces: &[Vec<usize>],
    t: f64,
) -> (Vec<Point3>, Vec<Vec<usize>>) {
    let mut new_points: Vec<Point3> = Vec::new();
    let mut cut_index: HashMap<(usize, usize), usize> = HashMap::new();
    let mut cut = |a: usize, b: usize, new_points: &mut Vec<Point3>| -> usize {
        *cut_index.entry((a, b)).or_insert_with(|| {
            let p = points[a] + (points[b] - points[a]) * t;
            new_points.push(p);
            new_points.len() - 1
        })
    };

    let mut new_faces: Vec<Vec<usize>> = Vec::new();
    for face in faces {
        let n = face.len();
        let mut poly = Vec::with_capacity(2 * n);
        for i in 0..n {
            let (a, b) = (face[i], face[(i + 1) % n]);
            poly.push(cut(a, b, &mut new_points));
            poly.push(cut(b, a, &mut new_points));
        }
        new_faces.push(poly);
    }

    // Corner facets via the neighbour cycle around each vertex.
    let mut succ: HashMap<(usize, usize), usize> = HashMap::new();
    for face in faces {
        let n = face.len();
        for i in 0..n {
            let a = face[i];
            let v = face[(i + 1) % n];
            let b = face[(i + 2) % n];
            succ.insert((v, a), b);
        }
    }
    let solid_centroid = points
        .iter()
        .fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| acc + (*p - origin()))
        * (1.0 / points.len() as f64);
    for v in 0..points.len() {
        let Some((&(_, first), _)) = succ
            .iter()
            .find(|((vv, _), _)| *vv == v)
            .map(|(k, s)| (k, s))
        else {
            continue;
        };
        let mut cycle = vec![first];
        loop {
            let Some(&next) = succ.get(&(v, *cycle.last().unwrap_or(&first))) else {
                break;
            };
            if next == first {
                break;
            }
            cycle.push(next);
            if cycle.len() > points.len() {
                break;
            }
        }
        let mut poly: Vec<usize> = cycle.iter().map(|&w| cut(v, w, &mut new_points)).collect();
        // Newell normal must point away from the solid centroid.
        let mut normal = Vec3::new(0.0, 0.0, 0.0);
        for i in 0..poly.len() {
            let p = new_points[poly[i]];
            let q = new_points[poly[(i + 1) % poly.len()]];
            normal.x += (p.y - q.y) * (p.z + q.z);
            normal.y += (p.z - q.z) * (p.x + q.x);
            normal.z += (p.x - q.x) * (p.y + q.y);
        }
        let outward = points[v] - origin() - solid_centroid;
        if normal.dot(&outward) < 0.0 {
            poly.reverse();
        }
        new_faces.push(poly);
    }
    (new_points, new_faces)
}

/// Build a closed manifold BRep from vertices and face cycles.
///
/// # Errors
///
/// [`PartError::NotApplicable`] for non-manifold input (a directed
/// edge used twice, or an unpaired edge);
/// [`PartError::Core`]/validation errors from geometry attachment.
pub fn build_polyhedron_brep(
    points: &[Point3],
    faces: &[PolyFace],
    tolerance: Tolerance,
) -> Result<BRep, PartError> {
    let core = |e: popkern_core::error::CoreError| PartError::Core {
        feature: "polyhedron".into(),
        source: e,
    };
    let mut brep = BRep::new(tolerance);

    let vhs: Vec<VertexHandle> = points
        .iter()
        .map(|p| {
            brep.vertices.insert(Vertex {
                position: *p,
                outgoing: Default::default(),
                name: None,
                tolerance: None,
            })
        })
        .collect();

    let solid_h = brep.solids.insert(Solid {
        shell: Default::default(),
    });
    let shell_h = brep.shells.insert(Shell {
        face: Default::default(),
        solid: solid_h,
    });
    brep.solids[solid_h].shell = shell_h;

    let mut edge_of: HashMap<(usize, usize), EdgeHandle> = HashMap::new();
    let ensure_edges =
        |cycle: &[usize], brep: &mut BRep, edge_of: &mut HashMap<(usize, usize), EdgeHandle>| {
            let n = cycle.len();
            for i in 0..n {
                let (a, b) = (cycle[i], cycle[(i + 1) % n]);
                let key = (a.min(b), a.max(b));
                edge_of.entry(key).or_insert_with(|| {
                    brep.insert_edge_tracked(Edge {
                        half_edge: Default::default(),
                        curve: None,
                        name: None,
                        tolerance: None,
                        canonical: Default::default(),
                        role: Default::default(),
                    })
                });
            }
        };
    for face in faces {
        ensure_edges(&face.outer, &mut brep, &mut edge_of);
        for inner in &face.inners {
            ensure_edges(inner, &mut brep, &mut edge_of);
        }
    }

    let mut he_of: HashMap<(usize, usize), HalfEdgeHandle> = HashMap::new();
    for face in faces {
        let fh = brep.faces.insert(Face {
            outer_loop: Default::default(),
            inner_loops: Vec::new(),
            shell: shell_h,
            surface: None,
            name: None,
            tolerance: None,
            role: Default::default(),
        });
        let wire_cycle = |cycle: &[usize],
                              brep: &mut BRep,
                              he_of: &mut HashMap<(usize, usize), HalfEdgeHandle>|
         -> Result<LoopHandle, PartError> {
            let lh = brep.loops.insert(Loop {
                half_edge: Default::default(),
                face: fh,
            });
            let n = cycle.len();
            let mut hes: Vec<HalfEdgeHandle> = Vec::with_capacity(n);
            for i in 0..n {
                let (a, b) = (cycle[i], cycle[(i + 1) % n]);
                let eh = edge_of[&(a.min(b), a.max(b))];
                let he = brep.half_edges.insert(HalfEdge {
                    origin: vhs[a],
                    twin: Default::default(),
                    next: Default::default(),
                    prev: Default::default(),
                    edge: eh,
                    loop_: lh,
                    pcurve: None,
                });
                if he_of.insert((a, b), he).is_some() {
                    return Err(PartError::NotApplicable(format!(
                        "polyhedron is not manifold: directed edge ({a},{b}) used twice"
                    )));
                }
                brep.edges[eh].half_edge = he;
                brep.vertices[vhs[a]].outgoing = he;
                hes.push(he);
            }
            for i in 0..n {
                brep.half_edges[hes[i]].next = hes[(i + 1) % n];
                brep.half_edges[hes[i]].prev = hes[(i + n - 1) % n];
            }
            brep.loops[lh].half_edge = hes[0];
            Ok(lh)
        };
        let outer_lh = wire_cycle(&face.outer, &mut brep, &mut he_of)?;
        brep.faces[fh].outer_loop = outer_lh;
        for inner in &face.inners {
            let lh = wire_cycle(inner, &mut brep, &mut he_of)?;
            brep.faces[fh].inner_loops.push(lh);
        }
    }
    brep.shells[shell_h].face = brep.faces.keys().next().unwrap_or_default();

    for (&(a, b), &he) in &he_of {
        let twin = *he_of.get(&(b, a)).ok_or_else(|| {
            PartError::NotApplicable(format!(
                "polyhedron is not closed: edge ({a},{b}) has no opposite face use"
            ))
        })?;
        brep.half_edges[he].twin = twin;
    }

    // Analytic line curves.
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for eh in edge_handles {
        let he = brep.edges[eh].half_edge;
        let a = brep.vertices[brep.half_edges[he].origin].position;
        let b = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
        let dir = b - a;
        let len = dir.length();
        if len < 1e-12 {
            return Err(PartError::NotApplicable(
                "polyhedron has a zero-length edge".into(),
            ));
        }
        brep.edges[eh].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: a,
                direction: dir * (1.0 / len),
            },
            t_start: 0.0,
            t_end: len,
        });
    }
    // Right-handed plane charts from Newell normals.
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    for fh in face_handles {
        let lh = brep.faces[fh].outer_loop;
        let start = brep.loops[lh].half_edge;
        let mut normal = Vec3::new(0.0, 0.0, 0.0);
        let mut face_origin = None;
        let mut he = start;
        loop {
            let d = &brep.half_edges[he];
            let p = brep.vertices[d.origin].position;
            let q = brep.vertices[brep.half_edges[d.next].origin].position;
            normal.x += (p.y - q.y) * (p.z + q.z);
            normal.y += (p.z - q.z) * (p.x + q.x);
            normal.z += (p.x - q.x) * (p.y + q.y);
            face_origin.get_or_insert(p);
            he = d.next;
            if he == start {
                break;
            }
        }
        let normal = normal.normalized().ok_or_else(|| {
            PartError::NotApplicable("polyhedron face has a degenerate normal".into())
        })?;
        let face_origin = face_origin.unwrap_or(origin());
        let reference = if normal.x.abs() <= normal.y.abs() && normal.x.abs() <= normal.z.abs() {
            Vec3::new(1.0, 0.0, 0.0)
        } else if normal.y.abs() <= normal.z.abs() {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let u_axis = (reference - normal * reference.dot(&normal))
            .normalized()
            .unwrap_or(Vec3::new(1.0, 0.0, 0.0));
        let v_axis = normal.cross(&u_axis);
        brep.faces[fh].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: face_origin,
                normal,
                u_axis,
                v_axis,
            },
            orientation: true,
        });
    }
    brep.recompute_all_pcurves().map_err(core)?;
    brep.validate().map_err(|source| PartError::Validation {
        feature: "polyhedron".into(),
        source,
    })?;
    Ok(brep)
}

/// Recover `(points, faces)` from a polyhedral BRep (every face
/// planar with straight edges). The inverse of
/// [`build_polyhedron_brep`], used so truncation and engraving can
/// operate on any polyhedral [`Part`](crate::Part), not just ones
/// built by this module.
///
/// # Errors
///
/// [`PartError::NotApplicable`] naming the offending face when a
/// surface is not a plane or an edge is not a line.
pub fn extract_polyhedron(brep: &BRep) -> Result<(Vec<Point3>, Vec<PolyFace>), PartError> {
    let mut index_of: HashMap<VertexHandle, usize> = HashMap::new();
    let mut points: Vec<Point3> = Vec::new();
    for (vh, v) in brep.vertices.iter() {
        index_of.insert(vh, points.len());
        points.push(v.position);
    }
    let mut faces = Vec::new();
    for (fh, face) in brep.faces.iter() {
        match face.surface.as_ref().map(|s| &s.geometry) {
            Some(SurfaceGeometry::Plane { .. }) => {}
            other => {
                return Err(PartError::NotApplicable(format!(
                    "face {fh:?} is not planar ({other:?}): polyhedral operations require planar faces"
                )));
            }
        }
        let cycle = |lh: LoopHandle| -> Result<Vec<usize>, PartError> {
            let start = brep
                .loops
                .get(lh)
                .map(|l| l.half_edge)
                .ok_or_else(|| PartError::NotApplicable("missing loop".into()))?;
            let mut out = Vec::new();
            let mut he = start;
            loop {
                let d = brep
                    .half_edges
                    .get(he)
                    .ok_or_else(|| PartError::NotApplicable("broken loop".into()))?;
                if let Some(c) = brep.edges.get(d.edge).and_then(|e| e.curve.as_ref())
                    && !matches!(c.geometry, CurveGeometry::Line { .. })
                {
                    return Err(PartError::NotApplicable(format!(
                        "edge {:?} is not a straight line: polyhedral operations require line edges",
                        d.edge
                    )));
                }
                out.push(index_of[&d.origin]);
                he = d.next;
                if he == start {
                    break;
                }
            }
            Ok(out)
        };
        faces.push(PolyFace {
            outer: cycle(face.outer_loop)?,
            inners: face
                .inner_loops
                .iter()
                .map(|&lh| cycle(lh))
                .collect::<Result<Vec<_>, _>>()?,
        });
    }
    Ok((points, faces))
}
