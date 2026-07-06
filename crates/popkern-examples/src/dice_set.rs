//! Polyhedral dice set: d4, d6, d8, d10, d12, d20 with truncated
//! (safety-trimmed) corners and ENGRAVED face numbers, exported as
//! watertight STEP solids.
//!
//! Each die is an exact polyhedron. Corner trimming is uniform
//! VERTEX TRUNCATION computed directly in the geometry — every
//! vertex is replaced by a flat corner facet whose points sit a
//! fixed fraction along each incident edge. Numbers are engraved as
//! 7-segment digits, each segment an independent rectangular-hexagon
//! POCKET: an inner loop on the face, four/six straight walls, and a
//! flat floor. The result is still a polyhedron with planar faces
//! and analytic line edges: no boolean or blend machinery, exact by
//! construction.
//!
//! Pipeline per die: canonical vertex coordinates → convex hull with
//! coplanar-face merging (one algorithm covers the Platonic solids
//! AND the d10 pentagonal trapezohedron's kite faces) → vertex
//! truncation → 7-segment number engraving on the original faces →
//! generic polyhedron BRep wiring (inner-loop aware) → validate →
//! watertight tessellation gate → STEP export (per die + full set).

use std::collections::HashMap;

use popkern_core::geometry::curve::{CurveGeometry, CurveOnEdge};
use popkern_core::geometry::surface::SurfaceGeometry;
use popkern_core::geometry::{Point3, SurfaceOnFace, Vec3};
use popkern_core::tolerance::Tolerance;
use popkern_core::topology::*;

/// Geometric epsilon for hull plane membership on unit-scale
/// canonical coordinates (exact symmetric constants; the epsilon
/// only absorbs f64 arithmetic noise).
const HULL_EPS: f64 = 1e-9;

/// A polyhedron face: outer boundary plus engraved inner loops, all
/// as vertex-index cycles. Outer cycles wind CCW viewed from
/// outside; inner cycles CW (material on the left of the walk).
struct FaceSpec {
    outer: Vec<usize>,
    inners: Vec<Vec<usize>>,
}

impl FaceSpec {
    fn simple(outer: Vec<usize>) -> Self {
        FaceSpec {
            outer,
            inners: Vec::new(),
        }
    }
}

fn origin() -> Point3 {
    Point3::new(0.0, 0.0, 0.0)
}

// ───────────────────────── convex hull ─────────────────────────

/// Faces of the convex hull of `points` (assumed in convex position,
/// as all dice coordinates are), with coplanar triangles merged into
/// polygons. Each face is a vertex-index cycle wound CCW viewed from
/// OUTSIDE. O(n⁴) plane enumeration — fine for n ≤ 32.
fn convex_hull_faces(points: &[Point3]) -> Vec<Vec<usize>> {
    let n = points.len();
    let mut faces: Vec<(Vec3, f64, Vec<usize>)> = Vec::new(); // (unit normal, offset, members)
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let normal = (points[j] - points[i]).cross(&(points[k] - points[i]));
                let Some(mut normal) = normal.normalized() else {
                    continue;
                };
                let mut offset = normal.dot(&(points[i] - origin()));
                // Orient outward (all other points on the negative side).
                let mut lo = f64::INFINITY;
                let mut hi = f64::NEG_INFINITY;
                for p in points {
                    let d = normal.dot(&(*p - origin())) - offset;
                    lo = lo.min(d);
                    hi = hi.max(d);
                }
                if hi > HULL_EPS && lo < -HULL_EPS {
                    continue; // points on both sides: not a hull plane
                }
                if hi > HULL_EPS {
                    normal = normal * -1.0;
                    offset = -offset;
                }
                if faces.iter().any(|(fn_, fo, _)| {
                    fn_.dot(&normal) > 1.0 - 1e-7 && (fo - offset).abs() < HULL_EPS
                }) {
                    continue; // plane already found
                }
                let members: Vec<usize> = (0..n)
                    .filter(|&m| (normal.dot(&(points[m] - origin())) - offset).abs() < HULL_EPS)
                    .collect();
                faces.push((normal, offset, members));
            }
        }
    }
    // Order each face's members CCW around the outward normal.
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

// ─────────────────────── vertex truncation ───────────────────────

/// Uniform vertex truncation: each vertex is cut by a corner facet
/// through points a fraction `t` along every incident edge. Original
/// n-gon faces become 2n-gons; each vertex of valence k yields a
/// k-gon corner facet. `t < 0.5` keeps the cut planes disjoint.
fn truncate_vertices(
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
    // Original faces → 2n-gons: on each edge (a→b) keep the segment
    // between the two cut points.
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

    // Corner facets: the neighbours of v, cycled via "in the face
    // [..., a, v, b, ...] the successor of a around v is b".
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
            let &next = succ.get(&(v, *cycle.last().unwrap())).expect("closed fan");
            if next == first {
                break;
            }
            cycle.push(next);
            assert!(cycle.len() <= points.len(), "corner fan did not close");
        }
        let mut poly: Vec<usize> = cycle.iter().map(|&w| cut(v, w, &mut new_points)).collect();
        // Ensure the corner facet winds CCW from outside (Newell
        // normal must point away from the solid centroid).
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

// ─────────────────── 7-segment digit engraving ───────────────────

/// The seven segments in a unit digit box (width 0.62, height 1.0),
/// as CCW polygons: hexagonal bars with 45° pointed ends and small
/// tip gaps so segments never touch (each is an independent pocket).
fn segment_polygon(horizontal: bool, cx: f64, cy: f64, half_len: f64, t: f64) -> Vec<(f64, f64)> {
    let h = t / 2.0;
    let raw: Vec<(f64, f64)> = if horizontal {
        vec![
            (cx - half_len, cy),
            (cx - half_len + h, cy - h),
            (cx + half_len - h, cy - h),
            (cx + half_len, cy),
            (cx + half_len - h, cy + h),
            (cx - half_len + h, cy + h),
        ]
    } else {
        vec![
            (cx, cy - half_len),
            (cx + h, cy - half_len + h),
            (cx + h, cy + half_len - h),
            (cx, cy + half_len),
            (cx - h, cy + half_len - h),
            (cx - h, cy - half_len + h),
        ]
    };
    // Enforce CCW via the shoelace sign.
    let area2: f64 = raw
        .iter()
        .zip(raw.iter().cycle().skip(1))
        .map(|(a, b)| a.0 * b.1 - b.0 * a.1)
        .take(raw.len())
        .sum();
    if area2 < 0.0 {
        raw.into_iter().rev().collect()
    } else {
        raw
    }
}

/// Segment polygons for one digit, in a box x ∈ [0, 0.62], y ∈ [0, 1].
fn digit_polygons(d: u8) -> Vec<Vec<(f64, f64)>> {
    const W: f64 = 0.62;
    const T: f64 = 0.15;
    const GAP: f64 = 0.028;
    let hl_h = (W - T) / 2.0 - GAP; // horizontal bar half-length
    let hl_v = 0.25 - GAP - T / 2.0; // vertical bar half-length
    let seg = |name: char| -> Vec<(f64, f64)> {
        match name {
            'A' => segment_polygon(true, W / 2.0, 1.0 - T / 2.0, hl_h, T),
            'G' => segment_polygon(true, W / 2.0, 0.5, hl_h, T),
            'D' => segment_polygon(true, W / 2.0, T / 2.0, hl_h, T),
            'B' => segment_polygon(false, W - T / 2.0, 0.75, hl_v, T),
            'C' => segment_polygon(false, W - T / 2.0, 0.25, hl_v, T),
            'F' => segment_polygon(false, T / 2.0, 0.75, hl_v, T),
            'E' => segment_polygon(false, T / 2.0, 0.25, hl_v, T),
            _ => unreachable!(),
        }
    };
    let segments: &str = match d {
        0 => "ABCDEF",
        1 => "BC",
        2 => "ABGED",
        3 => "ABGCD",
        4 => "FGBC",
        5 => "AFGCD",
        6 => "AFGECD",
        7 => "ABC",
        8 => "ABCDEFG",
        9 => "ABCDFG",
        _ => unreachable!(),
    };
    segments.chars().map(seg).collect()
}

/// Engrave `number` into face `face_idx` as 7-segment pockets of the
/// given `depth`. Each segment polygon becomes: an inner loop on the
/// face (reversed winding), one wall quad per polygon edge, and a
/// flat floor face — manifold by construction.
fn engrave_number(
    points: &mut Vec<Point3>,
    faces: &mut Vec<FaceSpec>,
    face_idx: usize,
    number: usize,
    depth: f64,
) {
    let outline: Vec<Point3> = faces[face_idx].outer.iter().map(|&i| points[i]).collect();
    let n = outline.len();
    let centroid = outline
        .iter()
        .fold(Vec3::new(0.0, 0.0, 0.0), |acc, p| acc + (*p - origin()))
        * (1.0 / n as f64);
    let centroid = origin() + centroid;
    // Face frame: outward normal (Newell), text-up = the projection
    // of global +Z (or +Y for horizontal faces), text-right = up × n
    // so digits read correctly from outside.
    let mut normal = Vec3::new(0.0, 0.0, 0.0);
    for i in 0..n {
        let p = outline[i];
        let q = outline[(i + 1) % n];
        normal.x += (p.y - q.y) * (p.z + q.z);
        normal.y += (p.z - q.z) * (p.x + q.x);
        normal.z += (p.x - q.x) * (p.y + q.y);
    }
    let normal = normal.normalized().expect("face normal");
    let up_hint = if normal.z.abs() > 0.9 {
        Vec3::new(0.0, 1.0, 0.0)
    } else {
        Vec3::new(0.0, 0.0, 1.0)
    };
    let up = (up_hint - normal * up_hint.dot(&normal))
        .normalized()
        .expect("in-plane up");
    let right = up.cross(&normal);

    // Inradius: nearest outline edge to the centroid.
    let mut inradius = f64::INFINITY;
    for i in 0..n {
        let a = outline[i];
        let b = outline[(i + 1) % n];
        let ab = b - a;
        let len = ab.length();
        let t = ((centroid - a).dot(&ab) / (len * len)).clamp(0.0, 1.0);
        let foot = a + ab * t;
        inradius = inradius.min(centroid.distance_to(&foot));
    }

    // Layout the digits, centred on the actual glyph bounding box.
    let digits: Vec<u8> = number.to_string().bytes().map(|b| b - b'0').collect();
    const DIGIT_ADVANCE: f64 = 0.62 + 0.18;
    let mut polys: Vec<Vec<(f64, f64)>> = Vec::new();
    for (i, &d) in digits.iter().enumerate() {
        let dx = i as f64 * DIGIT_ADVANCE;
        for poly in digit_polygons(d) {
            polys.push(poly.into_iter().map(|(x, y)| (x + dx, y)).collect());
        }
    }
    let (mut min_x, mut max_x) = (f64::INFINITY, f64::NEG_INFINITY);
    let (mut min_y, mut max_y) = (f64::INFINITY, f64::NEG_INFINITY);
    for poly in &polys {
        for &(x, y) in poly {
            min_x = min_x.min(x);
            max_x = max_x.max(x);
            min_y = min_y.min(y);
            max_y = max_y.max(y);
        }
    }
    let (bw, bh) = (max_x - min_x, max_y - min_y);
    let half_diag = ((bw / 2.0).powi(2) + (bh / 2.0).powi(2)).sqrt();
    let scale = 0.66 * inradius / half_diag;
    let (cx, cy) = ((min_x + max_x) / 2.0, (min_y + max_y) / 2.0);

    for poly in polys {
        // Top ring (on the face plane) and floor ring vertex indices.
        let top: Vec<usize> = poly
            .iter()
            .map(|&(x, y)| {
                let p = centroid + right * ((x - cx) * scale) + up * ((y - cy) * scale);
                points.push(p);
                points.len() - 1
            })
            .collect();
        let bottom: Vec<usize> = top
            .iter()
            .map(|&i| {
                let p = points[i] + normal * (-depth);
                points.push(p);
                points.len() - 1
            })
            .collect();
        let m = top.len();
        // Inner loop on the engraved face: the pocket outline
        // reversed (CW from outside — material on the left).
        let mut inner = top.clone();
        inner.reverse();
        faces[face_idx].inners.push(inner);
        // Walls: for each CCW outline edge (a→b), the wall cycle
        // [a, b, b_floor, a_floor] uses the directed edge (a→b) so
        // it twins the inner loop's (b→a), and (b_floor→a_floor)
        // twins the floor's (a_floor→b_floor).
        for i in 0..m {
            let (a, b) = (top[i], top[(i + 1) % m]);
            let (af, bf) = (bottom[i], bottom[(i + 1) % m]);
            faces.push(FaceSpec::simple(vec![a, b, bf, af]));
        }
        // Floor: the outline at depth, CCW (its outward normal
        // points out of the material, into the pocket).
        faces.push(FaceSpec::simple(bottom));
    }
}

// ─────────────────── generic polyhedron wiring ───────────────────

/// Build a closed manifold BRep from a polyhedron's vertices and
/// face cycles (outer CCW-from-outside, inners CW): half-edge
/// wiring, analytic line curves, right-handed plane charts
/// (v_axis = normal × u_axis, so STEP export needs no sense
/// compensation), and planar pcurves on every loop.
fn build_polyhedron(
    points: &[Point3],
    faces: &[FaceSpec],
    tolerance: Tolerance,
) -> Result<BRep, popkern_core::error::CoreError> {
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

    // Edges: one per unordered vertex pair present in any cycle.
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

    // Faces, loops, half-edges.
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
         -> LoopHandle {
            let lh = brep.loops.insert(Loop {
                half_edge: Default::default(),
                face: fh,
            });
            let n = cycle.len();
            let hes: Vec<HalfEdgeHandle> = (0..n)
                .map(|i| {
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
                    assert!(
                        he_of.insert((a, b), he).is_none(),
                        "directed edge ({a},{b}) used twice — input not manifold"
                    );
                    brep.edges[eh].half_edge = he;
                    brep.vertices[vhs[a]].outgoing = he;
                    he
                })
                .collect();
            for i in 0..n {
                brep.half_edges[hes[i]].next = hes[(i + 1) % n];
                brep.half_edges[hes[i]].prev = hes[(i + n - 1) % n];
            }
            brep.loops[lh].half_edge = hes[0];
            lh
        };
        let outer_lh = wire_cycle(&face.outer, &mut brep, &mut he_of);
        brep.faces[fh].outer_loop = outer_lh;
        for inner in &face.inners {
            let lh = wire_cycle(inner, &mut brep, &mut he_of);
            brep.faces[fh].inner_loops.push(lh);
        }
    }
    brep.shells[shell_h].face = brep.faces.keys().next().unwrap_or_default();

    // Twin wiring: (a,b) pairs with (b,a); a miss means the input
    // was not a closed manifold.
    for (&(a, b), &he) in &he_of {
        let twin = *he_of.get(&(b, a)).ok_or_else(|| {
            popkern_core::error::CoreError::PreconditionFailed(format!(
                "polyhedron not closed: edge ({a},{b}) has no opposite face use"
            ))
        })?;
        brep.half_edges[he].twin = twin;
    }

    // Geometry: line curves; right-handed plane charts.
    let edge_handles: Vec<EdgeHandle> = brep.edges.keys().collect();
    for eh in edge_handles {
        let he = brep.edges[eh].half_edge;
        let a = brep.vertices[brep.half_edges[he].origin].position;
        let b = brep.vertices[brep.half_edges[brep.half_edges[he].twin].origin].position;
        let dir = b - a;
        let len = dir.length();
        brep.edges[eh].curve = Some(CurveOnEdge {
            geometry: CurveGeometry::Line {
                origin: a,
                direction: dir * (1.0 / len),
            },
            t_start: 0.0,
            t_end: len,
        });
    }
    let face_handles: Vec<FaceHandle> = brep.faces.keys().collect();
    for fh in face_handles {
        let lh = brep.faces[fh].outer_loop;
        let start = brep.loops[lh].half_edge;
        // Newell normal over the outer loop (robust for any polygon).
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
            popkern_core::error::CoreError::PreconditionFailed("degenerate face normal".into())
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
    brep.recompute_all_pcurves()?;
    brep.validate()
        .map_err(popkern_core::error::CoreError::Validation)?;
    Ok(brep)
}

// ───────────────────────── dice data ─────────────────────────

/// Canonical (unit-scale) vertex coordinates for one die.
fn die_vertices(sides: usize) -> Vec<Point3> {
    let phi = (1.0 + 5.0_f64.sqrt()) / 2.0;
    match sides {
        4 => vec![
            Point3::new(1.0, 1.0, 1.0),
            Point3::new(1.0, -1.0, -1.0),
            Point3::new(-1.0, 1.0, -1.0),
            Point3::new(-1.0, -1.0, 1.0),
        ],
        6 => {
            let mut v = Vec::new();
            for &x in &[-1.0, 1.0] {
                for &y in &[-1.0, 1.0] {
                    for &z in &[-1.0, 1.0] {
                        v.push(Point3::new(x, y, z));
                    }
                }
            }
            v
        }
        8 => vec![
            Point3::new(1.0, 0.0, 0.0),
            Point3::new(-1.0, 0.0, 0.0),
            Point3::new(0.0, 1.0, 0.0),
            Point3::new(0.0, -1.0, 0.0),
            Point3::new(0.0, 0.0, 1.0),
            Point3::new(0.0, 0.0, -1.0),
        ],
        10 => {
            // Pentagonal trapezohedron: alternating zigzag ring plus
            // two apexes placed on the kite planes (each kite =
            // [apex, upper, lower, upper], planar by construction:
            // the apex is the z-axis intercept of the plane through
            // one upper-lower-upper triple).
            let h = 0.12;
            let mut v: Vec<Point3> = (0..10)
                .map(|k| {
                    let th = std::f64::consts::TAU * k as f64 / 10.0;
                    let z = if k % 2 == 0 { h } else { -h };
                    Point3::new(th.cos(), th.sin(), z)
                })
                .collect();
            let (a, b, c) = (v[0], v[1], v[2]);
            let n = (b - a).cross(&(c - a));
            // Plane n·(p − a) = 0 crossed with the z-axis.
            let z0 = (n.dot(&(a - origin()))) / n.z;
            v.push(Point3::new(0.0, 0.0, z0.abs()));
            v.push(Point3::new(0.0, 0.0, -z0.abs()));
            v
        }
        12 => {
            let mut v = Vec::new();
            for &x in &[-1.0, 1.0] {
                for &y in &[-1.0, 1.0] {
                    for &z in &[-1.0, 1.0] {
                        v.push(Point3::new(x, y, z));
                    }
                }
            }
            let inv = 1.0 / phi;
            for &(p, q) in &[(inv, phi)] {
                for &s1 in &[-1.0, 1.0] {
                    for &s2 in &[-1.0, 1.0] {
                        v.push(Point3::new(0.0, s1 * p, s2 * q));
                        v.push(Point3::new(s1 * p, s2 * q, 0.0));
                        v.push(Point3::new(s1 * q, 0.0, s2 * p));
                    }
                }
            }
            v
        }
        20 => {
            let mut v = Vec::new();
            for &s1 in &[-1.0, 1.0] {
                for &s2 in &[-1.0, 1.0] {
                    v.push(Point3::new(0.0, s1, s2 * phi));
                    v.push(Point3::new(s1, s2 * phi, 0.0));
                    v.push(Point3::new(s2 * phi, 0.0, s1));
                }
            }
            v
        }
        other => panic!("no die with {other} sides"),
    }
}

fn main() {
    let tol = Tolerance::default();
    let dice = [4usize, 6, 8, 10, 12, 20];
    // Corner-trim fraction: how far along each edge the cut points
    // sit. 0.18 reads as a generous, foot-friendly trim.
    let trim = 0.18;
    // Uniform circumradius (model units ~ mm), row spacing, and
    // engraving depth.
    let radius = 10.0;
    let spacing = 26.0;
    let depth = 0.5;

    let out_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../examples/output");
    std::fs::create_dir_all(&out_dir).expect("out dir");

    let mut products = Vec::new();
    let mut merged_obj = String::new();
    let mut obj_vertex_base = 1usize;

    for (i, &sides) in dice.iter().enumerate() {
        let raw = die_vertices(sides);
        // Normalize to the common circumradius.
        let max_r = raw
            .iter()
            .map(|p| (*p - origin()).length())
            .fold(0.0_f64, f64::max);
        let scaled: Vec<Point3> = raw
            .iter()
            .map(|p| {
                Point3::new(
                    p.x / max_r * radius,
                    p.y / max_r * radius,
                    p.z / max_r * radius,
                )
            })
            .collect();
        let faces = convex_hull_faces(&scaled);
        assert_eq!(
            faces.len(),
            sides,
            "d{sides}: hull must recover exactly {sides} faces, got {}",
            faces.len()
        );
        let (mut tp, tf) = truncate_vertices(&scaled, &faces, trim);
        let mut specs: Vec<FaceSpec> = tf.into_iter().map(FaceSpec::simple).collect();
        // Engrave the numbers on the original (first `sides`) faces:
        // d10 runs 0–9 by tradition, everything else 1–N.
        for f in 0..sides {
            let number = if sides == 10 { f } else { f + 1 };
            engrave_number(&mut tp, &mut specs, f, number, depth);
        }
        // Place along the row.
        let dx = i as f64 * spacing;
        for p in &mut tp {
            p.x += dx;
        }
        let brep = build_polyhedron(&tp, &specs, tol)
            .unwrap_or_else(|e| panic!("d{sides} build failed: {e:?}"));
        let mesh = popkern_mesh::tessellate(&brep, &popkern_mesh::TessellationOptions::default())
            .unwrap_or_else(|e| panic!("d{sides} tessellation failed: {e:?}"));
        let open = open_edge_count(&mesh);
        println!(
            "d{sides}: {} faces, {} verts, mesh {} tris, open edges {open}",
            specs.len(),
            tp.len(),
            mesh.triangles.len()
        );
        assert_eq!(open, 0, "d{sides} must tessellate watertight");

        // Merge into the combined OBJ for rendering.
        for v in &mesh.vertices {
            merged_obj.push_str(&format!(
                "v {} {} {}\n",
                v.position.x, v.position.y, v.position.z
            ));
        }
        for t in &mesh.triangles {
            merged_obj.push_str(&format!(
                "f {} {} {}\n",
                t[0] + obj_vertex_base,
                t[1] + obj_vertex_base,
                t[2] + obj_vertex_base
            ));
        }
        obj_vertex_base += mesh.vertices.len();

        // Individual STEP.
        let doc = single_product_doc(&brep, &format!("d{sides}"));
        let path = out_dir.join(format!("die_d{sides}.step"));
        popkern_io::export_step_file(&doc, &path, &popkern_io::StepExportOptions::default())
            .unwrap_or_else(|e| panic!("d{sides} STEP export failed: {e}"));

        products.push(product_node(&brep, &format!("d{sides}")));
    }

    // Combined set.
    let set_doc = popkern_io::StepDocument {
        roots: (0..products.len()).collect(),
        products,
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    };
    let set_path = out_dir.join("dice_set.step");
    popkern_io::export_step_file(
        &set_doc,
        &set_path,
        &popkern_io::StepExportOptions::default(),
    )
    .expect("dice set STEP export");
    println!("wrote {}", set_path.display());

    let obj_path = out_dir.join("dice_set.obj");
    std::fs::write(&obj_path, merged_obj).expect("obj write");
    println!("wrote {}", obj_path.display());
}

fn product_node(brep: &BRep, name: &str) -> popkern_io::StepProductNode {
    popkern_io::StepProductNode {
        name: name.into(),
        id: format!("DICE-{}", name.to_uppercase()),
        children: vec![],
        shapes: vec![std::sync::Arc::new(popkern_io::ShapeRep {
            kind: popkern_io::ShapeRepKind::AdvancedBrep,
            solids: vec![brep.clone()],
            context: popkern_io::RepresentationContext::default(),
            non_brep_items: vec![],
        })],
        metadata: popkern_io::ProductMetadata::default(),
    }
}

fn single_product_doc(brep: &BRep, name: &str) -> popkern_io::StepDocument {
    popkern_io::StepDocument {
        products: vec![product_node(brep, name)],
        roots: vec![0],
        units: popkern_io::UnitContext::default(),
        opaque_bag: popkern_io::OpaqueBag::default(),
    }
}

/// Welded open-edge count (the watertightness gate — same census the
/// flagship uses).
fn open_edge_count(mesh: &popkern_mesh::TriMesh) -> usize {
    let mut key: HashMap<(i64, i64, i64), usize> = HashMap::new();
    let id = |p: Point3, key: &mut HashMap<(i64, i64, i64), usize>| -> usize {
        let k = (
            (p.x * 1.0e7).round() as i64,
            (p.y * 1.0e7).round() as i64,
            (p.z * 1.0e7).round() as i64,
        );
        let next = key.len();
        *key.entry(k).or_insert(next)
    };
    let mut counts: HashMap<(usize, usize), usize> = HashMap::new();
    for t in &mesh.triangles {
        let w = [
            id(mesh.vertices[t[0]].position, &mut key),
            id(mesh.vertices[t[1]].position, &mut key),
            id(mesh.vertices[t[2]].position, &mut key),
        ];
        if w[0] == w[1] || w[1] == w[2] || w[0] == w[2] {
            continue;
        }
        for (a, b) in [(w[0], w[1]), (w[1], w[2]), (w[2], w[0])] {
            *counts.entry((a.min(b), a.max(b))).or_insert(0) += 1;
        }
    }
    counts.values().filter(|&&c| c == 1).count()
}
