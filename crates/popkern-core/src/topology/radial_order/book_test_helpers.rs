//! Test fixtures for radial order: a book of planar pages hinged on
//! one spine edge.

use crate::constructors::BRepBuilder;
use crate::geometry::surface::{SurfaceGeometry, SurfaceOnFace};
use crate::geometry::{Point3, Vec3};
use crate::tolerance::Tolerance;
use crate::topology::{BRep, EdgeHandle, HalfEdge, HalfEdgeHandle};
use popkern_math::PortableMath;

/// One degree in radians.
pub(crate) const DEG: f64 = std::f64::consts::PI / 180.0;

/// The default tolerance.
pub(crate) fn tol() -> Tolerance {
    Tolerance::default()
}

/// A book: square pages hinged on the spine x in [0, 1], page k
/// leaving it toward (0, cos θ_k, sin θ_k). Page k's spine fin runs
/// the spine forward when `forward[k]`. The ring is linked in the
/// order the pages are given. Returns the body, the spine and each
/// page's spine fin.
pub(crate) fn book(angles: &[f64], forward: &[bool]) -> (BRep, EdgeHandle, Vec<HalfEdgeHandle>) {
    let mut b = BRepBuilder::new(tol());
    let solid = b.add_solid();
    let shell = b.add_shell(solid);
    let s0 = b.add_vertex(Point3::new(0.0, 0.0, 0.0));
    let s1 = b.add_vertex(Point3::new(1.0, 0.0, 0.0));
    let (spine, along, against) = b.add_edge(s0, s1);
    let x = Vec3::new(1.0, 0.0, 0.0);
    let mut fins = Vec::new();
    for (&theta, &fwd) in angles.iter().zip(forward) {
        let w = Vec3::new(0.0, theta.pcos(), theta.psin());
        let a = b.add_vertex(Point3::new(0.0, 0.0, 0.0) + w);
        let c = b.add_vertex(Point3::new(1.0, 0.0, 0.0) + w);
        // The first page each way takes the spine's own fin that way,
        // the rest copies of it.
        let model = if fwd { along } else { against };
        let fin = if fins.contains(&model) {
            let copy: HalfEdge = b.brep().half_edges[model].clone();
            b.brep_mut().half_edges.insert(copy)
        } else {
            model
        };
        let face = b.add_face(shell);
        let loop_fins = if fwd {
            let (_, e1, _) = b.add_edge(s1, c);
            let (_, e2, _) = b.add_edge(c, a);
            let (_, e3, _) = b.add_edge(a, s0);
            [fin, e1, e2, e3]
        } else {
            let (_, e1, _) = b.add_edge(s0, a);
            let (_, e2, _) = b.add_edge(a, c);
            let (_, e3, _) = b.add_edge(c, s1);
            [fin, e1, e2, e3]
        };
        b.add_loop(&loop_fins, face).expect("loop");
        let (u_axis, v_axis) = if fwd { (x, w) } else { (w, x) };
        b.brep_mut().faces[face].surface = Some(SurfaceOnFace {
            geometry: SurfaceGeometry::Plane {
                origin: Point3::new(0.0, 0.0, 0.0),
                normal: u_axis.cross(&v_axis),
                u_axis,
                v_axis,
            },
            orientation: true,
        });
        fins.push(fin);
    }
    let mut brep = b.build_unchecked();
    brep.link_ring(&fins);
    for &f in &fins {
        brep.half_edges[f].edge = spine;
    }
    (brep, spine, fins)
}
