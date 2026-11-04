//! FaceInfo and resolution records (plan §4.3, §4.3a output).

use popkern_core::geometry::Point3;
use popkern_core::topology::{FaceHandle, VertexHandle};

use super::candidates::{CandidateVertex, Incidence, SourceKind};
use super::consequences::Collapse;
use super::face_info::face_infos;
use super::pave_blocks::{BlockOwner, pave_blocks};
use super::resolution_record::{VertexConsequence, resolution_json, vertex_records};
use super::resolve::{Incompatibility, ResolutionContext, resolve};
use super::tests::{TOL, cube, line, options, record, record_end, vertex};
use slotmap::Key;

/// A face of the cube and one of its boundary vertices.
fn face_and_corner(c: &popkern_core::topology::BRep) -> (FaceHandle, VertexHandle) {
    let (fh, f) = c.faces.iter().next().expect("face");
    let he = c.loops[f.outer_loop].half_edge;
    (fh, c.half_edges[he].origin)
}

/// A record on (face, face) with an end at the face's corner and a
/// point inside the face: the corner is on the face, the point in it,
/// and the record's pieces are its section blocks, on both operands'
/// faces; a trim residue's pieces only on its owner's.
#[test]
fn face_info_sorts_vertices_on_and_in_and_section_blocks_per_face() {
    let c = cube();
    let (face, corner) = face_and_corner(&c);
    let at = c.vertices[corner].position;
    let far = Point3::new(at.x + 0.5, at.y + 0.5, at.z + 0.5);
    for (residue, faces) in [(None, 2), (Some(true), 1)] {
        let mut r = record(line(at, far));
        r.face_a = face;
        r.face_b = face;
        r.residue_owner_a = residue;
        let records = [r];
        let inside = CandidateVertex {
            position: Point3::new(at.x + 0.25, at.y + 0.25, at.z + 0.25),
            radius: TOL,
            zone: TOL,
            source: SourceKind::EdgeFaceHit,
            incidence: vec![
                Incidence::Face {
                    operand: 0,
                    face,
                    uv: None,
                },
                Incidence::Record {
                    record: 0,
                    s: 0.5,
                    end: None,
                },
            ],
            key: (SourceKind::EdgeFaceHit, 0, 0, 0, 0),
        };
        let cands = vec![
            vertex(0, corner, at),
            record_end(0, 0.0, false, at),
            inside,
            record_end(0, 1.0, true, far),
        ];
        let ctx = ResolutionContext {
            operands: [&c, &c],
            records: &records,
        };
        let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
        let blocks = pave_blocks(&cands, &res, [&c, &c], &records, TOL);
        let infos = face_infos(&cands, &res, [&c, &c], &records, &blocks);
        assert_eq!(infos.len(), faces, "{infos:?}");
        let info = &infos[&(0, face.data().as_ffi())];
        assert!(info.on_vertices.contains(&res.group_of[0]), "{info:?}");
        assert_eq!(
            info.in_vertices.iter().copied().collect::<Vec<_>>(),
            [res.group_of[2]],
            "{info:?}"
        );
        assert_eq!(info.section_blocks.len(), 2, "{info:?}");
        assert!(
            info.section_blocks
                .iter()
                .all(|&k| blocks[k].owner == BlockOwner::Record(0))
        );
    }
}

/// Every vertex's record names its members, the rule behind each split
/// and its consequences: two vertices of one operand kept apart on a
/// record, the piece between them a micro block, and a record whose
/// ends met within a vertex collapsing into it.
#[test]
fn a_vertex_record_names_its_splits_and_consequences() {
    let c = cube();
    let vs: Vec<VertexHandle> = c.vertices.keys().take(2).collect();
    let records = [
        record(line(Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0))),
        record(line(
            Point3::new(5.0, 0.0, 0.0),
            Point3::new(5.0 + 0.1 * TOL, 0.0, 0.0),
        )),
    ];
    let on_record = |mut cand: CandidateVertex, s: f64| {
        cand.incidence.push(Incidence::Record {
            record: 0,
            s,
            end: None,
        });
        cand
    };
    let cands = vec![
        on_record(vertex(0, vs[0], Point3::new(0.5, 0.0, 0.0)), 0.5),
        on_record(
            vertex(0, vs[1], Point3::new(0.5 + 0.8 * TOL, 0.0, 0.0)),
            0.5 + 0.8 * TOL,
        ),
        record_end(1, 0.0, false, Point3::new(5.0, 0.0, 0.0)),
        record_end(1, 1.0, true, Point3::new(5.0 + 0.1 * TOL, 0.0, 0.0)),
    ];
    let ctx = ResolutionContext {
        operands: [&c, &c],
        records: &records,
    };
    let res = resolve(&cands, &ctx, options(10.0 * TOL)).expect("resolves");
    let blocks = pave_blocks(&cands, &res, [&c, &c], &records, TOL);
    let vertices = vertex_records(&cands, &res, &records, &blocks);
    assert_eq!(vertices.len(), res.groups.len());
    let first = &vertices[res.group_of[0]];
    assert_eq!(
        first.splits,
        [(0, 1, Incompatibility::DistinctOperandVertices)]
    );
    assert!(
        first.consequences.iter().any(|c| matches!(
            c,
            VertexConsequence::MicroBlock {
                owner: BlockOwner::Record(0)
            }
        )),
        "{first:?}"
    );
    let collapsed = &vertices[res.group_of[2]];
    assert_eq!(collapsed.members.len(), 2);
    assert!(collapsed.growth > 0.0);
    assert_eq!(
        collapsed.consequences,
        [VertexConsequence::CollapsedEnds {
            record: 1,
            kind: Collapse::MicroEdge
        }]
    );
    let json = resolution_json(&cands, &Ok(res.clone()), &records, &blocks);
    assert_eq!(json["candidates"].as_array().map(Vec::len), Some(4));
    assert_eq!(
        json["resolution"]["vertices"].as_array().map(Vec::len),
        Some(res.groups.len())
    );
}
