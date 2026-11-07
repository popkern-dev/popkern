//! The regions a glued complex of faces encloses (PK-1108, S9b).
//!
//! Faces assembled into one body edge by edge (a STEP shell whose edges
//! some face set uses more than twice, an internal face between two
//! cells, cavities touching each other) start as sheets in the void:
//! every side of every face in the outside's shells. Once their edges
//! are glued into rings, the sides a wedge connects form the pieces of
//! a cell structure, each closed piece one boundary of one cell of
//! space. Read with its sides signed (the region behind a back side, in
//! front of a front side), a piece enclosing positive volume is a
//! bounded cell's outer boundary: a solid when it holds a back side (the
//! faces' material lies behind them), else a bounded void (the inside of
//! a cavity, of a closed sheet). A piece enclosing negative volume is a
//! cell's inner boundary, seen from the cell around it: it bounds the
//! smallest bounded region containing it by winding number other than
//! the cell across its own faces (a cavity's material side is an inner
//! shell of the solid, not of the void inside the cavity; a solid's
//! outside in a closed sheet's void is that void's), or the outside when
//! none does. A face between two cells has its back in one piece and its
//! front in the other, so it bounds both. Free shells of the outside (a
//! sheet, wires, an acorn) lying in a bounded region join it.

use crate::geometry::Point3;
use crate::topology::{BRep, FaceHandle, Region, RegionHandle, RegionKind, ShellHandle};

/// A bounded region read by its outer boundary.
struct HostCell {
    region: RegionHandle,
    sides: Vec<(FaceHandle, f64)>,
    volume: f64,
}

/// One closed piece, its sides signed.
struct Piece {
    shell: ShellHandle,
    sides: Vec<(FaceHandle, f64)>,
    volume: f64,
    has_back: bool,
}

impl BRep {
    /// Give every closed piece of the outside's shells the region it
    /// bounds (see the module docs): a new solid or bounded void for a
    /// positive piece, the enclosing cell for a negative one. The
    /// outside's shells are first separated into their pieces (sides
    /// connected through wedges); free shells then join the smallest
    /// region holding them. Returns the regions made.
    pub fn fill_enclosed_cells(&mut self) -> Vec<RegionHandle> {
        let outside = self.outside;
        let shells: Vec<ShellHandle> = self
            .shells
            .iter()
            .filter(|(_, s)| s.region == outside)
            .map(|(h, _)| h)
            .collect();
        for shell in shells {
            self.separate_shell_pieces(shell);
        }
        let open = self.open_shells();
        let pieces: Vec<Piece> = self
            .shells
            .iter()
            .filter(|(h, s)| s.region == outside && !open.contains(h))
            .filter_map(|(shell, _)| {
                let sides: Vec<(FaceHandle, f64)> = self
                    .faces
                    .iter()
                    .filter_map(|(h, f)| match (f.back == shell, f.front == shell) {
                        (true, false) => Some((h, 1.0)),
                        (false, true) => Some((h, -1.0)),
                        _ => None,
                    })
                    .collect();
                let volume = self.enclosed_signed_volume(&sides)?;
                let has_back = sides.iter().any(|s| s.1 > 0.0);
                Some(Piece {
                    shell,
                    sides,
                    volume,
                    has_back,
                })
            })
            .collect();
        // Positive pieces: each the outer boundary of a new cell.
        let mut made = Vec::new();
        for piece in pieces.iter().filter(|p| p.volume > 0.0) {
            let kind = if piece.has_back {
                RegionKind::Solid
            } else {
                RegionKind::Void
            };
            let region = self.regions.insert(Region {
                kind,
                shell: piece.shell,
            });
            self.shells[piece.shell].region = region;
            made.push(region);
        }
        // Every bounded region is a candidate host, read by its outer
        // boundary: the smallest containing a point is the innermost.
        let hosts = self.host_cells();
        // Negative pieces: an inner boundary of the smallest cell holding
        // them, never the cell on the far side of their own faces (a
        // cavity's material side and the void inside it share faces).
        for piece in pieces.iter().filter(|p| p.volume < 0.0) {
            let Some(q) = self.piece_point(piece) else {
                continue;
            };
            let own: std::collections::BTreeSet<FaceHandle> =
                piece.sides.iter().map(|s| s.0).collect();
            let host = hosts
                .iter()
                .filter(|cell| !cell.sides.iter().any(|s| own.contains(&s.0)))
                .filter(|cell| self.winding_about(&cell.sides, q) > 0.5)
                .min_by(|a, b| a.volume.total_cmp(&b.volume))
                .map(|cell| cell.region);
            if let Some(region) = host {
                self.shells[piece.shell].region = region;
            }
        }
        // Free shells of the outside (sheets, wires, acorns) lying in a
        // cell are that cell's: a lamina or wire floating in a solid, an
        // acorn in a cavity.
        let free: Vec<ShellHandle> = self
            .shells
            .keys()
            .filter(|&s| self.is_free_shell(s))
            .collect();
        for shell in free {
            let Some(q) = self.free_shell_point(shell) else {
                continue;
            };
            let host = hosts
                .iter()
                .filter(|cell| self.winding_about(&cell.sides, q) > 0.5)
                .min_by(|a, b| a.volume.total_cmp(&b.volume))
                .map(|cell| cell.region);
            if let Some(region) = host {
                self.shells[shell].region = region;
            }
        }
        made
    }

    /// Every bounded region with its outer shell's sides signed and the
    /// volume they enclose.
    fn host_cells(&self) -> Vec<HostCell> {
        self.regions
            .iter()
            .filter(|&(r, _)| r != self.outside)
            .filter_map(|(region, data)| {
                let shell = data.shell;
                let sides: Vec<(FaceHandle, f64)> = self
                    .faces
                    .iter()
                    .filter_map(|(h, f)| match (f.back == shell, f.front == shell) {
                        (true, false) => Some((h, 1.0)),
                        (false, true) => Some((h, -1.0)),
                        _ => None,
                    })
                    .collect();
                let volume = self.enclosed_signed_volume(&sides)?;
                Some(HostCell {
                    region,
                    sides,
                    volume,
                })
            })
            .collect()
    }

    /// A point of free shell `shell` off any other surface: an interior
    /// sample of one of its faces, else a wire edge's midpoint, else its
    /// acorn.
    fn free_shell_point(&self, shell: ShellHandle) -> Option<Point3> {
        if let Some(face) = self
            .faces
            .iter()
            .find(|(_, f)| f.back == shell)
            .map(|(h, _)| h)
        {
            return self.face_interior_sample(face);
        }
        let data = self.shells.get(shell)?;
        if let Some(&edge) = data.wire_edges.first() {
            let (a, b) = self.edge_vertices(edge)?;
            if let Some(c) = self.edges.get(edge).and_then(|e| e.curve.as_ref()) {
                use crate::geometry::curve::CurveEval;
                return Some(c.geometry.point_at(0.5 * (c.t_start + c.t_end)));
            }
            let (p, q) = (
                self.vertices.get(a)?.position,
                self.vertices.get(b)?.position,
            );
            return Some(Point3::new(
                0.5 * (p.x + q.x),
                0.5 * (p.y + q.y),
                0.5 * (p.z + q.z),
            ));
        }
        data.acorn
            .and_then(|v| self.vertices.get(v))
            .map(|v| v.position)
    }

    /// A point inside the volume `piece` bounds: the mean of its faces'
    /// boundary samples (a closed surface's samples average to a point
    /// it encloses when the surface is not far from convex).
    fn piece_point(&self, piece: &Piece) -> Option<Point3> {
        let pts: Vec<Point3> = piece
            .sides
            .iter()
            .flat_map(|&(f, _)| self.face_boundary_sample_points(f, 2))
            .collect();
        if pts.is_empty() {
            return None;
        }
        let n = pts.len() as f64;
        Some(pts.iter().fold(Point3::new(0.0, 0.0, 0.0), |a, p| {
            Point3::new(a.x + p.x / n, a.y + p.y / n, a.z + p.z / n)
        }))
    }
}
