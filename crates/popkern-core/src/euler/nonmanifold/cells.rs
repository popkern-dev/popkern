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
//! smallest positive cell that contains it (a cavity's material side is
//! an inner shell of the solid; a piece of front sides inside a void
//! cell is that void's), or the outside when none does. A face between
//! two cells has its back in one piece and its front in the other, so
//! it bounds both.

use crate::geometry::Point3;
use crate::topology::{BRep, FaceHandle, Region, RegionHandle, RegionKind, ShellHandle};

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
    /// connected through wedges). Returns the regions made.
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
        let mut cells: Vec<(RegionHandle, &Piece)> = Vec::new();
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
            cells.push((region, piece));
        }
        // Negative pieces: an inner boundary of the smallest cell holding
        // them whose kind their sides face.
        for piece in pieces.iter().filter(|p| p.volume < 0.0) {
            let Some(q) = self.piece_point(piece) else {
                continue;
            };
            let want = if piece.has_back {
                RegionKind::Solid
            } else {
                RegionKind::Void
            };
            let host = cells
                .iter()
                .filter(|(r, _)| self.region_kind(*r) == Some(want))
                .filter(|(_, cell)| self.winding_about(&cell.sides, q) > 0.5)
                .min_by(|a, b| a.1.volume.total_cmp(&b.1.volume))
                .map(|(r, _)| *r);
            if let Some(region) = host {
                self.shells[piece.shell].region = region;
            }
        }
        made
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
