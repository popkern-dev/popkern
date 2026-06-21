//! Square pyramid constructor.
//!
//! Builds a square pyramid with an axis-aligned base and an apex
//! directly above the base centre. Used as a test vehicle for
//! non-trihedral vertex blending — the apex has 4 incident edges,
//! one per slant.

use crate::error::CoreError;
use crate::geometry::Point3;
use crate::tolerance::Tolerance;
use crate::topology::*;

impl BRep {
    /// Build an axis-aligned square pyramid with base side `side`
    /// and apex height `height`. The base is centred at the
    /// origin in the `z = 0` plane and the apex sits at
    /// `(0, 0, height)`. The result is a closed manifold solid
    /// with:
    ///
    /// - **5 vertices**: 4 base corners + 1 apex
    /// - **8 edges**: 4 base + 4 slant
    /// - **5 faces**: 1 base quad + 4 slant triangles
    ///
    /// The apex vertex has **4 incident edges**, providing the
    /// 4-valence test vehicle required for non-trihedral vertex
    /// blending (`§5.5`).
    ///
    /// All surfaces are analytic planes and all edge curves are
    /// straight lines. Pcurves are wired via
    /// [`BRep::attach_planar_pcurves`] after topology construction.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if `side` or `height`
    ///   is non-positive or non-finite.
    /// - [`CoreError::Validation`] if the downstream validate
    ///   pass fails (shouldn't happen for well-formed inputs).
    pub fn square_pyramid(side: f64, height: f64, tolerance: Tolerance) -> Result<Self, CoreError> {
        if !side.is_finite() || side <= 0.0 {
            return Err(CoreError::PreconditionFailed(format!(
                "square_pyramid: side must be positive and finite, got {side}"
            )));
        }
        if !height.is_finite() || height <= 0.0 {
            return Err(CoreError::PreconditionFailed(format!(
                "square_pyramid: height must be positive and finite, got {height}"
            )));
        }

        let h = side * 0.5;
        let mut brep = BRep::new(tolerance);

        // ── Vertices ──────────────────────────────────────────
        let v0 = brep.vertices.insert(Vertex {
            position: Point3::new(-h, -h, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let v1 = brep.vertices.insert(Vertex {
            position: Point3::new(h, -h, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let v2 = brep.vertices.insert(Vertex {
            position: Point3::new(h, h, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let v3 = brep.vertices.insert(Vertex {
            position: Point3::new(-h, h, 0.0),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });
        let apex = brep.vertices.insert(Vertex {
            position: Point3::new(0.0, 0.0, height),
            outgoing: Default::default(),
            name: None,
            tolerance: None,
        });

        // ── Solid / shell ─────────────────────────────────────
        let solid_h = brep.solids.insert(Solid {
            shell: Default::default(),
        });
        let shell_h = brep.shells.insert(Shell {
            face: Default::default(),
            solid: solid_h,
        });
        brep.solids[solid_h].shell = shell_h;

        // ── Faces & loops ─────────────────────────────────────
        // 1 base + 4 slant = 5 faces.
        let faces: [FaceHandle; 5] = std::array::from_fn(|_| {
            brep.faces.insert(Face {
                outer_loop: Default::default(),
                inner_loops: Vec::new(),
                shell: shell_h,
                surface: None,
                name: None,
                tolerance: None,
                role: Default::default(),
            })
        });
        let loops: [LoopHandle; 5] = std::array::from_fn(|i| {
            brep.loops.insert(Loop {
                half_edge: Default::default(),
                face: faces[i],
            })
        });
        for (i, lh) in loops.iter().enumerate() {
            brep.faces[faces[i]].outer_loop = *lh;
        }
        brep.shells[shell_h].face = faces[0];
        let [_f_base, _f_front, _f_right, _f_back, _f_left] = faces;
        let [l_base, l_front, l_right, l_back, l_left] = loops;

        // ── Edges ─────────────────────────────────────────────
        // 4 base + 4 slant = 8 edges.
        let edges: [EdgeHandle; 8] = std::array::from_fn(|_| {
            brep.insert_edge_tracked(Edge {
                half_edge: Default::default(),
                curve: None,
                name: None,
                tolerance: None,
                canonical: Default::default(),
                role: Default::default(),
            })
        });
        let [e_01, e_12, e_23, e_30, e_0a, e_1a, e_2a, e_3a] = edges;

        // ── Half-edges ────────────────────────────────────────
        // Base loop walks v0 → v3 → v2 → v1 → v0 (outward normal
        // pointing in −z, so we walk CCW as seen from below).
        let mk_he = |brep: &mut BRep, origin, edge, loop_| {
            brep.half_edges.insert(HalfEdge {
                origin,
                twin: Default::default(),
                next: Default::default(),
                prev: Default::default(),
                edge,
                loop_,
                pcurve: None,
            })
        };

        // Base loop half-edges (origin, edge, loop).
        let hb_v0_v3 = mk_he(&mut brep, v0, e_30, l_base);
        let hb_v3_v2 = mk_he(&mut brep, v3, e_23, l_base);
        let hb_v2_v1 = mk_he(&mut brep, v2, e_12, l_base);
        let hb_v1_v0 = mk_he(&mut brep, v1, e_01, l_base);

        // Front loop (v0 → v1 → apex → v0), normal has −y.
        let hf_v0_v1 = mk_he(&mut brep, v0, e_01, l_front);
        let hf_v1_a = mk_he(&mut brep, v1, e_1a, l_front);
        let hf_a_v0 = mk_he(&mut brep, apex, e_0a, l_front);

        // Right loop (v1 → v2 → apex → v1), normal has +x.
        let hr_v1_v2 = mk_he(&mut brep, v1, e_12, l_right);
        let hr_v2_a = mk_he(&mut brep, v2, e_2a, l_right);
        let hr_a_v1 = mk_he(&mut brep, apex, e_1a, l_right);

        // Back loop (v2 → v3 → apex → v2), normal has +y.
        let hk_v2_v3 = mk_he(&mut brep, v2, e_23, l_back);
        let hk_v3_a = mk_he(&mut brep, v3, e_3a, l_back);
        let hk_a_v2 = mk_he(&mut brep, apex, e_2a, l_back);

        // Left loop (v3 → v0 → apex → v3), normal has −x.
        let hl_v3_v0 = mk_he(&mut brep, v3, e_30, l_left);
        let hl_v0_a = mk_he(&mut brep, v0, e_0a, l_left);
        let hl_a_v3 = mk_he(&mut brep, apex, e_3a, l_left);

        // ── Twin wiring ────────────────────────────────────────
        // Base edges: each base half-edge pairs with a slant face's
        // reverse half-edge.
        let pairs = [
            (hb_v1_v0, hf_v0_v1), // e_01
            (hb_v2_v1, hr_v1_v2), // e_12
            (hb_v3_v2, hk_v2_v3), // e_23
            (hb_v0_v3, hl_v3_v0), // e_30
            // Slant edges: each slant half-edge pairs between
            // two adjacent slant faces.
            (hl_v0_a, hf_a_v0), // e_0a (left ↔ front)
            (hf_v1_a, hr_a_v1), // e_1a (front ↔ right)
            (hr_v2_a, hk_a_v2), // e_2a (right ↔ back)
            (hk_v3_a, hl_a_v3), // e_3a (back ↔ left)
        ];
        for &(a, b) in &pairs {
            brep.half_edges[a].twin = b;
            brep.half_edges[b].twin = a;
        }

        // ── Loop next/prev wiring ──────────────────────────────
        let wire_loop = |brep: &mut BRep, hes: &[HalfEdgeHandle]| {
            let n = hes.len();
            for i in 0..n {
                brep.half_edges[hes[i]].next = hes[(i + 1) % n];
                brep.half_edges[hes[i]].prev = hes[(i + n - 1) % n];
            }
        };
        wire_loop(&mut brep, &[hb_v0_v3, hb_v3_v2, hb_v2_v1, hb_v1_v0]);
        wire_loop(&mut brep, &[hf_v0_v1, hf_v1_a, hf_a_v0]);
        wire_loop(&mut brep, &[hr_v1_v2, hr_v2_a, hr_a_v1]);
        wire_loop(&mut brep, &[hk_v2_v3, hk_v3_a, hk_a_v2]);
        wire_loop(&mut brep, &[hl_v3_v0, hl_v0_a, hl_a_v3]);

        // ── Loop → half_edge / edge → half_edge / vertex → outgoing ──
        brep.loops[l_base].half_edge = hb_v0_v3;
        brep.loops[l_front].half_edge = hf_v0_v1;
        brep.loops[l_right].half_edge = hr_v1_v2;
        brep.loops[l_back].half_edge = hk_v2_v3;
        brep.loops[l_left].half_edge = hl_v3_v0;

        brep.edges[e_01].half_edge = hf_v0_v1;
        brep.edges[e_12].half_edge = hr_v1_v2;
        brep.edges[e_23].half_edge = hk_v2_v3;
        brep.edges[e_30].half_edge = hl_v3_v0;
        brep.edges[e_0a].half_edge = hl_v0_a;
        brep.edges[e_1a].half_edge = hf_v1_a;
        brep.edges[e_2a].half_edge = hr_v2_a;
        brep.edges[e_3a].half_edge = hk_v3_a;

        brep.vertices[v0].outgoing = hf_v0_v1;
        brep.vertices[v1].outgoing = hr_v1_v2;
        brep.vertices[v2].outgoing = hk_v2_v3;
        brep.vertices[v3].outgoing = hl_v3_v0;
        brep.vertices[apex].outgoing = hf_a_v0;

        // ── Geometry: straight-line curves + plane surfaces ────
        brep.attach_geometry()?;
        brep.attach_planar_pcurves()?;

        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }

    /// Build an axis-aligned `n`-gonal pyramid with a regular-polygon
    /// base inscribed in a circle of the given `base_radius` (z = 0)
    /// and an apex at `(0, 0, height)`.
    ///
    /// Generalises [`BRep::square_pyramid`] to N ≥ 3 base sides —
    /// primarily a test vehicle for §5.3 N-valence vertex blending
    /// where N > 4 (e.g., pentagonal / hexagonal hub-and-spoke
    /// geometries). Built via public Euler operators (`mvfs` / `mev`
    /// / `mef`) then decorated with analytic plane surfaces + line
    /// curves via the shared constructor plumbing.
    ///
    /// # Errors
    ///
    /// - [`CoreError::PreconditionFailed`] if `n < 3`, or if
    ///   `base_radius` / `height` is non-positive or non-finite.
    /// - [`CoreError::Validation`] if the downstream validate
    ///   pass fails.
    pub fn n_gonal_pyramid(
        n: usize,
        base_radius: f64,
        height: f64,
        tolerance: Tolerance,
    ) -> Result<Self, CoreError> {
        if n < 3 {
            return Err(CoreError::PreconditionFailed(format!(
                "n_gonal_pyramid: n must be ≥ 3, got {n}"
            )));
        }
        if !base_radius.is_finite() || base_radius <= 0.0 {
            return Err(CoreError::PreconditionFailed(format!(
                "n_gonal_pyramid: base_radius must be positive and finite, got {base_radius}"
            )));
        }
        if !height.is_finite() || height <= 0.0 {
            return Err(CoreError::PreconditionFailed(format!(
                "n_gonal_pyramid: height must be positive and finite, got {height}"
            )));
        }

        // Base vertices on a circle of `base_radius` at z=0.
        let base_pts: Vec<Point3> = (0..n)
            .map(|i| {
                let theta = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
                Point3::new(base_radius * theta.cos(), base_radius * theta.sin(), 0.0)
            })
            .collect();
        let apex_pt = Point3::new(0.0, 0.0, height);

        let mut brep = BRep::new(tolerance);

        // Euler ops: build the base polygon first.
        let mvfs = brep.mvfs(base_pts[0])?;
        let first_mev = brep.mev_at_vertex(mvfs.vertex, mvfs.loop_, base_pts[1])?;
        let mut last_twin = first_mev.half_edge_twin;
        let mut base_vertices: Vec<VertexHandle> = Vec::with_capacity(n);
        base_vertices.push(mvfs.vertex);
        base_vertices.push(first_mev.vertex);
        for pt in &base_pts[2..] {
            let m = brep.mev(last_twin, *pt)?;
            last_twin = m.half_edge_twin;
            base_vertices.push(m.vertex);
        }
        // Close the base polygon: now we have 2 faces (base + pent "top").
        let mef_close = brep.mef(last_twin, first_mev.half_edge)?;

        // The apex is added inside the "top" face (created by mef_close).
        // Find a half-edge in the top face's loop whose origin is base[0].
        let top_loop = brep.faces[mef_close.face].outer_loop;
        let he_b0_in_top = {
            let start = brep.loops[top_loop].half_edge;
            let mut cursor = start;
            loop {
                if brep.half_edges[cursor].origin == base_vertices[0] {
                    break cursor;
                }
                cursor = brep.half_edges[cursor].next;
                if cursor == start {
                    return Err(CoreError::PreconditionFailed(
                        "n_gonal_pyramid: base[0] not in top loop".into(),
                    ));
                }
            }
        };
        // mev inserts the apex spur before he_b0_in_top.
        let apex_mev = brep.mev(he_b0_in_top, apex_pt)?;
        let apex_vh = apex_mev.vertex;

        // Triangulate the top face by connecting apex to each base[i]
        // for i ∈ 1..n. After each mef the "remaining large face"
        // shrinks; we always locate the apex half-edge in the loop
        // with size ≥ 4 (the large face), then find the target base
        // half-edge in that same loop.
        for target_base in base_vertices.iter().copied().skip(1) {
            // Locate the apex half-edge in the large top-face loop.
            let he_at_apex = {
                let mut found = None;
                for (heh, he) in &brep.half_edges {
                    if he.origin != apex_vh {
                        continue;
                    }
                    let loop_h = he.loop_;
                    let start = brep.loops[loop_h].half_edge;
                    let mut size = 0usize;
                    let mut cursor = start;
                    loop {
                        size += 1;
                        cursor = brep.half_edges[cursor].next;
                        if cursor == start || size > 1 + 2 * n {
                            break;
                        }
                    }
                    if size >= 4 {
                        found = Some(heh);
                        break;
                    }
                }
                found.ok_or_else(|| {
                    CoreError::PreconditionFailed(
                        "n_gonal_pyramid: no apex half-edge in large loop".into(),
                    )
                })?
            };
            let he_at_target = {
                let loop_h = brep.half_edges[he_at_apex].loop_;
                let start = brep.loops[loop_h].half_edge;
                let mut cursor = start;
                let mut found = None;
                loop {
                    if brep.half_edges[cursor].origin == target_base {
                        found = Some(cursor);
                        break;
                    }
                    cursor = brep.half_edges[cursor].next;
                    if cursor == start {
                        break;
                    }
                }
                found.ok_or_else(|| {
                    CoreError::PreconditionFailed(
                        "n_gonal_pyramid: target base vertex not in loop".into(),
                    )
                })?
            };
            brep.mef(he_at_apex, he_at_target)?;
        }

        brep.attach_geometry()?;
        brep.attach_planar_pcurves()?;
        brep.validate().map_err(CoreError::Validation)?;

        Ok(brep)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn square_pyramid_topology_counts() {
        let p = BRep::square_pyramid(2.0, 1.5, Tolerance::default()).unwrap();
        assert_eq!(p.vertices.len(), 5);
        assert_eq!(p.edges.len(), 8);
        assert_eq!(p.faces.len(), 5);
    }

    #[test]
    fn square_pyramid_apex_has_four_incident_edges() {
        let p = BRep::square_pyramid(2.0, 1.5, Tolerance::default()).unwrap();
        let apex = p
            .vertices
            .iter()
            .find(|(_, v)| v.position.x == 0.0 && v.position.y == 0.0 && v.position.z > 0.0)
            .map(|(h, _)| h)
            .unwrap();
        // Walk around apex via twin.next and count unique edges.
        let mut edges = std::collections::HashSet::new();
        let start = p.vertices[apex].outgoing;
        let mut cursor = start;
        loop {
            let he = &p.half_edges[cursor];
            edges.insert(he.edge);
            cursor = p.half_edges[he.twin].next;
            if cursor == start {
                break;
            }
        }
        assert_eq!(edges.len(), 4, "apex should have 4 incident edges");
    }

    #[test]
    fn square_pyramid_passes_validation() {
        let p = BRep::square_pyramid(2.0, 1.5, Tolerance::default()).unwrap();
        p.validate().unwrap();
    }

    #[test]
    fn square_pyramid_rejects_invalid_dimensions() {
        assert!(BRep::square_pyramid(-1.0, 1.0, Tolerance::default()).is_err());
        assert!(BRep::square_pyramid(1.0, 0.0, Tolerance::default()).is_err());
        assert!(BRep::square_pyramid(f64::NAN, 1.0, Tolerance::default()).is_err());
    }
}
