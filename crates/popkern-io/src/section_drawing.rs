//! 2-D section drawings (PK-856): the closed loops a plane cuts from a
//! part or an assembly, per body, written as SVG (hatched regions, holes
//! by the even-odd rule) or DXF (closed polylines, a layer per body).
//! The loops arrive already sampled in the section plane's 2-D frame
//! (`popkern_bool::ops::plane_section`); this module only draws them.
//!
//! The SVG is readable back ([`SectionDrawing::from_svg`]): coordinates
//! are written in Rust's shortest round-trip form, so a drawing read
//! from its own SVG is the drawing that wrote it.

use std::fmt::Write as _;

/// A section drawing: the regions each body leaves on the plane.
#[derive(Debug, Clone, PartialEq, Default)]
pub struct SectionDrawing {
    /// One entry per body (a part's body, an assembly's occurrence).
    pub bodies: Vec<SectionBody>,
}

/// One body's section.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionBody {
    /// The body's name (an occurrence path, `body N`).
    pub name: String,
    /// Its regions.
    pub regions: Vec<SectionArea>,
}

/// One region: an outer loop and its holes, in the plane's 2-D frame.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionArea {
    /// The outer boundary, first point not repeated.
    pub outer: Vec<[f64; 2]>,
    /// The holes.
    pub holes: Vec<Vec<[f64; 2]>>,
}

/// Why an SVG could not be read as a section drawing.
#[derive(Debug, thiserror::Error, PartialEq)]
pub enum SectionDrawingError {
    /// A path's data was not `M x y L x y … Z` subpaths.
    #[error("section SVG path {path}: {reason}")]
    Path {
        /// The path's index in the document.
        path: usize,
        /// What was wrong.
        reason: String,
    },
}

impl SectionDrawing {
    /// The drawing's 2-D bounds `(min, max)`, `None` when empty.
    pub fn bounds(&self) -> Option<([f64; 2], [f64; 2])> {
        let mut out: Option<([f64; 2], [f64; 2])> = None;
        let points = self.bodies.iter().flat_map(|b| {
            b.regions
                .iter()
                .flat_map(|r| r.outer.iter().chain(r.holes.iter().flatten()))
        });
        for p in points {
            let (lo, hi) = out.get_or_insert((*p, *p));
            for k in 0..2 {
                lo[k] = lo[k].min(p[k]);
                hi[k] = hi[k].max(p[k]);
            }
        }
        out
    }

    /// The drawing as SVG: one group per body, each region one path
    /// (outer loop then holes, even-odd filled) hatched in the body's
    /// own direction and colour, y up as in the plane's frame.
    pub fn to_svg(&self) -> String {
        let ([x0, y0], [x1, y1]) = self.bounds().unwrap_or(([0.0, 0.0], [1.0, 1.0]));
        let pad = 0.05 * (x1 - x0).max(y1 - y0).max(1e-9);
        let (w, h) = (x1 - x0 + 2.0 * pad, y1 - y0 + 2.0 * pad);
        let stroke = 0.002 * w.max(h);
        let mut svg = String::new();
        let _ = writeln!(
            svg,
            r#"<svg xmlns="http://www.w3.org/2000/svg" viewBox="{} {} {} {}">"#,
            x0 - pad,
            -(y1 + pad),
            w,
            h
        );
        svg.push_str("<defs>\n");
        for (k, _) in self.bodies.iter().enumerate() {
            let (hue, angle) = (k * 137 % 360, 45 + (k % 4) * 30);
            let gap = 0.02 * w.max(h);
            let _ = writeln!(
                svg,
                r#"<pattern id="hatch{k}" patternUnits="userSpaceOnUse" width="{gap}" height="{gap}" patternTransform="rotate({angle})"><rect width="{gap}" height="{gap}" fill="hsl({hue},60%,92%)"/><line x1="0" y1="0" x2="0" y2="{gap}" stroke="hsl({hue},60%,35%)" stroke-width="{stroke}"/></pattern>"#
            );
        }
        svg.push_str("</defs>\n<g transform=\"scale(1,-1)\">\n");
        for (k, body) in self.bodies.iter().enumerate() {
            let _ = writeln!(
                svg,
                r#"<g data-body="{}" fill="url(#hatch{k})" fill-rule="evenodd" stroke="black" stroke-width="{stroke}">"#,
                escape(&body.name)
            );
            for region in &body.regions {
                let mut d = String::new();
                for l in std::iter::once(&region.outer).chain(&region.holes) {
                    for (i, p) in l.iter().enumerate() {
                        let _ = write!(d, "{}{} {} ", if i == 0 { "M" } else { "L" }, p[0], p[1]);
                    }
                    d.push_str("Z ");
                }
                let _ = writeln!(svg, r#"<path d="{}"/>"#, d.trim_end());
            }
            svg.push_str("</g>\n");
        }
        svg.push_str("</g>\n</svg>\n");
        svg
    }

    /// Read a drawing back from [`SectionDrawing::to_svg`]'s output.
    ///
    /// # Errors
    ///
    /// [`SectionDrawingError`] for path data that is not closed `M`/`L`
    /// subpaths of number pairs.
    pub fn from_svg(svg: &str) -> Result<SectionDrawing, SectionDrawingError> {
        let mut drawing = SectionDrawing::default();
        let mut path_index = 0;
        for line in svg.lines() {
            let line = line.trim();
            if let Some(rest) = line.strip_prefix(r#"<g data-body=""#) {
                let name = rest.split('"').next().unwrap_or_default();
                drawing.bodies.push(SectionBody {
                    name: unescape(name),
                    regions: Vec::new(),
                });
            } else if let Some(rest) = line.strip_prefix(r#"<path d=""#) {
                let d = rest.split('"').next().unwrap_or_default();
                let region = parse_region(d).map_err(|reason| SectionDrawingError::Path {
                    path: path_index,
                    reason,
                })?;
                path_index += 1;
                match drawing.bodies.last_mut() {
                    Some(body) => body.regions.push(region),
                    None => {
                        return Err(SectionDrawingError::Path {
                            path: path_index - 1,
                            reason: "a path outside any body group".into(),
                        });
                    }
                }
            }
        }
        Ok(drawing)
    }

    /// The drawing as an ASCII DXF (R12): every loop a closed POLYLINE on
    /// its body's layer (`body name`, made DXF-safe). DXF readers fill
    /// regions themselves; the hatching of the SVG is not written.
    pub fn to_dxf(&self) -> String {
        let mut dxf = String::from("0\nSECTION\n2\nENTITIES\n");
        for body in &self.bodies {
            let layer = dxf_layer(&body.name);
            for region in &body.regions {
                for l in std::iter::once(&region.outer).chain(&region.holes) {
                    let _ = write!(dxf, "0\nPOLYLINE\n8\n{layer}\n66\n1\n70\n1\n");
                    for p in l {
                        let _ = write!(dxf, "0\nVERTEX\n8\n{layer}\n10\n{}\n20\n{}\n", p[0], p[1]);
                    }
                    dxf.push_str("0\nSEQEND\n");
                }
            }
        }
        dxf.push_str("0\nENDSEC\n0\nEOF\n");
        dxf
    }
}

fn parse_region(d: &str) -> Result<SectionArea, String> {
    let mut loops: Vec<Vec<[f64; 2]>> = Vec::new();
    let mut current: Option<Vec<[f64; 2]>> = None;
    let mut tokens = d.split_whitespace();
    while let Some(token) = tokens.next() {
        if token == "Z" {
            loops.push(current.take().ok_or("Z before M")?);
            continue;
        }
        let (command, x) = token.split_at(1);
        let y = tokens.next().ok_or("a coordinate pair cut short")?;
        let p = [
            x.parse::<f64>().map_err(|e| format!("{x}: {e}"))?,
            y.parse::<f64>().map_err(|e| format!("{y}: {e}"))?,
        ];
        match command {
            "M" if current.is_none() => current = Some(vec![p]),
            "L" => current.as_mut().ok_or("L before M")?.push(p),
            other => return Err(format!("unexpected command {other}")),
        }
    }
    if current.is_some() {
        return Err("a subpath without Z".into());
    }
    let mut loops = loops.into_iter();
    let outer = loops.next().ok_or("no subpath")?;
    Ok(SectionArea {
        outer,
        holes: loops.collect(),
    })
}

fn escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('"', "&quot;")
        .replace('<', "&lt;")
}

fn unescape(s: &str) -> String {
    s.replace("&lt;", "<")
        .replace("&quot;", "\"")
        .replace("&amp;", "&")
}

/// A DXF layer name: letters, digits, `-` and `_` only.
fn dxf_layer(name: &str) -> String {
    let layer: String = name
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if layer.is_empty() {
        "SECTION".into()
    } else {
        layer
    }
}

#[cfg(test)]
#[path = "section_drawing_tests.rs"]
mod tests;
