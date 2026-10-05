//! Text as exact polygon sets (feature `text`).
//!
//! TrueType glyph outlines are flattened to simple polygons with
//! counter holes classified by containment parity, cleaned to a
//! minimum edge length (~1.5% of the cap height — micro edges read
//! as degenerate to downstream tessellation), and laid out with the
//! font's advance widths, scaled so the cap height is 1.0.

use crate::engrave::HoledPolygon;
use crate::error::PartError;
use popkern_math::PortableMath;

/// A loaded font.
pub struct Font {
    data: Vec<u8>,
}

impl Font {
    /// Load a font from a TTF/OTF file.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when the file cannot be read or
    /// parsed.
    pub fn from_file(path: impl AsRef<std::path::Path>) -> Result<Font, PartError> {
        let path = path.as_ref();
        let data = std::fs::read(path).map_err(|e| {
            PartError::NotApplicable(format!("cannot read font '{}': {e}", path.display()))
        })?;
        ttf_parser::Face::parse(&data, 0).map_err(|e| {
            PartError::NotApplicable(format!("cannot parse font '{}': {e}", path.display()))
        })?;
        Ok(Font { data })
    }

    /// Load the first available common system sans font (override
    /// with an explicit [`Font::from_file`]).
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] listing the tried paths.
    pub fn system_sans() -> Result<Font, PartError> {
        let candidates = [
            "/usr/share/fonts/truetype/noto/NotoSans-Bold.ttf",
            "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
            "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        ];
        for c in candidates {
            if std::path::Path::new(c).exists() {
                return Font::from_file(c);
            }
        }
        Err(PartError::NotApplicable(format!(
            "no system font found (tried {candidates:?}); use Font::from_file"
        )))
    }

    /// Lay out `text` as glyph polygons: outer contours with
    /// parity-classified counter holes, cap height 1.0, pen
    /// advancing along +x.
    ///
    /// # Errors
    ///
    /// [`PartError::NotApplicable`] when the font fails to parse or
    /// produces no outlines for the text.
    pub fn polygons(&self, text: &str) -> Result<Vec<HoledPolygon>, PartError> {
        let face = ttf_parser::Face::parse(&self.data, 0)
            .map_err(|e| PartError::NotApplicable(format!("font parse: {e}")))?;
        let upem = face.units_per_em() as f64;
        let cap = face
            .capital_height()
            .map(|h| h as f64)
            .unwrap_or(upem * 0.7);
        let scale = 1.0 / cap;
        let mut pen_x = 0.0_f64;
        let mut contours_all: Vec<Vec<(f64, f64)>> = Vec::new();
        for ch in text.chars() {
            let Some(gid) = face.glyph_index(ch) else {
                continue;
            };
            let mut collector = OutlineCollector {
                contours: Vec::new(),
                current: Vec::new(),
                start: (0.0, 0.0),
            };
            face.outline_glyph(gid, &mut collector);
            ttf_parser::OutlineBuilder::close(&mut collector);
            for contour in collector.contours {
                let scaled: Vec<(f64, f64)> = contour
                    .into_iter()
                    .map(|(x, y)| ((x + pen_x) * scale, y * scale))
                    .collect();
                if let Some(cleaned) = clean_contour(scaled) {
                    contours_all.push(cleaned);
                }
            }
            pen_x += face
                .glyph_hor_advance(gid)
                .map(|a| a as f64)
                .unwrap_or(upem * 0.5);
        }
        if contours_all.is_empty() {
            return Err(PartError::NotApplicable(format!(
                "text '{text}' produced no glyph outlines"
            )));
        }
        Ok(classify_contours(contours_all))
    }
}

/// Minimum edge length in cap-height units.
const MIN_EDGE: f64 = 0.015;

/// Drop micro edges and enforce ≥ 3 points.
fn clean_contour(mut scaled: Vec<(f64, f64)>) -> Option<Vec<(f64, f64)>> {
    let mut cleaned: Vec<(f64, f64)> = Vec::with_capacity(scaled.len());
    for p in scaled.drain(..) {
        if cleaned
            .last()
            .map(|q| (p.0 - q.0).phypot(p.1 - q.1) >= MIN_EDGE)
            .unwrap_or(true)
        {
            cleaned.push(p);
        }
    }
    while cleaned.len() > 3 {
        let (first, last) = (cleaned[0], *cleaned.last()?);
        if (first.0 - last.0).phypot(first.1 - last.1) < MIN_EDGE {
            cleaned.pop();
        } else {
            break;
        }
    }
    (cleaned.len() >= 3).then_some(cleaned)
}

/// Containment-parity classification: even depth = outer, odd =
/// hole of its innermost containing outer. All output CCW.
fn classify_contours(contours: Vec<Vec<(f64, f64)>>) -> Vec<HoledPolygon> {
    let depth_of = |idx: usize| -> usize {
        let pt = contours[idx][0];
        contours
            .iter()
            .enumerate()
            .filter(|(j, poly)| *j != idx && point_in_poly(pt, poly))
            .count()
    };
    let mut outers: Vec<HoledPolygon> = Vec::new();
    let mut holes: Vec<Vec<(f64, f64)>> = Vec::new();
    for (i, c) in contours.iter().enumerate() {
        let mut c = c.clone();
        if shoelace(&c) < 0.0 {
            c.reverse();
        }
        if depth_of(i) % 2 == 0 {
            outers.push(HoledPolygon {
                outer: c,
                holes: Vec::new(),
            });
        } else {
            holes.push(c);
        }
    }
    for hole in holes {
        let pt = hole[0];
        if let Some(owner) = outers
            .iter_mut()
            .filter(|o| point_in_poly(pt, &o.outer))
            .min_by(|a, b| {
                shoelace(&a.outer)
                    .abs()
                    .partial_cmp(&shoelace(&b.outer).abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
        {
            owner.holes.push(hole);
        }
    }
    outers
}

/// Signed area × 2 of a polygon.
pub(crate) fn shoelace(poly: &[(f64, f64)]) -> f64 {
    poly.iter()
        .zip(poly.iter().cycle().skip(1))
        .take(poly.len())
        .map(|(a, b)| a.0 * b.1 - b.0 * a.1)
        .sum::<f64>()
        / 2.0
}

fn point_in_poly(pt: (f64, f64), poly: &[(f64, f64)]) -> bool {
    let mut inside = false;
    let n = poly.len();
    for i in 0..n {
        let (a, b) = (poly[i], poly[(i + 1) % n]);
        if (a.1 > pt.1) != (b.1 > pt.1) {
            let x = a.0 + (pt.1 - a.1) / (b.1 - a.1) * (b.0 - a.0);
            if x > pt.0 {
                inside = !inside;
            }
        }
    }
    inside
}

/// TrueType outline collector: flattens Béziers to polylines.
struct OutlineCollector {
    contours: Vec<Vec<(f64, f64)>>,
    current: Vec<(f64, f64)>,
    start: (f64, f64),
}

impl ttf_parser::OutlineBuilder for OutlineCollector {
    fn move_to(&mut self, x: f32, y: f32) {
        self.start = (x as f64, y as f64);
        self.current.push(self.start);
    }
    fn line_to(&mut self, x: f32, y: f32) {
        self.current.push((x as f64, y as f64));
    }
    fn quad_to(&mut self, x1: f32, y1: f32, x: f32, y: f32) {
        let p0 = *self.current.last().unwrap_or(&self.start);
        let (c, p) = ((x1 as f64, y1 as f64), (x as f64, y as f64));
        for k in 1..=4 {
            let t = k as f64 / 4.0;
            let mt = 1.0 - t;
            self.current.push((
                mt * mt * p0.0 + 2.0 * mt * t * c.0 + t * t * p.0,
                mt * mt * p0.1 + 2.0 * mt * t * c.1 + t * t * p.1,
            ));
        }
    }
    fn curve_to(&mut self, x1: f32, y1: f32, x2: f32, y2: f32, x: f32, y: f32) {
        let p0 = *self.current.last().unwrap_or(&self.start);
        let (c1, c2, p) = (
            (x1 as f64, y1 as f64),
            (x2 as f64, y2 as f64),
            (x as f64, y as f64),
        );
        for k in 1..=5 {
            let t = k as f64 / 5.0;
            let mt = 1.0 - t;
            self.current.push((
                mt * mt * mt * p0.0
                    + 3.0 * mt * mt * t * c1.0
                    + 3.0 * mt * t * t * c2.0
                    + t * t * t * p.0,
                mt * mt * mt * p0.1
                    + 3.0 * mt * mt * t * c1.1
                    + 3.0 * mt * t * t * c2.1
                    + t * t * t * p.1,
            ));
        }
    }
    fn close(&mut self) {
        let mut contour = std::mem::take(&mut self.current);
        if let (Some(&first), Some(&last)) = (contour.first(), contour.last())
            && (first.0 - last.0).phypot(first.1 - last.1) < 1e-6
        {
            contour.pop();
        }
        contour.dedup_by(|a, b| (a.0 - b.0).phypot(a.1 - b.1) < 1e-6);
        if contour.len() >= 3 {
            self.contours.push(contour);
        }
    }
}
