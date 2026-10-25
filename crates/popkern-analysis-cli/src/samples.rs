//! Surface samples as CSV (`--samples`): every clipped piece's centroid
//! and corners, where the assessment searches for peaks, with the
//! recovered stress and the displacement there. For comparing against
//! another solver point by point: display vertices can be far coarser
//! than a cell.

use popkern_analysis::AnalysisError;
use popkern_analysis::structural::StructuralResult;

/// CSV header of [`samples_csv`].
pub const HEADER: &str = "solid,face,x,y,z,area,von_mises,max_principal,ux,uy,uz";

/// The surface samples of `result` as CSV (mm, MPa); displacement
/// columns are empty where the field has no value at the sample.
pub fn samples_csv(result: &StructuralResult) -> Result<(String, usize), AnalysisError> {
    let mut out = String::from(HEADER);
    out.push('\n');
    let mut count = 0;
    for (solid, stresses) in result.stresses.iter().enumerate() {
        for s in &stresses.surface {
            let [x, y, z] = s.position;
            let u = match result.field.displacement_at(s.position)? {
                Some([a, b, c]) => format!("{a:e},{b:e},{c:e}"),
                None => ",,".to_string(),
            };
            out.push_str(&format!(
                "{solid},{},{x:e},{y:e},{z:e},{:e},{:e},{:e},{u}\n",
                s.face.0, s.area, s.stress.von_mises, s.stress.max_principal
            ));
            count += 1;
        }
    }
    Ok((out, count))
}
