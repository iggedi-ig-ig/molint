use ndarray::Array4;

pub fn shell_norm(result: &Array4<f64>) -> f64 {
    let mut norm = 0.0;

    for abab in result.iter() {
        norm = abab.max(norm);
    }

    norm.sqrt()
}
