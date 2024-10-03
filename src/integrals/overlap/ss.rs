use nalgebra::{DMatrix, Vector3};

use crate::basis::ContractedGaussian;

pub(crate) fn ss_overlap(
    diff: Vector3<f64>,
    basis_a: &[ContractedGaussian],
    basis_b: &[ContractedGaussian],
) -> DMatrix<f64> {
    let length_squared = diff.norm_squared();

    // Use a matrix to organize results. result[(i, j)] = S_ij
    let mut result = DMatrix::zeros(basis_a.len(), basis_b.len());

    for (i, a) in basis_a.iter().enumerate() {
        // TODO: Is symmetry exploitable here aswell? it feels like that should be the case, but
        // I am way too unsure to just do it lol
        for (j, b) in basis_b.iter().enumerate() {
            // here, we are iterating through pairs of _basis functions_, not primitives.
            // This means that each inner loop here has a unique result which we are interested in

            let mut sum = 0.0;

            for (coeff_a, exp_a) in a.iter() {
                for (coeff_b, exp_b) in b.iter() {
                    let p = exp_a + exp_b;
                    let q = exp_a * exp_b / p;

                    // exp(-q * r^2) is a simplification of the hermite expansion function where
                    // angular terms are zero.
                    sum += coeff_a
                        * coeff_b
                        * f64::exp(-q * length_squared)
                        * (std::f64::consts::PI / (exp_a + exp_b)).powi(3).sqrt()
                }
            }

            result[(i, j)] = sum;
        }
    }
    result
}
