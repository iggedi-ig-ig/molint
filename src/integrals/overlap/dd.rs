use nalgebra::{DMatrix, Vector3};

use crate::basis::ContractedGaussian;

fn hermite_dd([i, j]: [i32; 2], dist: f64, a: f64, b: f64) -> f64 {
    let p = a + b;
    let q = a * b / p;

    match [i, j] {
        [0, 0] => f64::exp(-q * dist.powi(2)),
        [1, 0] => (2.0 * p).recip() - q * dist / a,
        [0, 1] => (2.0 * p).recip() + q * dist / b,
        [1, 1] => (2.0 * p).recip().powi(2) - q * dist.powi(2) / (a * b),
        _ => panic!("Invalid input for PP hermite expansion: expected [1, 0] or [0, 1]"),
    }
}

pub(super) fn dd_overlap(
    diff: Vector3<f64>,
    basis_a: &[ContractedGaussian],
    basis_b: &[ContractedGaussian],
) -> DMatrix<f64> {
    // Use a matrix to organize results. result[(i, j)] = S_ij
    let mut result = DMatrix::zeros(basis_a.len(), basis_b.len());

    for (i, a) in basis_a.iter().enumerate() {
        // TODO: Is symmetry exploitable here aswell? it feels like that should be the case, but
        // I am way too unsure to just do it lol
        for (j, b) in basis_b.iter().enumerate() {
            // here, we are iterating through pairs of _basis functions_, not primitives.
            // This means that each inner loop here has a unique result which we are interested in

            let mut sum = 0.0;

            let [l1, m1, n1] = a.angular.map(|n| n as i32);
            let [l2, m2, n2] = a.angular.map(|n| n as i32);

            for (coeff_a, exp_a) in a.iter() {
                for (coeff_b, exp_b) in b.iter() {
                    let p = exp_a + exp_b;
                    let q = exp_a * exp_b / p;

                    sum += coeff_a
                        * coeff_b
                        * hermite_dd([l1, l2], diff.x, exp_a, exp_b)
                        * hermite_dd([m1, m2], diff.y, exp_a, exp_b)
                        * hermite_dd([n1, n2], diff.z, exp_a, exp_b)
                        * (std::f64::consts::PI / (exp_a + exp_b)).powi(3).sqrt()
                }
            }

            result[(i, j)] = sum;
        }
    }
    result
}
