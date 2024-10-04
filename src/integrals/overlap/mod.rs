mod dd;
mod pp;
mod ss;

use dd::dd_overlap;
use pp::pp_overlap;
use ss::ss_overlap;

use nalgebra::{DMatrix, Vector3};

use crate::{basis::ContractedGaussian, system::ShellType, utils};

/// Generic function to compute the overlap between two electron shells
pub(crate) fn compute_overlap(
    shell_types: (ShellType, ShellType),
    diff: Vector3<f64>,
    basis_a: &[&ContractedGaussian],
    basis_b: &[&ContractedGaussian],
) -> DMatrix<f64> {
    // TODO: are there other combinations that simplify the calculation?
    match shell_types {
        (ShellType::S, ShellType::S) => ss_overlap(diff, basis_a, basis_b),
        (ShellType::P, ShellType::P) => pp_overlap(diff, basis_a, basis_b),
        (ShellType::D, ShellType::D) => dd_overlap(diff, basis_a, basis_b),
        shell_types => gen_overlap(shell_types, diff, basis_a, basis_b),
    }
}

/// Generic overlap integral between two electron shells
fn gen_overlap(
    shell_types: (ShellType, ShellType),
    diff: Vector3<f64>,
    basis_a: &[&ContractedGaussian],
    basis_b: &[&ContractedGaussian],
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
                        * utils::hermite_expansion([l1, l2, 0], diff.x, exp_a, exp_b)
                        * utils::hermite_expansion([m1, m2, 0], diff.y, exp_a, exp_b)
                        * utils::hermite_expansion([n1, n2, 0], diff.z, exp_a, exp_b)
                        * (std::f64::consts::PI / (exp_a + exp_b)).powi(3).sqrt()
                }
            }

            result[(i, j)] = sum;
        }
    }
    result
}
