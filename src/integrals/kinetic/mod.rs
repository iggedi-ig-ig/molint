use nalgebra::DMatrix;

use crate::{integrals::utils::hermite_expansion, system::ShellBasis};

/// Function to compute the kinetic energy integrals between two electron shells of arbitrary type
pub(crate) fn compute_kinetic(
    basis_a @ ShellBasis {
        shell_type: type_a, ..
    }: ShellBasis,
    basis_b @ ShellBasis {
        shell_type: type_b, ..
    }: ShellBasis,
) -> DMatrix<f64> {
    match (type_a, type_b) {
        _ => gen_kinetic(basis_a, basis_b),
    }
}

fn gen_kinetic(
    ShellBasis {
        shell_type: type_a,
        center: pos_a,
        basis: basis_a,
        start_index: start_a,
        count: count_a,
    }: ShellBasis,
    ShellBasis {
        shell_type: type_b,
        center: pos_b,
        basis: basis_b,
        start_index: start_b,
        count: count_b,
    }: ShellBasis,
) -> DMatrix<f64> {
    // a lot of stuff that is helpful to understand this function is documented in the gen_overlap
    // function in 'src/integrals/overlap/mod.rs'.

    let diff = pos_b - pos_a;

    let mut result = DMatrix::zeros(count_a, count_b);

    let angular_magnitude_1 = type_a.magnitude();
    let angular_magnitude_2 = type_b.magnitude();

    // there unfortunately just aren't nice names for coeffients like these :/
    let overlap_coeffient = 2.0 * angular_magnitude_2 as f64;

    let overlaps_1d = ();
    for global_a in start_a..start_a + basis_a.len() {
        let start_j = global_a.max(start_b);
        for global_b in start_j..start_b + count_b {
            let a = basis_a[global_a - start_a];
            let b = basis_b[global_b - start_b];
        }
    }

    for (i, a) in basis_a.iter().enumerate() {
        let global_index_a = start_a + i;
        for (j, b) in basis_b.iter().enumerate() {
            let global_index_b = start_b + j;

            if global_index_b < global_index_a {
                continue;
            }

            let mut sum = 0.0;

            for (coeff_a, exp_a) in a.iter() {
                for (coeff_b, exp_b) in b.iter() {
                    let term_b = exp_b * 2.0 * angular_magnitude_2 as f64;
                }
            }

            result[(i, j)] = sum;
        }
    }

    result
}

// unnormalized overlap between two 1d hermite gaussians
fn overlap_1d(exp_a: f64, exp_b: f64, l_a: i32, l_b: i32, d: f64) -> f64 {
    hermite_expansion([l_a, l_b, 0], d, exp_a, exp_b)
}
