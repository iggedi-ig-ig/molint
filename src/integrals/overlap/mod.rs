mod dd;
mod pp;
mod ss;

use nalgebra::DMatrix;

use crate::{system::ShellBasis, utils};

/// Generic function to compute the overlap between two electron shells
pub(crate) fn compute_overlap(
    basis_a @ ShellBasis {
        shell_type: type_a, ..
    }: ShellBasis,
    basis_b @ ShellBasis {
        shell_type: type_b, ..
    }: ShellBasis,
) -> DMatrix<f64> {
    // TODO(perf): some combinations of shell types can be simplifed. There are already modules for some
    // of them which are left out
    match (type_a, type_b) {
        shell_types => gen_overlap(basis_a, basis_b),
    }
}

/// Generic overlap integral between two electron shells.
fn gen_overlap(
    ShellBasis {
        center: pos_a,
        basis: basis_a,
        start_index: start_a,
        count: count_a,
        ..
    }: ShellBasis,
    ShellBasis {
        center: pos_b,
        basis: basis_b,
        start_index: start_b,
        count: count_b,
        ..
    }: ShellBasis,
) -> DMatrix<f64> {
    let diff = pos_b - pos_a;

    // Use a matrix to organize results. result[(i, j)] = S_ij
    let mut result = DMatrix::zeros(count_a, count_b);

    // basis_a contains all basis functions that are part of shell A
    // basis_b contains all basis functions taht are part of shell B
    // They may be equal to each other.
    // Now, we compute S_AB which is a block in the total overlap matrix S
    for (i, a) in basis_a.iter().enumerate() {
        let global_index_a = start_a + i;

        for (j, b) in basis_b.iter().enumerate() {
            let global_index_b = start_b + j;

            // we only want to compute overlaps between basis functions k and l if k <= l
            // (because S_kl = S_lk) which in this case means we must ensure that

            // here, we are iterating through pairs of _basis functions_, not primitives.
            // This means that each inner loop here has a unique result which we are interested in.
            // k and l are global_index_a and global_index_b respectively.
            // thus:
            if global_index_b < global_index_a {
                continue;
            }

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
