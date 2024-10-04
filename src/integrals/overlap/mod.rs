mod dd;
mod pp;
mod ss;

use nalgebra::DMatrix;

use crate::system::ShellBasis;

use super::utils::hermite_expansion;

/// Function to compute the overlap integrals between two electron shells of arbitrary type
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
        _ => gen_overlap(basis_a, basis_b),
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
    for global_a in start_a..start_a + basis_a.len() {
        // We only want to compute overlaps S_ij where i <= j (as S_ij = S_ji).
        // In this case, this means that we only want to compute S_ab
        // if global_a <= global_b, meaning we can start the loop at
        // *at least* global_a (but at start_b if start_b > global_a)
        for global_b in global_a.max(start_b)..start_b + count_b {
            let i = global_a - start_a;
            let j = global_b - start_b;

            let a = basis_a[i];
            let b = basis_b[j];

            let mut sum = 0.0;

            let [l1, m1, n1] = a.angular.map(|n| n as i32);
            let [l2, m2, n2] = b.angular.map(|n| n as i32);

            for (coeff_a, exp_a) in a.iter() {
                for (coeff_b, exp_b) in b.iter() {
                    let p = exp_a + exp_b;
                    let q = exp_a * exp_b / p;

                    sum += coeff_a
                        * coeff_b
                        * hermite_expansion([l1, l2, 0], diff.x, exp_a, exp_b)
                        * hermite_expansion([m1, m2, 0], diff.y, exp_a, exp_b)
                        * hermite_expansion([n1, n2, 0], diff.z, exp_a, exp_b)
                        * (std::f64::consts::PI / (exp_a + exp_b)).powi(3).sqrt()
                }
            }

            result[(i, j)] = sum;
        }
    }
    result
}
