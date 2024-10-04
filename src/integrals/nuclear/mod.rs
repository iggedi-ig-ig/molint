use nalgebra::DMatrix;

use crate::system::{MolecularSystem, ShellBasis};

use super::utils::{coulomb_auxiliary, hermite_expansion, product_center};

pub(crate) fn compute_nuclear(
    basis_a @ ShellBasis {
        shell_type: type_a, ..
    }: ShellBasis,
    basis_b @ ShellBasis {
        shell_type: type_b, ..
    }: ShellBasis,
    system: &MolecularSystem,
) -> DMatrix<f64> {
    // TODO(perf): specific implementations for simple shell types
    gen_nuclear(basis_a, basis_b, system)
}

fn gen_nuclear(
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
    system: &MolecularSystem,
) -> DMatrix<f64> {
    let diff = pos_b - pos_a;

    let mut result = DMatrix::zeros(count_a, count_b);

    // TODO(style): deeply nested loop

    for global_a in start_a..start_a + basis_a.len() {
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
                    let product_center = product_center(exp_a, pos_a, exp_b, pos_b);

                    for atom in &system.atoms {
                        let nuclear_charge = atom.ordinal as f64;
                        let diff_nucl = atom.position - product_center;

                        let mut atom_sum = 0.0;
                        for t in 0..=l1 + l2 {
                            for u in 0..=m1 + m2 {
                                for v in 0..=n1 + n2 {
                                    let e1 = hermite_expansion([l1, l2, t], diff.x, exp_a, exp_b);
                                    let e2 = hermite_expansion([m1, m2, u], diff.y, exp_a, exp_b);
                                    let e3 = hermite_expansion([n1, n2, v], diff.z, exp_a, exp_b);

                                    atom_sum +=
                                        e1 * e2 * e3 * coulomb_auxiliary(t, u, v, 0, p, diff_nucl);
                                }
                            }
                        }

                        sum += coeff_a
                            * coeff_b
                            * atom_sum
                            * (-nuclear_charge * std::f64::consts::TAU / p);
                    }
                }

                result[(i, j)] = sum;
            }
        }
    }

    result
}
