use crate::system::{MolecularSystem, ShellBasis};
use nalgebra::DMatrix;
use ndarray::Array4;

// TODO(perf): for all integrals, think about row / column majorness of matrices
//  (and the effects indexing order has on performance based on that)
mod eri;
mod kinetic;
mod nuclear;
mod overlap;
mod utils;

// TODO(style): probably remove this macro
macro_rules! one_electron_integral {
    ($name:ident, $function:path) => {
        pub fn $name(system: &MolecularSystem) -> DMatrix<f64> {
            let mut output = DMatrix::zeros(system.n_basis(), system.n_basis());

            for a in 0..system.shells.len() {
                let basis_a @ ShellBasis {
                    start_index: start_a,
                    count: count_a,
                    ..
                } = system.shell_basis(a);

                for b in a..system.shells.len() {
                    let basis_b @ ShellBasis {
                        start_index: start_b,
                        count: count_b,
                        ..
                    } = system.shell_basis(b);

                    let result = $function(basis_a, basis_b);

                    // copy result of this particular integral into the total result matrix
                    for (i, a) in (start_a..start_a + count_a).enumerate() {
                        for (j, b) in (start_b..start_b + count_b).enumerate() {
                            output[(a, b)] = result[(i, j)];
                        }
                    }
                }
            }

            output
        }
    };
}

one_electron_integral!(overlap, overlap::compute_overlap);
one_electron_integral!(kinetic, kinetic::compute_kinetic);

pub fn nuclear(system: &MolecularSystem) -> DMatrix<f64> {
    let mut output = DMatrix::zeros(system.n_basis(), system.n_basis());
    for a in 0..system.shells.len() {
        let basis_a @ ShellBasis {
            start_index: start_a,
            count: count_a,
            ..
        } = system.shell_basis(a);
        for b in a..system.shells.len() {
            let basis_b @ ShellBasis {
                start_index: start_b,
                count: count_b,
                ..
            } = system.shell_basis(b);
            let result = nuclear::compute_nuclear(basis_a, basis_b, system);
            for (i, a) in (start_a..start_a + count_a).enumerate() {
                for (j, b) in (start_b..start_b + count_b).enumerate() {
                    output[(a, b)] = result[(i, j)];
                }
            }
        }
    }
    output
}

/// TODO: write electron tensor type. Or maybe ndarray?
pub fn eri(system: &MolecularSystem) -> Array4<f64> {
    let n_basis = system.basis.len();
    let n_shells = system.shells.len();

    let mut output = Array4::zeros([n_basis; 4]);

    for a in 0..n_shells {
        for b in a..n_shells {
            for c in 0..n_shells {
                for d in c..n_shells {
                    let basis_a @ ShellBasis {
                        start_index: start_a,
                        count: count_a,
                        ..
                    } = system.shell_basis(a);
                    let basis_b @ ShellBasis {
                        start_index: start_b,
                        count: count_b,
                        ..
                    } = system.shell_basis(b);
                    let basis_c @ ShellBasis {
                        start_index: start_c,
                        count: count_c,
                        ..
                    } = system.shell_basis(c);
                    let basis_d @ ShellBasis {
                        start_index: start_d,
                        count: count_d,
                        ..
                    } = system.shell_basis(d);

                    let result = eri::compute_eri(basis_a, basis_b, basis_c, basis_d);
                    for (i, a) in (start_a..start_a + count_a).enumerate() {
                        for (j, b) in (start_b..start_b + count_b).enumerate() {
                            for (k, c) in (start_c..start_c + count_c).enumerate() {
                                for (l, d) in (start_d..start_d + count_d).enumerate() {
                                    output[(a, b, c, d)] = result[(i, j, k, l)]
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    output
}

fn int_template() {
    let (start_a, start_b) = (0, 0);
    let (count_a, count_b) = (1, 1);
    let (basis_a, basis_b) = ([false], [false]);
    let mut result = DMatrix::zeros(count_a, count_b);

    // TEMPLATE BELOW HERE

    for global_a in start_a..start_a + basis_a.len() {
        for global_b in global_a.max(start_b)..start_b + count_b {
            let i = global_a - start_a;
            let j = global_b - start_b;

            let a = basis_a[i];
            let b = basis_b[j];

            let mut sum = 0.0;

            result[(i, j)] = sum;
        }
    }

    // TEMPLATE ABOVE HERE
}
