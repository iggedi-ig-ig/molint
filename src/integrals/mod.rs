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
                    output
                        .view_mut((start_a, start_b), (count_a, count_b))
                        .copy_from(&result);
                }
            }

            log::debug!("{}: {output:2.4}", stringify!($name));
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

            output
                .view_mut((start_a, start_b), (count_a, count_b))
                .copy_from(&result);
        }
    }
    log::debug!("nuclear: {output:2.4}");
    output
}

// TODO: custom tensor type to save on storage (only 1/8 should be needed)
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

                    output
                        .slice_mut(ndarray::s![
                            start_a..start_a + count_a,
                            start_b..start_b + count_b,
                            start_c..start_c + count_c,
                            start_d..start_d + count_d,
                        ])
                        .assign(&result);
                }
            }
        }
    }

    output
}
