use crate::{
    system::{MolecularSystem, ShellBasis},
    EriTensor, SymmetricMatrix,
};
use nalgebra::DMatrix;

mod eri;
mod kinetic;
mod nuclear;
mod overlap;
mod utils;

// TODO(style): probably remove this macro
macro_rules! one_electron_integral {
    ($name:ident, $function:path) => {
        pub fn $name(system: &MolecularSystem) -> SymmetricMatrix {
            let mut output = SymmetricMatrix::zeros(system.n_basis());

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
                        for (j, b) in (start_b..start_b + count_b)
                            .enumerate()
                            .skip_while(|&(_, b)| a > b)
                        {
                            output[(a, b)] = result[(i, j)];
                        }
                    }
                }
            }

            let log_level = log::Level::Debug;
            if log::log_enabled!(log_level) {
                log::log!(
                    log_level,
                    "{}: {:2.4}",
                    stringify!($name),
                    DMatrix::from(&output)
                );
            }
            output
        }
    };
}

one_electron_integral!(overlap, overlap::compute_overlap);
one_electron_integral!(kinetic, kinetic::compute_kinetic);

pub fn nuclear(system: &MolecularSystem) -> SymmetricMatrix {
    let mut output = SymmetricMatrix::zeros(system.n_basis());
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
                for (j, b) in (start_b..start_b + count_b)
                    .enumerate()
                    .skip_while(|&(_, b)| a > b)
                {
                    *output.index_unchecked_mut((a, b)) = result[(i, j)];
                }
            }
        }
    }
    let log_level = log::Level::Debug;
    if log::log_enabled!(log::Level::Debug) {
        log::log!(log_level, "nuclear: {:2.4}", DMatrix::from(&output));
    }
    output
}

pub fn eri(system: &MolecularSystem) -> EriTensor {
    let n_shells = system.shells.len();

    let mut output = EriTensor::zeros(system.n_basis());

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
                        for (j, b) in (start_b..start_b + count_b)
                            .enumerate()
                            .skip_while(|&(_, b)| a > b)
                        {
                            let ab = a * (a + 1) / 2 + b;
                            for (k, c) in (start_c..start_c + count_c).enumerate() {
                                for (l, d) in (start_d..start_d + count_d)
                                    .enumerate()
                                    .skip_while(|&(_, d)| c > d || ab > c * (c + 1) / 2 + d)
                                {
                                    *output.index_unchecked_mut((a, b, c, d)) = result[(i, j, k, l)]
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
