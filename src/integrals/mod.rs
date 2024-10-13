use std::time::Instant;

use crate::{
    storage::{hermite::HermiteCache, EriTensor, SymmetricMatrix},
    system::{MolecularSystem, ShellBasis},
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
                    output.copy_from(&result, (start_a, start_b), (count_a, count_b));
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
            output.copy_from(&result, (start_a, start_b), (count_a, count_b))
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

    let start = Instant::now();
    let hermite_cache = HermiteCache::new(system);
    log::debug!(
        "computing hermite expansion coefficient cache took {:3.3?}",
        start.elapsed()
    );

    let start = Instant::now();
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

                    let result =
                        eri::compute_eri(basis_a, basis_b, basis_c, basis_d, &hermite_cache);

                    output.copy_from(
                        &result,
                        (start_a, start_b, start_c, start_d),
                        (count_a, count_b, count_c, count_d),
                    );
                }
            }
        }
    }

    log::debug!("computing full ERI tensor took {:3.3?}", start.elapsed());
    output
}
