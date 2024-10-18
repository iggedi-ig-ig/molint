//! This module contains the definitions of all logic associated with integral evaluation

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
mod screening;
mod utils;

/// Computes and returns the overlap integral matrix for the given [MolecularSystem] as a [SymmetricMatrix].
pub fn overlap(system: &MolecularSystem) -> SymmetricMatrix {
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
            let result = overlap::compute_overlap(basis_a, basis_b);
            output.copy_from(&result, (start_a, start_b), (count_a, count_b));
        }
    }
    let log_level = log::Level::Trace;
    if log::log_enabled!(log_level) {
        log::log!(
            log_level,
            "{}: {:2.4}",
            stringify!(overlap),
            DMatrix::from(&output)
        );
    }
    output
}

/// Returns the kinetic energy integral matrix for the given [MolecularSystem] as a [SymmetricMatrix]
pub fn kinetic(system: &MolecularSystem) -> SymmetricMatrix {
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
            let result = kinetic::compute_kinetic(basis_a, basis_b);
            output.copy_from(&result, (start_a, start_b), (count_a, count_b));
        }
    }
    let log_level = log::Level::Trace;
    if log::log_enabled!(log_level) {
        log::log!(
            log_level,
            "{}: {:2.4}",
            stringify!(kinetic),
            DMatrix::from(&output)
        );
    }
    output
}

/// Returns the electron-nuclear attraction energy integral matrix for the given [MolecularSystem] as a
/// [SymmetricMatrix]
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
    let log_level = log::Level::Trace;
    if log::log_enabled!(log_level) {
        log::log!(log_level, "nuclear: {:2.4}", DMatrix::from(&output));
    }
    output
}

/// Returns the electron-electron repulsion energy integral tensor for the given [MolecularSystem]
/// as an [EriTensor]
pub fn eri(system: &MolecularSystem) -> EriTensor {
    let n_shells = system.n_shells();

    let start = Instant::now();
    let hermite_cache = HermiteCache::new(system);
    log::debug!(
        "computing hermite expansion coefficient cache took {:3.3?}",
        start.elapsed()
    );

    let start = Instant::now();
    let mut output = EriTensor::zeros(system.n_basis());
    let mut shell_norms = SymmetricMatrix::zeros(n_shells);
    for a in 0..n_shells {
        for b in a..n_shells {
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

            let result = eri::compute_eri(basis_a, basis_b, basis_a, basis_b, &hermite_cache);

            shell_norms[(a, b)] = screening::shell_norm(&result);

            output.copy_from(
                &result,
                (start_a, start_b, start_a, start_b),
                (count_a, count_b, count_a, count_b),
            );
        }
    }
    let diagonal_duration = start.elapsed();
    let prediction = diagonal_duration * n_shells.pow(2) as u32 / 8;
    println!("done precomputing diagonal for screening. took {diagonal_duration:3.3?}. Full is thus estimated to take {prediction:3.3?}");

    const SUFFICIENTLY_SMALL_THRESHOLD: f64 = 1e-6;

    let mut screened = 0;
    let mut total = 0;
    for a in 0..n_shells {
        for b in a..n_shells {
            let norm_ab = shell_norms[(a, b)];

            for c in 0..n_shells {
                for d in c..n_shells {
                    if (a, b) == (c, d) {
                        // skip diagonal entries as they are precomputed for screening
                        continue;
                    }

                    let norm_cd = shell_norms[(c, d)];
                    let norm_abcd = norm_ab * norm_cd;

                    total += 1;
                    if norm_abcd < SUFFICIENTLY_SMALL_THRESHOLD {
                        screened += 1;
                        continue;
                    }

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

    println!(
        "screened {screened} out of {total} shell quartets ({:3.1}%)",
        screened as f32 / total as f32 * 100.0
    );
    println!("computing full ERI tensor took {:3.3?}", start.elapsed());
    output
}
