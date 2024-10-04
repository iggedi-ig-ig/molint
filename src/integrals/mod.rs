use crate::system::{MolecularSystem, ShellBasis};
use nalgebra::DMatrix;

mod eri;
mod kinetic;
mod nuclear;
mod overlap;
mod utils;

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
one_electron_integral!(nuclear, nuclear::compute_nuclear);

/// TODO: write electron tensor type. Or maybe ndarray?
pub fn eri(system: &MolecularSystem) -> () {
    todo!()
}
