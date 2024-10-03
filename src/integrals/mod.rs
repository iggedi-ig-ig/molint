#![allow(unused)]
#![deny(unused_imports)]

use crate::system::{MolecularSystem, ShellType};
use nalgebra::DMatrix;

mod overlap;

pub fn overlap(system: &MolecularSystem) -> DMatrix<f64> {
    let mut overlap = DMatrix::zeros(system.n_basis(), system.n_basis());

    for a in 0..system.shells.len() {
        let (shell_type_a, pos_a, basis_indices_a) = system.shell(a);
        let basis_a = system.shell_basis(a);

        for b in a + 1..system.shells.len() {
            let (shell_type_b, pos_b, basis_indices_b) = system.shell(b);
            let basis_b = system.shell_basis(b);
            let diff = pos_b - pos_a;

            let result = match (shell_type_a, shell_type_b) {
                (ShellType::S, ShellType::S) => overlap::ss_overlap(diff, &basis_a, &basis_b),
                (ShellType::P, ShellType::P) => overlap::pp_overlap(diff, &basis_a, &basis_b),
                (ShellType::D, ShellType::D) => overlap::dd_overlap(diff, &basis_a, &basis_b),
                _ => panic!(),
            };

            // copy result of this particular overlap into the total overlap matrix
            for (i, &a) in basis_indices_a.iter().enumerate() {
                for (j, &b) in basis_indices_b.iter().enumerate() {
                    overlap[(a, b)] = result[(i, j)];
                }
            }
        }
    }

    overlap
}

pub fn kinetic(system: &MolecularSystem) -> DMatrix<f64> {
    todo!()
}

pub fn nuclear(system: &MolecularSystem) -> DMatrix<f64> {
    todo!()
}

/// TODO: write electron tensor type. Or maybe ndarray?
pub fn eri(system: &MolecularSystem) -> () {
    todo!()
}
