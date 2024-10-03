#![allow(unused)]
#![deny(unused_imports)]

use crate::system::{MolecularSystem, Shell, ShellType};
use nalgebra::DMatrix;

mod overlap;

pub fn overlap(system: &MolecularSystem) -> DMatrix<f64> {
    let mut overlap = DMatrix::zeros(system.n_basis(), system.n_basis());

    let shell = |shell: usize| {
        let &Shell {
            shell_type,
            atom_index,
            ref basis_indices,
        } = &system.shells[shell];

        (shell_type, system.atoms[atom_index].position, basis_indices)
    };

    for i in 0..system.shells.len() {
        let (shell_type_a, pos_a, basis_indices_a) = shell(i);

        for j in i + 1..system.shells.len() {
            let (shell_type_b, pos_b, basis_indices_b) = shell(j);

            match (shell_type_a, shell_type_b) {
                (ShellType::S, ShellType::S) => overlap::ss_overlap(),
                (ShellType::P, ShellType::P) => overlap::pp_overlap(),
                (ShellType::D, ShellType::D) => overlap::dd_overlap(),
                (a, b) => {
                    unimplemented!("integral between shell types {a:?} and {b:?}")
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
