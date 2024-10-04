use nalgebra::DMatrix;

use crate::system::{MolecularSystem, ShellBasis};

pub(crate) fn compute_nuclear(
    basis_a @ ShellBasis {
        shell_type: type_a, ..
    }: ShellBasis,
    basis_b @ ShellBasis {
        shell_type: type_b, ..
    }: ShellBasis,
    system: &MolecularSystem,
) -> DMatrix<f64> {
    // TODO(perf): some combinations of shell types can be simplifed. There are already modules for some
    // of them which are left out
    match (type_a, type_b) {
        _ => gen_nuclear(basis_a, basis_b, system),
    }
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
    let mut result = DMatrix::zeros(count_a, count_b);

    for atom in &system.atoms {
        let charge = atom.ordinal;
    }

    result
}
