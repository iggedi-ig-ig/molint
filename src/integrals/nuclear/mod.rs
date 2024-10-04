use nalgebra::{DMatrix, Vector3};

use crate::{basis::ContractedGaussian, system::ShellType};

pub fn compute_nuclear(
    (shell_a, shell_b): (ShellType, ShellType),
    diff: Vector3<f64>,
    basis_a: &[&ContractedGaussian],
    basis_b: &[&ContractedGaussian],
) -> DMatrix<f64> {
    todo!()
}
