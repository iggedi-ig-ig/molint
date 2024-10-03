use nalgebra::Vector3;

use crate::{basis::ContractedGaussian, system::ShellType};

pub(crate) fn compute_electron_repulsion(
    shell_types: (ShellType, ShellType, ShellType, ShellType),
    diff_ab: Vector3<f64>,
    diff_cd: Vector3<f64>,
    basis_a: &[ContractedGaussian],
    basis_b: &[ContractedGaussian],
    basis_c: &[ContractedGaussian],
    basis_d: &[ContractedGaussian],
) -> () {
    todo!()
}
