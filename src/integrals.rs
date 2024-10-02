use nalgebra::DMatrix;

use crate::system::MolecularSystem;

/// This is only used for convenient passing around of values
/// (thus it is defined here, not in crate::basis)
pub(crate) struct GaussianPrimitive {
    pub coefficient: f64,
    pub exponent: f64,
    pub angular: (i32, i32, i32),
}

pub fn overlap(system: &MolecularSystem) -> DMatrix<f64> {
    todo!()
}

pub fn kinetic(system: &MolecularSystem) -> DMatrix<f64> {
    todo!()
}

pub fn nuclear(system: &MolecularSystem) -> DMatrix<f64> {
    todo!()
}

/// TODO: write electron tensor type
pub fn eri(system: &MolecularSystem) -> () {
    todo!()
}
