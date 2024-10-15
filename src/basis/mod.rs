//! This module contains types that are associated with either basis functions or with basis sets.
//! It exposes exactly two types directly, namely, [BasisSet] and [ContractedGaussian] which is all

mod basis_set;
pub(super) mod bse_basis_set;
mod contracted_gaussian;

pub use basis_set::BasisSet;
pub use contracted_gaussian::ContractedGaussian;
