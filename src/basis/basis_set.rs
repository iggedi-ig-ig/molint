use std::{collections::HashMap, fs::File, path::Path};

use crate::{periodic_table::ElementType, system::Atom};

use super::{bse_basis_set::BseBasisSet, ContractedGaussian};

/// This type represents a basis set that can be used as a basis in the integral evaluation.
pub struct BasisSet(pub(crate) HashMap<ElementType, Vec<ContractedGaussian>>);

impl BasisSet {
    /// Given a path, this function tries to load a basis set from a json file with the format that
    /// basissetexchange.org provides by default.
    pub fn load(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        let basis_set: BseBasisSet = serde_json::from_reader(File::open(path)?)?;
        basis_set.try_into()
    }

    /// Returns the basis functions for the given atom.
    pub(crate) fn atomic_basis(&self, atom: &Atom) -> &[ContractedGaussian] {
        let element_type = ElementType::from_ordinal(atom.ordinal)
            .unwrap_or_else(|| panic!("failed to convert ordinal {} to ElementType", atom.ordinal));
        &self.0[&element_type]
    }
}
