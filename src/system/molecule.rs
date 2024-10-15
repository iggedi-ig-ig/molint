use indexmap::IndexMap;
use itertools::Itertools;
use std::{fs::File, path::Path};

use crate::{
    basis::{BasisSet, ContractedGaussian},
    system::ShellType,
};

use super::{config_atom::ConfigAtom, shell::Shell, Atom, ShellBasis};

#[derive(Debug)]
/// Represents the quantum system of a molecule, represented in some [BasisSet]
pub struct MolecularSystem<'b> {
    /// The constituent atoms of this system
    pub atoms: Vec<Atom>,
    /// references to the [ContractedGaussian]s of the [BasisSet] that is used to represent this
    /// system.
    pub basis: Vec<&'b ContractedGaussian>,
    /// The [Shell]s that this system has.
    pub(crate) shells: Vec<Shell>,
}

impl<'a> MolecularSystem<'a> {
    pub fn load(path: impl AsRef<Path>, basis_set: &'a BasisSet) -> anyhow::Result<Self> {
        let config_atoms: Vec<ConfigAtom> = serde_json::from_reader(File::open(path)?)?;
        let atoms: Vec<Atom> = config_atoms.into_iter().map(Atom::try_from).try_collect()?;

        Ok(Self::from_atoms(&atoms, basis_set))
    }
}

impl<'b> MolecularSystem<'b> {
    /// Create a molecular system given the atom types and positons and a basis set.
    /// The basis set must outlive this object.
    pub fn from_atoms(atoms: &[Atom], basis_set: &'b BasisSet) -> Self {
        let mut shell_map = IndexMap::new();

        // collect basis functions based on shells
        for (i, atom) in atoms.iter().enumerate() {
            let atomic_basis = basis_set.atomic_basis(atom);

            for basis_function in atomic_basis {
                shell_map
                    .entry((i, basis_function.angular))
                    .and_modify(|v: &mut Vec<_>| v.push(basis_function))
                    .or_insert_with(|| vec![basis_function]);
            }
        }

        let mut shells = Vec::with_capacity(shell_map.len());
        let mut basis = Vec::new();

        for ((atom_index, [i, j, k]), basis_functions) in shell_map {
            let angular_magnitude = i + j + k;

            let shell = Shell {
                shell_type: ShellType(angular_magnitude),
                atom_index,
                basis_start_index: basis.len(),
                basis_size: basis_functions.len(),
            };
            shells.push(shell);
            basis.extend_from_slice(&basis_functions);

            // invariant: basis[start_index..start_index+size] == basis_functions
            assert_eq!(
                &basis_functions,
                &basis[shell.basis_start_index..shell.basis_start_index + shell.basis_size],
            );
        }

        log::info!("loaded molecular system with {} atoms and {} basis functions, which were decomposed into {} shells", atoms.len(), basis.len(), shells.len());

        Self {
            atoms: atoms.to_vec(),
            basis,
            shells,
        }
    }

    pub fn n_basis(&self) -> usize {
        self.basis.len()
    }

    /// Get the concrete shell basis of a shell in this system  
    pub(crate) fn shell_basis(&self, shell_index: usize) -> ShellBasis {
        let Shell {
            shell_type,
            atom_index,
            basis_start_index,
            basis_size,
        } = self.shells[shell_index];

        ShellBasis {
            shell_type,
            center: self.atoms[atom_index].position,
            basis: &self.basis[basis_start_index..basis_start_index + basis_size],
            start_index: basis_start_index,
            count: basis_size,
        }
    }
}
