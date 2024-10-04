use std::collections::HashMap;

use nalgebra::Point3;

use crate::basis::{BasisSet, ContractedGaussian};

#[derive(Copy, Clone, Debug)]
pub struct Atom {
    /// also the charge of the nucleus
    pub ordinal: usize,
    pub position: Point3<f64>,
}

#[derive(Debug)]
/// Represents the quantum system of a molecule.
pub struct MolecularSystem<'b> {
    pub atoms: Vec<Atom>,
    pub basis: Vec<&'b ContractedGaussian>,
    pub(crate) shells: Vec<Shell>,
}

/// Represents the basis of a specific shell in the basis of some molecular system.
/// It's sort of a view into the molecular system, indexed by [Shell].
//
// TODO(style): should this even store shell_type? because it can technically be deduced from the
//  angular terms of the contracted gaussians. Or maybe contracted gaussians shouldn't even store
//  their angular terms, because they are implied by the shell type? This might also reduce memory
//  usage
#[derive(Copy, Clone, Debug)]
pub(crate) struct ShellBasis<'b> {
    pub(crate) shell_type: ShellType,
    pub(crate) center: Point3<f64>,
    pub(crate) basis: &'b [&'b ContractedGaussian],
    pub(crate) start_index: usize,
    pub(crate) count: usize,
}

impl<'b> MolecularSystem<'b> {
    /// Create a molecular system given the atom types and positons and a basis set.
    /// The basis set must outlive this object.
    pub fn from_atoms(atoms: &[Atom], basis_set: &'b BasisSet) -> Self {
        let mut shell_map = HashMap::new();

        // collect basis functions based on shells
        for (i, atom) in atoms.into_iter().enumerate() {
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
    pub(crate) fn shell_basis<'a>(&'a self, shell_index: usize) -> ShellBasis {
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

#[derive(Copy, Clone, Debug)]
pub struct Shell {
    pub(crate) shell_type: ShellType,
    pub(crate) atom_index: usize,
    /// where in the list of basis functions does this shell "start" (i.e., where is the first of
    /// this shells basis functions in that list)
    basis_start_index: usize,
    /// how many basis functions does this shell contain
    basis_size: usize,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ShellType(u32);
