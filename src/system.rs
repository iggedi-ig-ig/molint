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
                shell_type: ShellType::from_angular(angular_magnitude),
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

    /// Get the data of a single shell
    pub(crate) fn shell(&self, shell_index: usize) -> (ShellType, Point3<f64>, (usize, usize)) {
        let &Shell {
            shell_type,
            atom_index,
            basis_start_index,
            basis_size,
        } = &self.shells[shell_index];

        (
            shell_type,
            self.atoms[atom_index].position,
            (basis_start_index, basis_size),
        )
    }

    /// Get the basis functions that are referenced by a shell
    pub(crate) fn shell_basis<'a>(&'a self, shell_index: usize) -> &[&ContractedGaussian] {
        // TODO: we are cloning here. There are multiple ways of going about avoiding this clone.
        //  The one I have in mind is to sort ContractedGaussians such that Shells only have to
        //  store a start- and end index. Then, slice indices can be used instead of collecting a
        //  vec (because in the end, we want a &[ContractedGaussiasn] anyway)

        let Shell {
            basis_start_index,
            basis_size,
            ..
        } = self.shells[shell_index];
        &self.basis[basis_start_index..basis_start_index + basis_size]
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

#[derive(Copy, Clone, Debug)]
pub enum ShellType {
    /// 0
    S,
    /// 1
    P,
    /// 2
    D,
    /// 3
    F,
    /// 4
    G,
    /// n
    N(u32),
}

impl ShellType {
    pub fn from_angular(i: u32) -> Self {
        match i {
            0 => Self::S,
            1 => Self::P,
            2 => Self::D,
            3 => Self::F,
            4 => Self::G,
            n => Self::N(n),
        }
    }

    pub fn angular(&self) -> u32 {
        match *self {
            Self::S => 0,
            Self::P => 1,
            Self::D => 2,
            Self::F => 3,
            Self::G => 4,
            Self::N(n) => n,
        }
    }
}
