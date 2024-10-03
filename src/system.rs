use nalgebra::Point3;

use crate::basis::{ConfigBasisSet, ContractedGaussian};

#[derive(Copy, Clone, Debug)]
pub struct Atom {
    /// also the charge of the nucleus
    pub ordinal: usize,
    pub position: Point3<f64>,
}

#[derive(Debug)]
pub struct MolecularSystem {
    pub atoms: Vec<Atom>,
    pub basis: Vec<ContractedGaussian>,
    pub(crate) shells: Vec<Shell>,
    pub charge: i32,
    pub spin_multiplicity: i32,
}

impl MolecularSystem {
    pub fn from_atoms(atoms: &[Atom], basis_set: &ConfigBasisSet) {
        todo!()
    }

    pub fn n_basis(&self) -> usize {
        self.basis.len()
    }

    /// Get the data of a single shell
    pub(crate) fn shell(&self, shell_index: usize) -> (ShellType, Point3<f64>, &[usize]) {
        let &Shell {
            shell_type,
            atom_index,
            ref basis_indices,
        } = &self.shells[shell_index];

        (shell_type, self.atoms[atom_index].position, basis_indices)
    }

    /// Get the basis functions that are referenced by a shell
    pub(crate) fn shell_basis<'a>(&'a self, shell_index: usize) -> Vec<ContractedGaussian> {
        // TODO: we are cloning here. There are multiple ways of going about avoiding this clone.
        //  The one I have in mind is to sort ContractedGaussians such that Shells only have to
        //  store a start- and end index. Then, slice indices can be used instead of collecting a
        //  vec (because in the end, we want a &[ContractedGaussiasn] anyway)
        self.shells[shell_index]
            .basis_indices
            .iter()
            .map(|&i| self.basis[i].clone())
            .collect()
    }
}

#[derive(Debug)]
pub struct Shell {
    pub(crate) shell_type: ShellType,
    pub(crate) atom_index: usize,
    pub(crate) basis_indices: Vec<usize>,
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
