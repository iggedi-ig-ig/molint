use nalgebra::Vector3;

use crate::basis::ContractedGaussian;

#[derive(Copy, Clone, Debug)]
pub struct Atom {
    /// also the charge of the nucleus
    pub ordinal: usize,
    pub position: Vector3<f64>,
}

#[derive(Debug)]
pub struct MolecularSystem {
    pub atoms: Vec<Atom>,
    pub basis: Vec<ContractedGaussian>,
    pub(crate) shells: Vec<Shell>,
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
