pub struct Atom {
    /// also the charge of the nucleus
    pub ordinal: usize,
    pub position: Vector3<f64>,
}

pub struct MolecularSystem {
    atoms: Vec<Atom>,
    basis: Vec<ContractedGaussian>,
    shells: Vec<Shell>,
}

pub struct Shell {
    shell_type: ShellType,
    atom_index: usize,
    basis_indices: Vec<usize>,
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
