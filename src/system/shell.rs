use nalgebra::Point3;

use crate::basis::ContractedGaussian;

#[derive(Copy, Clone, Debug)]
pub struct Shell {
    pub(crate) shell_type: ShellType,
    pub(crate) atom_index: usize,
    /// where in the list of basis functions does this shell "start" (i.e., where is the first of
    /// this shells basis functions in that list)
    pub(super) basis_start_index: usize,
    /// how many basis functions does this shell contain
    pub(super) basis_size: usize,
}

/// Represents the basis of a specific shell in the basis of some molecular system.
/// It's sort of a view into the molecular system, indexed by [Shell].
//
// TODO(style): shell_type and the angular term of contracted gaussians technically imply each
//  other, so it's not technically necessary to store both.
#[derive(Copy, Clone, Debug)]
pub struct ShellBasis<'b> {
    pub(crate) shell_type: ShellType,
    pub(crate) center: Point3<f64>,
    pub(crate) basis: &'b [&'b ContractedGaussian],
    pub(crate) start_index: usize,
    pub(crate) count: usize,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct ShellType(pub(crate) i32);
