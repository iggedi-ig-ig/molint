//! This module contains types that represent the system for which matrix elements should be
//! calculated.

mod atom;
mod config_atom;
mod molecule;
mod shell;

pub use atom::Atom;
pub use molecule::MolecularSystem;
pub(crate) use shell::{ShellBasis, ShellType};
