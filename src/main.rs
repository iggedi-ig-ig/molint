#![allow(dependency_on_unit_never_type_fallback)]
#![allow(unused)]
#![deny(unused_imports)]

mod basis;
mod integrals;
mod periodic_table;
pub mod system;
pub mod types;
mod utils;

use basis::BasisSet;
pub use integrals::{eri, kinetic, nuclear, overlap};
use nalgebra::{Point3, Vector3};
use system::{Atom, MolecularSystem};

// TODO: rename file to lib.rs and remove main.rs
fn main() -> anyhow::Result<()> {
    let basis_set: BasisSet = BasisSet::load("data/basis/STO-3G.json")?;

    let system = MolecularSystem::from_atoms(
        &[
            Atom {
                ordinal: 1,
                position: Point3::origin(),
            },
            Atom {
                ordinal: 1,
                position: Point3::origin() + Vector3::x() * 1.2,
            },
        ],
        &basis_set,
    );

    println!("{}", overlap(&system));
    Ok(())
}
