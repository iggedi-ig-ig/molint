#![allow(dependency_on_unit_never_type_fallback)]
#![allow(unused)]
#![deny(unused_imports)]

mod basis;
mod integrals;
mod periodic_table;
pub mod system;
pub mod types;

use basis::BasisSet;
pub use integrals::{eri, kinetic, nuclear, overlap};
use nalgebra::Point3;
use system::{Atom, MolecularSystem};

// TODO: rename file to lib.rs and remove main.rs
fn main() -> anyhow::Result<()> {
    let basis_set: BasisSet = BasisSet::load("data/basis/STO-3G.json")?;

    let system = MolecularSystem::from_atoms(
        &[
            Atom {
                ordinal: 1,
                position: Point3::new(0.4175, 0.0, 0.83),
            },
            Atom {
                ordinal: 8,
                position: Point3::new(0.0, 0.0, -0.31),
            },
            Atom {
                ordinal: 1,
                position: Point3::new(-0.4175, 0.0, 0.83),
            },
        ],
        &basis_set,
    );

    println!("{:0.4}", overlap(&system));
    println!("{:0.4}", kinetic(&system));
    println!("{:0.4}", nuclear(&system));
    Ok(())
}
