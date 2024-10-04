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
use nalgebra::{DMatrix, Point3};
use system::{Atom, MolecularSystem};

// TODO: rename file to lib.rs and remove main.rs
fn main() -> anyhow::Result<()> {
    pretty_env_logger::init();

    let basis_set: BasisSet = BasisSet::load("data/basis/6-31G.json")?;
    let system = MolecularSystem::from_atoms(
        &[
            Atom {
                ordinal: 1,
                position: Point3::origin(),
            },
            Atom {
                ordinal: 1,
                position: Point3::new(0.0, 0.0, 1.4),
            },
        ],
        &basis_set,
    );

    let reflect = |m: DMatrix<f64>| {
        DMatrix::from_fn(m.nrows(), m.ncols(), |i, j| {
            let (i, j) = if i <= j { (i, j) } else { (j, i) };
            m[(i, j)]
        })
    };

    let overlap = overlap(&system);
    log::debug!("{:0.4}", reflect(overlap).trace());
    let kinetic = kinetic(&system);
    log::debug!("{:0.4}", reflect(kinetic).trace());
    let nuclear = nuclear(&system);
    log::debug!("{:0.4}", reflect(nuclear).trace());
    let electron = eri(&system);
    log::debug!("{:0.8}", electron);
    Ok(())
}
