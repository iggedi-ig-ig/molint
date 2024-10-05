use std::hint::black_box;

use criterion::{criterion_group, criterion_main, Criterion};
use molint::{
    basis::BasisSet,
    system::{Atom, MolecularSystem},
};
use nalgebra::Point3;

const HYDROGEN_ATOMS: &[Atom] = &[
    Atom {
        ordinal: 1,
        position: Point3::new(0.0, 0.0, 0.0),
    },
    Atom {
        ordinal: 1,
        position: Point3::new(0.0, 0.0, 1.4),
    },
];
const WATER_ATOMS: &[Atom] = &[
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
];

macro_rules! integral_bench {
    ($name:ident, $function:path) => {
        fn $name(c: &mut Criterion) {
            let basis_sto_3g =
                BasisSet::load("data/basis/STO-3G.json").expect("couldn't load sto-3g basis set'");
            let basis_631g =
                BasisSet::load("data/basis/6-31G.json").expect("couldn't load 6-31g basis set'");

            let hydrogen_sto_3g = MolecularSystem::from_atoms(HYDROGEN_ATOMS, &basis_sto_3g);
            let hydrogen_631g = MolecularSystem::from_atoms(HYDROGEN_ATOMS, &basis_631g);
            let water_631g = MolecularSystem::from_atoms(WATER_ATOMS, &basis_631g);

            let mut group = c.benchmark_group(stringify!($name));
            group.bench_function("H2 STO-3G", |b| {
                b.iter_with_large_drop(|| $function(black_box(&hydrogen_sto_3g)))
            });
            group.bench_function("H2 6-31G", |b| {
                b.iter_with_large_drop(|| $function(black_box(&hydrogen_631g)))
            });
            group.bench_function("H20 6-31G", |b| {
                b.iter_with_large_drop(|| $function(black_box(&water_631g)))
            });
            group.finish();
        }
    };
}

integral_bench!(overlap, molint::overlap);
integral_bench!(kinetic, molint::kinetic);
integral_bench!(nuclear, molint::nuclear);
integral_bench!(eri, molint::eri);

criterion_group!(benches, overlap, kinetic, nuclear, eri);
criterion_main!(benches);
