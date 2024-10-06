# molint
A crate to compute the integrals necessary for quantum chemistry calculations, written in pure rust. 
I would like to call it "blazingly fast", but I don't think that's accurate. Maybe some day

# Usage
```rs
use molint::basis::BasisSet;
use molint::system::MolecularSystem;

// uses the json format basissetexchange.org uses for basis sets
let basis_set = BasisSet::open("path/to/basis_set.json").unwrap();

// H2, for example:
// [
//     {
//         "element": "1",
//         "position": [0.0, 0.0, 0.0]
//     },
//     {
//         "element": "1",
//         "position": [0.0, 0.0, 1.4]
//     }
// ]
let system = MolecularSystem::open("path/to/molecule.json", &basis_set).unwrap();

let overlap = molint::overlap(&system);
let kinetic = molint::kinetic(&system);
let nuclear = molint::nuclear(&system);
let eris = molint::eri(&system);
```
