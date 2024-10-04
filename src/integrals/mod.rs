use crate::system::MolecularSystem;
use nalgebra::DMatrix;

mod eri;
mod kinetic;
mod nuclear;
mod overlap;

macro_rules! one_electron_integral {
    ($name:ident, $function:path) => {
        pub fn $name(system: &MolecularSystem) -> DMatrix<f64> {
            let mut output = DMatrix::zeros(system.n_basis(), system.n_basis());

            for a in 0..system.shells.len() {
                let (shell_type_a, pos_a, (basis_index_start_a, basis_size_a)) = system.shell(a);
                let basis_a = system.shell_basis(a);

                for b in a..system.shells.len() {
                    let (shell_type_b, pos_b, (basis_index_start_b, basis_size_b)) =
                        system.shell(b);

                    let basis_b = system.shell_basis(b);
                    let diff = pos_b - pos_a;

                    let result = $function((shell_type_a, shell_type_b), diff, basis_a, basis_b);

                    // copy result of this particular integral into the total result matrix
                    for (i, a) in
                        (basis_index_start_a..basis_index_start_a + basis_size_a).enumerate()
                    {
                        for (j, b) in
                            (basis_index_start_b..basis_index_start_b + basis_size_b).enumerate()
                        {
                            output[(a, b)] = result[(i, j)];
                        }
                    }
                }
            }

            output
        }
    };
}

one_electron_integral!(overlap, overlap::compute_overlap);
one_electron_integral!(kinetic, kinetic::compute_kinetic);
one_electron_integral!(nuclear, nuclear::compute_nuclear);

/// TODO: write electron tensor type. Or maybe ndarray?
pub fn eri(system: &MolecularSystem) -> () {
    todo!()
}
