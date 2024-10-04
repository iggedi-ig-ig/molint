use nalgebra::{DMatrix, Vector3};

use crate::{integrals::utils::hermite_expansion, system::ShellBasis};

/// Function to compute the kinetic energy integrals between two electron shells of arbitrary type
pub(crate) fn compute_kinetic(
    basis_a @ ShellBasis {
        shell_type: type_a, ..
    }: ShellBasis,
    basis_b @ ShellBasis {
        shell_type: type_b, ..
    }: ShellBasis,
) -> DMatrix<f64> {
    match (type_a, type_b) {
        _ => gen_kinetic(basis_a, basis_b),
    }
}

fn gen_kinetic(
    ShellBasis {
        shell_type: type_a,
        center: pos_a,
        basis: basis_a,
        start_index: start_a,
        count: count_a,
    }: ShellBasis,
    ShellBasis {
        shell_type: type_b,
        center: pos_b,
        basis: basis_b,
        start_index: start_b,
        count: count_b,
    }: ShellBasis,
) -> DMatrix<f64> {
    // a lot of stuff that is helpful to understand this function is documented in the gen_overlap
    // function in 'src/integrals/overlap/mod.rs'.

    let diff = pos_b - pos_a;

    let mut result = DMatrix::zeros(count_a, count_b);

    // this neseted loop is weird - for better explanation, see the comments in gen_overlap in
    // 'integrals/overlap/mod.rs'
    for global_a in start_a..start_a + basis_a.len() {
        for global_b in global_a.max(start_b)..start_b + count_b {
            let i = global_a - start_a;
            let j = global_b - start_b;

            let a = basis_a[i];
            let b = basis_b[j];

            let [l1, m1, n1] = a.angular.map(|n| n as i32);
            let [l2, m2, n2] = b.angular.map(|n| n as i32);

            let mut sum = 0.0;

            for (coeff_a, exp_a) in a.iter() {
                for (coeff_b, exp_b) in b.iter() {
                    // don't know a good name to call this
                    let angular_step = |i: i32, j: i32, k: i32| {
                        primitive_overlap(
                            exp_a,
                            exp_b,
                            (l1, m1, n1),
                            (l2 + i, m2 + j, n2 + k),
                            diff,
                        )
                    };

                    // TODO(perf): maybe some of the overlap / angular step calculations can be
                    // factored out?
                    let term0 = exp_b
                        * 2.0
                        * (l2 + m2 + n2 + 3) as f64
                        * primitive_overlap(exp_a, exp_b, (l1, m1, n1), (l2, m2, n2), diff);
                    let term1 = -2.0
                        * exp_b.powi(2)
                        * (angular_step(2, 0, 0) + angular_step(0, 2, 0) + angular_step(0, 0, 2));
                    let term2 = -0.5
                        * ((l2 * (l2 - 1)) as f64 * angular_step(-2, 0, 0)
                            + (m2 * (m2 - 1)) as f64 * angular_step(0, -2, 0)
                            + (n2 * (n2 - 1)) as f64 * angular_step(0, 0, -2));

                    sum += coeff_a * coeff_b * (term0 + term1 + term2);
                }
            }

            result[(i, j)] = sum;
        }
    }

    result
}

fn primitive_overlap(
    exp_a: f64,
    exp_b: f64,
    (l1, m1, n1): (i32, i32, i32),
    (l2, m2, n2): (i32, i32, i32),
    diff: Vector3<f64>,
) -> f64 {
    (std::f64::consts::PI / (exp_a + exp_b)).powi(3).sqrt()
        * hermite_expansion([l1, l2, 0], diff.x, exp_a, exp_b)
        * hermite_expansion([m1, m2, 0], diff.y, exp_a, exp_b)
        * hermite_expansion([n1, n2, 0], diff.z, exp_a, exp_b)
}
