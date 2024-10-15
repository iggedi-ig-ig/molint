mod ssss;

use nalgebra::{Point3, Vector3};
use ndarray::Array4;

use crate::{
    basis::ContractedGaussian,
    storage::hermite::{ExpansionCoefficients, HermiteCache},
    system::{ShellBasis, ShellType},
};

use super::utils::coulomb_auxiliary;

/// Computes the electron-electron repulsion energy integral between four [ShellBasis].
/// Branches for potentially simplifying conditions based on shell types
/// (for example, (SS|SS) integrals are way simpler to compute than (DD|DD))
pub(crate) fn compute_eri(
    basis_a @ ShellBasis {
        shell_type: type_a, ..
    }: ShellBasis,
    basis_b @ ShellBasis {
        shell_type: type_b, ..
    }: ShellBasis,
    basis_c @ ShellBasis {
        shell_type: type_c, ..
    }: ShellBasis,
    basis_d @ ShellBasis {
        shell_type: type_d, ..
    }: ShellBasis,
    hermite_cache: &HermiteCache,
) -> Array4<f64> {
    match (type_a, type_b, type_c, type_d) {
        (ShellType(0), ShellType(0), ShellType(0), ShellType(0)) => {
            ssss::ssss_eri(basis_a, basis_b, basis_c, basis_d)
        }
        _ => gen_eri(basis_a, basis_b, basis_c, basis_d, hermite_cache),
    }
}

/// Generic eri integral between four electron shells.
fn gen_eri(
    ShellBasis {
        center: pos_a,
        basis: basis_a,
        start_index: start_a,
        count: count_a,
        ..
    }: ShellBasis,
    ShellBasis {
        center: pos_b,
        basis: basis_b,
        start_index: start_b,
        count: count_b,
        ..
    }: ShellBasis,
    ShellBasis {
        center: pos_c,
        basis: basis_c,
        start_index: start_c,
        count: count_c,
        ..
    }: ShellBasis,
    ShellBasis {
        center: pos_d,
        basis: basis_d,
        start_index: start_d,
        count: count_d,
        ..
    }: ShellBasis,
    hermite_cache: &HermiteCache,
) -> Array4<f64> {
    let mut result = Array4::zeros((count_a, count_b, count_c, count_d));

    // Symmetry:
    // when calculating eri (ij|kl), we're ensuring that:
    //  1. i <= j
    //  2. k <= l
    //  3. ij <= kl (hyperindices ij and kl with ab = a * (a + 1) / 2 + b)
    for global_a in start_a..start_a + basis_a.len() {
        for global_b in start_b.max(global_a)..start_b + basis_b.len() {
            let ab = global_a * (global_a + 1) / 2 + global_b;

            let i = global_a - start_a;
            let j = global_b - start_b;

            let a = basis_a[i];
            let b = basis_b[j];

            let expansion_ab = hermite_cache.basis_pair(global_a, global_b);

            for global_c in start_c..start_c + basis_c.len() {
                for global_d in start_d.max(global_c)..start_d + basis_d.len() {
                    let cd = global_c * (global_c + 1) / 2 + global_d;

                    if ab > cd {
                        continue;
                    }

                    let k = global_c - start_c;
                    let l = global_d - start_d;
                    let c = basis_c[k];
                    let d = basis_d[l];

                    let expansion_cd = hermite_cache.basis_pair(global_c, global_d);
                    result[(i, j, k, l)] = contracted_gaussian_eri(
                        [a, b, c, d],
                        [pos_a, pos_b, pos_c, pos_d],
                        [expansion_ab, expansion_cd],
                    );
                }
            }
        }
    }

    result
}

fn contracted_gaussian_eri(
    [a, b, c, d]: [&ContractedGaussian; 4],
    [pos_a, pos_b, pos_c, pos_d]: [Point3<f64>; 4],
    [expansion_ab, expansion_cd]: [&ExpansionCoefficients; 2],
) -> f64 {
    let mut sum = 0.0;

    let angular_vector_ab = Vector3::from(a.angular) + Vector3::from(b.angular);
    let angular_vector_cd = Vector3::from(c.angular) + Vector3::from(d.angular);

    for (i, (coeff_a, exp_a)) in a.iter().enumerate() {
        for (j, (coeff_b, exp_b)) in b.iter().enumerate() {
            let p = exp_a + exp_b;

            // inlined utils::product_center to reuse p
            let product_center_ab = (exp_a * pos_a.coords + exp_b * pos_b.coords) / p;

            for (k, (coeff_c, exp_c)) in c.iter().enumerate() {
                for (l, (coeff_d, exp_d)) in d.iter().enumerate() {
                    let q = exp_c + exp_d;

                    // inlined utils::product_center to reuse q
                    let product_center_cd = (exp_c * pos_c.coords + exp_d * pos_d.coords) / q;

                    let diff_product = product_center_cd - product_center_ab;

                    sum += coeff_a
                        * coeff_b
                        * coeff_c
                        * coeff_d
                        * primitive_eri(
                            [expansion_ab, expansion_cd],
                            [i, j, k, l],
                            angular_vector_ab,
                            angular_vector_cd,
                            [p, q],
                            diff_product,
                        );
                }
            }
        }
    }
    sum
}

#[allow(clippy::too_many_arguments)]
fn primitive_eri(
    [expansion_ab, expansion_cd]: [&ExpansionCoefficients; 2],
    [i, j, k, l]: [usize; 4],
    angular_ab: Vector3<i32>,
    angular_cd: Vector3<i32>,
    [p, q]: [f64; 2],
    diff_product: Vector3<f64>,
) -> f64 {
    let alpha = p * q / (p + q);

    let mut sum = 0.0;
    for t1 in 0..=angular_ab.x {
        for u1 in 0..=angular_ab.y {
            for v1 in 0..=angular_ab.z {
                for t2 in 0..=angular_cd.x {
                    for u2 in 0..=angular_cd.y {
                        for v2 in 0..=angular_cd.z {
                            sum += expansion_ab.coefficient(0, i, j, t1 as usize)
                                * expansion_ab.coefficient(1, i, j, u1 as usize)
                                * expansion_ab.coefficient(2, i, j, v1 as usize)
                                * expansion_cd.coefficient(0, k, l, t2 as usize)
                                * expansion_cd.coefficient(1, k, l, u2 as usize)
                                * expansion_cd.coefficient(2, k, l, v2 as usize)
                                * coulomb_auxiliary(
                                    t1 + t2,
                                    u1 + u2,
                                    v1 + v2,
                                    0,
                                    alpha,
                                    diff_product,
                                )
                                * if (t2 + u2 + v2) % 2 == 0 { 1.0 } else { -1.0 }
                        }
                    }
                }
            }
        }
    }

    2.0 * std::f64::consts::PI.powi(5).sqrt() * (p * q * (p + q).sqrt()).recip() * sum
}
