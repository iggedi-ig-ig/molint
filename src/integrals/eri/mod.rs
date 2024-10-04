use nalgebra::Vector3;
use ndarray::Array4;

use crate::system::ShellBasis;

use super::utils::{coulomb_auxiliary, hermite_expansion};

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
) -> Array4<f64> {
    gen_eri(basis_a, basis_b, basis_c, basis_d)
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
) -> Array4<f64> {
    let diff_ab = pos_b - pos_a;
    let diff_cd = pos_d - pos_c;

    let mut result = Array4::zeros((count_a, count_b, count_c, count_d));

    for global_a in start_a..start_a + basis_a.len() {
        for global_b in start_b..start_b + basis_b.len() {
            for global_c in start_c..start_c + basis_c.len() {
                for global_d in start_d..start_d + basis_d.len() {
                    let i = global_a - start_a;
                    let j = global_b - start_b;
                    let k = global_c - start_c;
                    let l = global_d - start_d;

                    let a = basis_a[i];
                    let b = basis_b[j];
                    let c = basis_c[k];
                    let d = basis_d[l];

                    let mut sum = 0.0;

                    let angular_a @ [l1, m1, n1] = a.angular.map(|n| n as i32);
                    let angular_b @ [l2, m2, n2] = b.angular.map(|n| n as i32);
                    let angular_c @ [l3, m3, n3] = c.angular.map(|n| n as i32);
                    let angular_d @ [l4, m4, n4] = d.angular.map(|n| n as i32);

                    for (coeff_a, exp_a) in a.iter() {
                        for (coeff_b, exp_b) in b.iter() {
                            let p = exp_a + exp_b;

                            // inlined utils::product_center
                            let product_center_ab =
                                (exp_a * pos_a.coords + exp_b * pos_b.coords) / p;

                            let e1s: Vec<_> = (0..=l1 + l2)
                                .map(|k| hermite_expansion([l1, l2, k], diff_ab.x, exp_a, exp_b))
                                .collect();
                            let e2s: Vec<_> = (0..=m1 + m2)
                                .map(|k| hermite_expansion([m1, m2, k], diff_ab.y, exp_a, exp_b))
                                .collect();
                            let e3s: Vec<_> = (0..=n1 + n2)
                                .map(|k| hermite_expansion([n1, n2, k], diff_ab.z, exp_a, exp_b))
                                .collect();

                            for (coeff_c, exp_c) in c.iter() {
                                for (coeff_d, exp_d) in d.iter() {
                                    let q = exp_c + exp_d;

                                    // inlined utils::product_center
                                    let product_center_cd =
                                        (exp_c * pos_c.coords + exp_d * pos_d.coords) / q;

                                    let diff_product = product_center_cd - product_center_ab;

                                    let e4s: Vec<_> = (0..=l3 + l4)
                                        .map(|k| {
                                            hermite_expansion([l3, l4, k], diff_cd.x, exp_c, exp_d)
                                        })
                                        .collect();
                                    let e5s: Vec<_> = (0..=m3 + m4)
                                        .map(|k| {
                                            hermite_expansion([m3, m4, k], diff_cd.y, exp_c, exp_d)
                                        })
                                        .collect();
                                    let e6s: Vec<_> = (0..=n3 + n4)
                                        .map(|k| {
                                            hermite_expansion([n3, n4, k], diff_cd.z, exp_c, exp_d)
                                        })
                                        .collect();

                                    sum += coeff_a
                                        * coeff_b
                                        * coeff_c
                                        * coeff_d
                                        * inner_sum(
                                            &[&e1s, &e2s, &e3s, &e4s, &e5s, &e6s],
                                            angular_a,
                                            angular_b,
                                            angular_c,
                                            angular_d,
                                            p,
                                            q,
                                            diff_product,
                                        );
                                }
                            }
                        }
                    }

                    result[(i, j, k, l)] = sum;
                }
            }
        }
    }

    result
}

fn inner_sum(
    expansion_coeffs: &[&[f64]],
    [l1, m1, n1]: [i32; 3],
    [l2, m2, n2]: [i32; 3],
    [l3, m3, n3]: [i32; 3],
    [l4, m4, n4]: [i32; 3],
    p: f64,
    q: f64,
    diff_product: Vector3<f64>,
) -> f64 {
    let alpha = p * q / (p + q);

    let mut sum = 0.0;
    for t1 in 0..=l1 + l2 {
        for u1 in 0..=m1 + m2 {
            for v1 in 0..=n1 + n2 {
                for t2 in 0..=l3 + l4 {
                    for u2 in 0..=m3 + m4 {
                        for v2 in 0..=n3 + n4 {
                            sum += expansion_coeffs[0][t1 as usize]
                                * expansion_coeffs[1][u1 as usize]
                                * expansion_coeffs[2][v1 as usize]
                                * expansion_coeffs[3][t2 as usize]
                                * expansion_coeffs[4][u2 as usize]
                                * expansion_coeffs[5][v2 as usize]
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
