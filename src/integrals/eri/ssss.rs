use nalgebra::{Point3, Vector3};
use ndarray::Array4;

use crate::{basis::ContractedGaussian, system::ShellBasis};

/// (SS|SS) eri
pub(super) fn ssss_eri(
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

    // TODO(perf): symmetry
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

                    result[(i, j, k, l)] = contracted_gaussian_eri(
                        [a, b, c, d],
                        [pos_a, pos_b, pos_c, pos_d],
                        [diff_ab, diff_cd],
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
    [diff_ab, diff_cd]: [Vector3<f64>; 2],
) -> f64 {
    let mut sum = 0.0;

    for (coeff_a, exp_a) in a.iter() {
        for (coeff_b, exp_b) in b.iter() {
            let p1 = exp_a + exp_b;
            let q1 = exp_a * exp_b / p1;

            // inlined utils::product_center to reuse p
            let product_center_ab = (exp_a * pos_a.coords + exp_b * pos_b.coords) / p1;

            for (coeff_c, exp_c) in c.iter() {
                for (coeff_d, exp_d) in d.iter() {
                    let p2 = exp_c + exp_d;
                    let q2 = exp_c * exp_d / p2;

                    // inlined utils::product_center to reuse q
                    let product_center_cd = (exp_c * pos_c.coords + exp_d * pos_d.coords) / p2;
                    let diff_product = product_center_cd - product_center_ab;
                    let alpha = p1 * p2 / (p1 + p2);

                    sum += coeff_a
                        * coeff_b
                        * coeff_c
                        * coeff_d
                        * f64::exp(-q1 * diff_ab.norm_squared())
                        * f64::exp(-q2 * diff_cd.norm_squared())
                        * boys::micb25::boys(0, alpha * diff_product.norm_squared())
                        * 2.0
                        * std::f64::consts::PI.powi(5).sqrt()
                        * (p1 * p2 * (p1 + p2).sqrt()).recip()
                }
            }
        }
    }
    sum
}
