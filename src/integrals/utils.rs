//! Module for helper functions thate are used for different types of integralsc

use nalgebra::{Point3, Vector3};

/// Returns the hermite expansion coefficients as commonly used in molecular integrals
///
/// # References
///
/// [1] Goings, J. Integrals. https://joshuagoings.com/2017/04/28/integrals/
#[allow(clippy::many_single_char_names)]
pub(super) fn hermite_expansion([i, j, t]: [i32; 3], dist: f64, a: f64, b: f64) -> f64 {
    let p = a + b;
    let q = a * b / p;

    if t < 0 || t > i + j {
        0.0
    } else if i == j && j == t && t == 0 {
        f64::exp(-q * dist.powi(2))
    } else if j == 0 {
        (2.0 * p).recip() * hermite_expansion([i - 1, j, t - 1], dist, a, b)
            - (q * dist / a) * hermite_expansion([i - 1, j, t], dist, a, b)
            + (t + 1) as f64 * hermite_expansion([i - 1, j, t + 1], dist, a, b)
    } else {
        (2.0 * p).recip() * hermite_expansion([i, j - 1, t - 1], dist, a, b)
            + (q * dist / b) * hermite_expansion([i, j - 1, t], dist, a, b)
            + (t + 1) as f64 * hermite_expansion([i, j - 1, t + 1], dist, a, b)
    }
}

/// Returns the auxiliary coulomb integral as needed for ERIs
///
/// # References
///
/// [1] Goings, J. Integrals. https://joshuagoings.com/2017/04/28/integrals/
#[allow(clippy::many_single_char_names)]
pub(super) fn coulomb_auxiliary(t: i32, u: i32, v: i32, n: i32, p: f64, diff: Vector3<f64>) -> f64 {
    if t == u && u == v && v == 0 {
        (-2.0 * p).powi(n) * boys::micb25::boys(n as u64, p * diff.norm_squared())
    } else if t == u && u == 0 {
        diff.z * coulomb_auxiliary(t, u, v - 1, n + 1, p, diff)
            + if v > 1 {
                (v - 1) as f64 * coulomb_auxiliary(t, u, v - 2, n + 1, p, diff)
            } else {
                0.0
            }
    } else if t == 0 {
        diff.y * coulomb_auxiliary(t, u - 1, v, n + 1, p, diff)
            + if u > 1 {
                (u - 1) as f64 * coulomb_auxiliary(t, u - 2, v, n + 1, p, diff)
            } else {
                0.0
            }
    } else {
        diff.x * coulomb_auxiliary(t - 1, u, v, n + 1, p, diff)
            + if t > 1 {
                coulomb_auxiliary(t - 2, u, v, n + 1, p, diff)
            } else {
                0.0
            }
    }
}

/// Retruns the product center of two gaussians with the given positions and exponents
pub(super) fn product_center(
    exp_a: f64,
    pos_a: Point3<f64>,
    exp_b: f64,
    pos_b: Point3<f64>,
) -> Point3<f64> {
    Point3::from((exp_a * pos_a.coords + exp_b * pos_b.coords) / (exp_a + exp_b))
}
