use nalgebra::{DMatrix, Vector3};

use crate::basis::ContractedGaussian;

pub struct ExpansionCoefficients {
    expansion_x: Vec<DMatrix<f64>>,
    expansion_y: Vec<DMatrix<f64>>,
    expansion_z: Vec<DMatrix<f64>>,
}

impl ExpansionCoefficients {
    pub fn compute_for(
        basis_a: &ContractedGaussian,
        basis_b: &ContractedGaussian,
        diff: Vector3<f64>,
    ) -> Self {
        let [l1, m1, n1] = basis_a.angular;
        let [l2, m2, n2] = basis_b.angular;

        let a = &basis_a.exponents;
        let b = &basis_b.exponents;

        let p = DMatrix::from_fn(a.len(), b.len(), |i, j| a[i] + b[j]);
        let q = DMatrix::from_fn(a.len(), b.len(), |i, j| a[i] * b[j] / (a[i] + b[j]));
        let [expansion_x, expansion_y, expansion_z] =
            [(l1, l2, diff.x), (m1, m2, diff.y), (n1, n2, diff.z)].map(|(i, j, diff)| {
                (0..=i + j)
                    .map(|t| hermite_expansion_simd(i, j, t, diff, a, b, &p, &q))
                    .collect()
            });

        Self {
            expansion_x,
            expansion_y,
            expansion_z,
        }
    }

    pub fn coefficient(&self, axis: usize, i: usize, j: usize, k: usize) -> f64 {
        let index_into = match axis {
            0 => &self.expansion_x,
            1 => &self.expansion_y,
            2 => &self.expansion_z,
            _ => panic!(),
        };

        index_into[k][(i, j)]
    }
}

fn hermite_expansion_simd(
    i: i32,
    j: i32,
    t: i32,
    dist: f64,
    a: &[f64],
    b: &[f64],
    p: &DMatrix<f64>,
    q: &DMatrix<f64>,
) -> DMatrix<f64> {
    match (i, j, t) {
        (i, j, t) if t < 0 || t > i + j => DMatrix::zeros(a.len(), b.len()),
        (0, 0, 0) => q.map(|q| f64::exp(-q * dist.powi(2))),
        (_, 0, _) => {
            let im1tm1 = hermite_expansion_simd(i - 1, j, t - 1, dist, a, b, p, q);
            let im1t = hermite_expansion_simd(i - 1, j, t, dist, a, b, p, q);
            let im1tp1 = hermite_expansion_simd(i - 1, j, t + 1, dist, a, b, p, q);

            p.map(|p| (2.0 * p).recip()).component_mul(&im1tm1)
                - q.map_with_location(|i, _, q| q * dist / a[i])
                    .component_mul(&im1t)
                + (t + 1) as f64 * im1tp1
        }
        (_, _, _) => {
            let jm1tm1 = hermite_expansion_simd(i, j - 1, t - 1, dist, a, b, p, q);
            let jm1t = hermite_expansion_simd(i, j - 1, t, dist, a, b, p, q);
            let jm1tp1 = hermite_expansion_simd(i, j - 1, t + 1, dist, a, b, p, q);

            p.map(|p| (2.0 * p).recip()).component_mul(&jm1tm1)
                + q.map_with_location(|_, j, q| q * dist / b[j])
                    .component_mul(&jm1t)
                + (t + 1) as f64 * jm1tp1
        }
    }
}
