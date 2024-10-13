use itertools::Itertools;
use nalgebra::{DMatrix, Vector3};
use smallvec::SmallVec;

use crate::{
    basis::ContractedGaussian,
    system::{MolecularSystem, ShellBasis},
};

// TODO: improve cache friendlyness of this struct.
//  The access pattern is:
//  Accesses are localized in k (see parameters of `ExpansionCoefficients::coeffificient`)
//  So we should be storing something that acts like DMatrix<Vec<f64>> instead.
#[derive(Debug)]
pub struct ExpansionCoefficients {
    // Note: Though the 6 here is the same magic value as in the exponents in `ContractedGaussian`,
    //  here this value was chosen for a different reason.
    //  In `ContractedGaussian` it was chosen because 6 is a typical value for the largest
    //  contraction degree, and here it was chosen because the length of these lists is equal to
    //  the sum of angular momentum magnitudes on certain axes.
    //  6 thus seems like a reasonable choice, because anything with higher angular momentum than d
    //  orbitals is rather rare, so it can be a bit slower for these cases
    expansion_x: SmallVec<[DMatrix<f64>; 6]>,
    expansion_y: SmallVec<[DMatrix<f64>; 6]>,
    expansion_z: SmallVec<[DMatrix<f64>; 6]>,
}

impl ExpansionCoefficients {
    pub fn from_basis_pair(
        &ContractedGaussian {
            exponents: ref a,
            angular: [l1, m1, n1],
            ..
        }: &ContractedGaussian,
        &ContractedGaussian {
            exponents: ref b,
            angular: [l2, m2, n2],
            ..
        }: &ContractedGaussian,
        diff: Vector3<f64>,
    ) -> Self {
        let p = DMatrix::from_fn(a.len(), b.len(), |i, j| a[i] + b[j]);
        let q = DMatrix::from_fn(a.len(), b.len(), |i, j| a[i] * b[j] / (a[i] + b[j]));

        let [expansion_x, expansion_y, expansion_z] =
            [(l1, l2, diff.x), (m1, m2, diff.y), (n1, n2, diff.z)].map(|(i, j, diff)| {
                (0..=i + j)
                    .map(|t| hermite_expansion_simd([i, j, t], diff, a, b, &p, &q))
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

/// Essentially a symmetric matrix of expansion coefficients for pairs of basis functions
pub struct HermiteCache {
    data: Vec<ExpansionCoefficients>,
    n: usize,
}

impl HermiteCache {
    pub fn new(system: &MolecularSystem) -> Self {
        let n_basis = system.n_basis();
        let n_shells = system.shells.len();

        let mut data: Vec<_> = (0..n_basis * n_basis).map(|_| None).collect();

        for a in 0..n_shells {
            for b in 0..n_shells {
                let ShellBasis {
                    start_index: start_a,
                    count: count_a,
                    center: pos_a,
                    ..
                } = system.shell_basis(a);

                let ShellBasis {
                    start_index: start_b,
                    count: count_b,
                    center: pos_b,
                    ..
                } = system.shell_basis(b);

                let diff_ab = pos_b - pos_a;

                for i in start_a..start_a + count_a {
                    for j in start_b..start_b + count_b {
                        data[n_basis * i + j] = Some(ExpansionCoefficients::from_basis_pair(
                            system.basis[i],
                            system.basis[j],
                            diff_ab,
                        ));
                    }
                }
            }
        }

        Self {
            data: data.into_iter().map(Option::unwrap).collect_vec(),
            n: n_basis,
        }
    }

    pub fn basis_pair(&self, i: usize, j: usize) -> &ExpansionCoefficients {
        let linear = self.n * i + j;
        &self.data[linear]
    }
}

fn hermite_expansion_simd(
    [i, j, t]: [i32; 3],
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
            let im1tm1 = hermite_expansion_simd([i - 1, j, t - 1], dist, a, b, p, q);
            let im1t = hermite_expansion_simd([i - 1, j, t], dist, a, b, p, q);
            let im1tp1 = hermite_expansion_simd([i - 1, j, t + 1], dist, a, b, p, q);

            p.map(|p| (2.0 * p).recip()).component_mul(&im1tm1)
                - q.map_with_location(|i, _, q| q * dist / a[i])
                    .component_mul(&im1t)
                + (t + 1) as f64 * im1tp1
        }
        (_, _, _) => {
            let jm1tm1 = hermite_expansion_simd([i, j - 1, t - 1], dist, a, b, p, q);
            let jm1t = hermite_expansion_simd([i, j - 1, t], dist, a, b, p, q);
            let jm1tp1 = hermite_expansion_simd([i, j - 1, t + 1], dist, a, b, p, q);

            p.map(|p| (2.0 * p).recip()).component_mul(&jm1tm1)
                + q.map_with_location(|_, j, q| q * dist / b[j])
                    .component_mul(&jm1t)
                + (t + 1) as f64 * jm1tp1
        }
    }
}
