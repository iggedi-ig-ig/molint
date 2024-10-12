use nalgebra::DMatrix;

use super::*;

pub struct SymmetricMatrix {
    data: Vec<f64>,
    n: usize,
}

impl SymmetricMatrix {
    pub(crate) fn zeros(n: usize) -> Self {
        Self {
            data: vec![0.0; n * (n + 1) / 2],
            n,
        }
    }

    pub(crate) fn index_unchecked_mut(&mut self, index: (usize, usize)) -> &mut f64 {
        &mut self.data[linearize_upper_triangular(self.n, index)]
    }
}

impl std::ops::Index<(usize, usize)> for SymmetricMatrix {
    type Output = f64;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        let index = canonicalize_2d_index(index);
        &self.data[linearize_upper_triangular(self.n, index)]
    }
}

impl std::ops::IndexMut<(usize, usize)> for SymmetricMatrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        let index = canonicalize_2d_index(index);
        &mut self.data[linearize_upper_triangular(self.n, index)]
    }
}

impl From<SymmetricMatrix> for DMatrix<f64> {
    fn from(value: SymmetricMatrix) -> Self {
        DMatrix::from(&value)
    }
}

impl From<&SymmetricMatrix> for DMatrix<f64> {
    fn from(value: &SymmetricMatrix) -> Self {
        DMatrix::from_fn(value.n, value.n, |i, j| value[(i, j)])
    }
}
