use nalgebra::DMatrix;

use super::*;

/// Represents a symmetric matrix, and efficiently stores it by only storing it's upper triangle
pub struct SymmetricMatrix {
    data: Vec<f64>,
    n: usize,
}

impl SymmetricMatrix {
    /// Create and allocate a [SymmetricMatrix] filled with zeros, of the given size.
    pub(crate) fn zeros(n: usize) -> Self {
        Self {
            data: vec![0.0; n * (n + 1) / 2],
            n,
        }
    }

    /// Copy all entires from the given [DMatrix] into the block starting at
    /// (start_a, start_b) containing exactly (count_a, count_b) elements in the respective
    /// dimensions, while respecting symmetry.
    pub fn copy_from(
        &mut self,
        from: &DMatrix<f64>,
        (start_a, start_b): (usize, usize),
        (count_a, count_b): (usize, usize),
    ) {
        for (i, a) in (start_a..start_a + count_a).enumerate() {
            for (j, b) in (start_b..start_b + count_b)
                .enumerate()
                .skip_while(|&(_, b)| a > b)
            {
                *self.index_unchecked_mut((a, b)) = from[(i, j)];
            }
        }
    }

    fn index_unchecked_mut(&mut self, index: (usize, usize)) -> &mut f64 {
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
