use nalgebra::DMatrix;
use ndarray::Array4;

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

/// Typesafe wrapper for a tensor respecting the symmetries of ERIs
pub struct EriTensor(pub(crate) Array4<f64>);

impl EriTensor {
    /// Index into the underlying [Array4] directly, without canonicalizing the index
    ///
    /// # Safety
    /// The caller must ensure that the index is canonical
    pub unsafe fn index_unchecked_mut(&mut self, index: (usize, usize, usize, usize)) -> &mut f64 {
        &mut self.0[index]
    }

    /// Index into the underlying [Array4] directly, without canonicalizing the index
    ///
    /// # Safety
    /// The caller must ensure that the index is canonical
    pub unsafe fn index_unchecked(&self, index: (usize, usize, usize, usize)) -> f64 {
        self.0[index]
    }
}

impl std::ops::Index<(usize, usize, usize, usize)> for EriTensor {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize, usize)) -> &Self::Output {
        &self.0[canonicalize_4d_index(index)]
    }
}

impl std::fmt::Display for EriTensor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

const fn linearize_upper_triangular(n: usize, (i, j): (usize, usize)) -> usize {
    n * i - i * (i + 1) / 2 + j
}

/// if necessary, permute (i, j) such that
///  1. i <= j
const fn canonicalize_2d_index((i, j): (usize, usize)) -> (usize, usize) {
    if i <= j {
        (i, j)
    } else {
        (j, i)
    }
}

/// if necessary, permute (i, j, k, l) such that
///  1. i <= j
///  2. k <= l
///  3. i(i+1)/2+j <= k(k+1)/2+l
const fn canonicalize_4d_index(
    (i, j, k, l): (usize, usize, usize, usize),
) -> (usize, usize, usize, usize) {
    let (i, j) = canonicalize_2d_index((i, j));
    let (k, l) = canonicalize_2d_index((k, l));

    let ij = i * (i + 1) / 2 + j;
    let kl = k * (k + 1) / 2 + l;

    if ij <= kl {
        (i, j, k, l)
    } else {
        (k, l, i, j)
    }
}
