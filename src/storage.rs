use nalgebra::DMatrix;
use ndarray::Array4;

pub struct SymmetricMatrix(DMatrix<f64>);

impl SymmetricMatrix {
    pub fn zeros(n: usize) -> Self {
        Self(DMatrix::zeros(n, n))
    }

    /// Index into the underlying [Array2] directly, without canonicalizing the index
    ///
    /// # Safety
    /// The caller must ensure that the index is canonical
    pub unsafe fn index_unchecked_mut(&mut self, index: (usize, usize)) -> &mut f64 {
        &mut self.0[index]
    }

    /// Index into the underlying [Array2] directly, without canonicalizing the index
    ///
    /// # Safety
    /// The caller must ensure that the index is canonical
    pub unsafe fn index_unchecked(&self, index: (usize, usize)) -> f64 {
        self.0[index]
    }
}

impl std::ops::Index<(usize, usize)> for SymmetricMatrix {
    type Output = f64;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.0[canonicalize_2d_index(index)]
    }
}

impl std::ops::IndexMut<(usize, usize)> for SymmetricMatrix {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.0[canonicalize_2d_index(index)]
    }
}

impl std::fmt::Display for SymmetricMatrix {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl From<SymmetricMatrix> for DMatrix<f64> {
    fn from(value: SymmetricMatrix) -> Self {
        DMatrix::from_fn(value.0.nrows(), value.0.ncols(), |i, j| value[(i, j)])
    }
}

pub struct EriTensor(Array4<f64>);

impl EriTensor {
    pub fn zeros(n: usize) -> Self {
        Self(Array4::zeros([n; 4]))
    }

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

impl std::ops::IndexMut<(usize, usize, usize, usize)> for EriTensor {
    fn index_mut(&mut self, index: (usize, usize, usize, usize)) -> &mut Self::Output {
        &mut self.0[canonicalize_4d_index(index)]
    }
}

impl std::fmt::Display for EriTensor {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
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

#[cfg(test)]
mod tests {}
