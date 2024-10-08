use nalgebra::DMatrix;

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

pub struct EriTensor {
    data: Vec<f64>,
    n: usize,
}

impl EriTensor {
    pub(crate) fn zeros(n: usize) -> Self {
        // This capacity is an upper bound, as there is no closed form expression for the exact
        // number of integrals to compute.
        // It could technically be precomputed, but this version allows for simplified index math.
        Self {
            data: vec![0.0; n.pow(2) * (n + 1).pow(2) / 4],
            n,
        }
    }

    pub(crate) fn index_unchecked_mut(&mut self, index: (usize, usize, usize, usize)) -> &mut f64 {
        &mut self.data[linearize_symmetric_4d(self.n, index)]
    }
}

impl std::ops::Index<(usize, usize, usize, usize)> for EriTensor {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize, usize)) -> &Self::Output {
        let index = canonicalize_4d_index(index);
        &self.data[linearize_symmetric_4d(self.n, index)]
    }
}

impl std::ops::IndexMut<(usize, usize, usize, usize)> for EriTensor {
    fn index_mut(&mut self, index: (usize, usize, usize, usize)) -> &mut Self::Output {
        let index = canonicalize_4d_index(index);
        &mut self.data[linearize_symmetric_4d(self.n, index)]
    }
}

const fn linearize_upper_triangular(n: usize, (i, j): (usize, usize)) -> usize {
    n * i + j - i * (i + 1) / 2
}

const fn linearize_symmetric_4d(n: usize, (i, j, k, l): (usize, usize, usize, usize)) -> usize {
    let block_index_ij = linearize_upper_triangular(n, (i, j));
    let block_index_kl = linearize_upper_triangular(n, (k, l));

    block_index_ij * n * (n + 1) / 2 + block_index_kl
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
