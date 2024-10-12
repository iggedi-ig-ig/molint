mod eri_tensor;
mod symmetric_matrix;

pub use eri_tensor::EriTensor;
pub use symmetric_matrix::SymmetricMatrix;

pub(super) const fn linearize_upper_triangular(n: usize, (i, j): (usize, usize)) -> usize {
    n * i + j - i * (i + 1) / 2
}

pub(super) const fn linearize_symmetric_4d(
    n: usize,
    (i, j, k, l): (usize, usize, usize, usize),
) -> usize {
    let block_index_ij = linearize_upper_triangular(n, (i, j));
    let block_index_kl = linearize_upper_triangular(n, (k, l));

    block_index_ij * n * (n + 1) / 2 + block_index_kl
}

/// if necessary, permute (i, j) such that
///  1. i <= j
pub(super) const fn canonicalize_2d_index((i, j): (usize, usize)) -> (usize, usize) {
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
pub(super) const fn canonicalize_4d_index(
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
