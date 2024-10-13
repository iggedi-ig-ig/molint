use ndarray::Array4;

use super::*;

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

    pub(crate) fn copy_from(
        &mut self,
        from: &Array4<f64>,
        (start_a, start_b, start_c, start_d): (usize, usize, usize, usize),
        (count_a, count_b, count_c, count_d): (usize, usize, usize, usize),
    ) {
        for (i, a) in (start_a..start_a + count_a).enumerate() {
            for (j, b) in (start_b..start_b + count_b)
                .enumerate()
                .skip_while(|&(_, b)| a > b)
            {
                let ab = a * (a + 1) / 2 + b;
                for (k, c) in (start_c..start_c + count_c).enumerate() {
                    for (l, d) in (start_d..start_d + count_d)
                        .enumerate()
                        .skip_while(|&(_, d)| c > d || ab > c * (c + 1) / 2 + d)
                    {
                        *self.index_unchecked_mut((a, b, c, d)) = from[(i, j, k, l)];
                    }
                }
            }
        }
    }

    fn index_unchecked_mut(&mut self, index: (usize, usize, usize, usize)) -> &mut f64 {
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
