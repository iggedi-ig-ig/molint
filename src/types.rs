use std::ops::Index;

#[derive(Copy, Clone, Debug, Eq, PartialEq, Hash)]
struct IntegralIndex(usize, usize, usize, usize);

impl IntegralIndex {
    /// Creates a new integral index with the given indices.
    pub(crate) const fn new(index: (usize, usize, usize, usize)) -> Self {
        let (i, j, k, l) = Self::correct_order(index);
        Self(i, j, k, l)
    }

    /// Returns the indices with the correct order, such that xy <= zw.
    #[inline(always)]
    const fn correct_order(
        (i, j, k, l): (usize, usize, usize, usize),
    ) -> (usize, usize, usize, usize) {
        let (i, j) = if i < j { (i, j) } else { (j, i) };
        let (k, l) = if k < l { (k, l) } else { (l, k) };

        let ij = i * (i + 1) / 2 + j;
        let kl = k * (k + 1) / 2 + l;

        if ij < kl {
            (i, j, k, l)
        } else {
            (k, l, i, j)
        }
    }

    pub fn linear(&self, size: usize) -> usize {
        let &Self(i, j, k, l) = self;
        l * size.pow(3) + k * size.pow(2) + j * size + i
    }
}

pub struct ElectronTensor {
    data: Vec<f64>,
    size: usize,
}

impl Index<(usize, usize, usize, usize)> for ElectronTensor {
    type Output = f64;

    fn index(&self, index: (usize, usize, usize, usize)) -> &Self::Output {
        let index = IntegralIndex::new(index);
        let linear = index.linear(self.size);
        &self.data[linear]
    }
}

impl Index<IntegralIndex> for ElectronTensor {
    type Output = f64;

    fn index(&self, index: IntegralIndex) -> &Self::Output {
        let linear = index.linear(self.size);
        &self.data[linear]
    }
}
