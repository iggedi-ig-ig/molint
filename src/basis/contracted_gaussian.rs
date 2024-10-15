use smallvec::SmallVec;

/// Represents a single contracted gaussian basis function.
// TODO(style): 6 is kind of a magic number here. I think I'm probably fine with this as an
// arbitary amount of exponents / coefficients can still be used, but it's worth thinking about
// again at some later time
#[derive(Clone, Debug, PartialEq, PartialOrd)]
pub struct ContractedGaussian {
    pub coefficients: SmallVec<[f64; 6]>,
    pub exponents: SmallVec<[f64; 6]>,
    pub angular: [i32; 3],
}

impl ContractedGaussian {
    /// Returns an iterator over all coefficients with their corresponding exponents
    pub(crate) fn iter(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.coefficients
            .iter()
            .copied()
            .zip(self.exponents.iter().copied())
    }
}
