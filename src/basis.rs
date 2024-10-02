/// Data associated with a contracted gaussian. stored as a struct of lists.
#[derive(Clone, Debug)]
pub struct ContractedGaussian {
    pub coefficients: Vec<f64>,
    pub exponents: Vec<f64>,
    pub angular: [u32; 3],
}

/// Normalization constant of a [GaussianPrimitive] with the given parameters
pub(crate) fn gaussian_norm(exponent: f64, angular: (i32, i32, i32)) -> f64 {
    let (i, j, k) = angular;

    (std::f64::consts::FRAC_2_PI * exponent)
        .powi(3)
        .sqrt()
        .sqrt()
        * f64::sqrt(
            (8.0 * exponent).powi(i + j + k)
                / ((i + 1..=2 * i).product::<i32>()
                    * (j + 1..=2 * j).product::<i32>()
                    * (k + 1..=2 * k).product::<i32>()) as f64,
        )
}
