use std::{collections::HashMap, error::Error};

use serde::Deserialize;

use crate::{periodic_table::ElementType, utils};

/// Data associated with a contracted gaussian. stored as a struct of lists.
#[derive(Clone, Debug)]
pub struct ContractedGaussian {
    pub coefficients: Vec<f64>,
    pub exponents: Vec<f64>,
    pub angular: [u32; 3],
}

impl ContractedGaussian {
    pub(crate) fn iter(&self) -> impl Iterator<Item = (f64, f64)> + '_ {
        self.coefficients
            .iter()
            .copied()
            .zip(self.exponents.iter().copied())
    }
}

pub struct BasisSet(HashMap<ElementType, Vec<ContractedGaussian>>);

#[derive(Deserialize)]
pub struct ConfigBasisSet {
    elements: HashMap<ElementType, ConfigElectronicConfiguration>,
}

#[derive(Deserialize)]
struct ConfigElectronicConfiguration {
    electron_shells: Vec<ConfigElectronShell>,
}

#[derive(Deserialize)]
#[allow(unused)]
struct ConfigElectronShell {
    function_type: String,
    angular_momentum: Vec<i32>,
    exponents: Vec<String>,
    coefficients: Vec<Vec<String>>,
}

impl TryFrom<ConfigBasisSet> for BasisSet {
    // TODO: use a "better" type for error
    type Error = Box<dyn Error>;

    fn try_from(value: ConfigBasisSet) -> Result<Self, Self::Error> {
        let mut atomic_mapping = HashMap::with_capacity(value.elements.len());

        // TODO: this is pretty deeply nested, this can definitely be improved somehow
        for (element, configuration) in value.elements {
            let mut element_basis = Vec::new();

            for electron_shell in &configuration.electron_shells {
                for (index, &angular_magnitude) in
                    electron_shell.angular_momentum.iter().enumerate()
                {
                    let angular_vectors = generate_angular_vectors(angular_magnitude);

                    for angular @ (i, j, k) in angular_vectors {
                        assert_eq!(
                            electron_shell.exponents.len(),
                            electron_shell.coefficients.len(),
                            "have to have the same amount of exponents and coefficients"
                        );

                        let mut exponents = Vec::with_capacity(electron_shell.exponents.len());
                        let mut coefficients =
                            Vec::with_capacity(electron_shell.coefficients.len());

                        for (exponent, coefficient) in electron_shell
                            .exponents
                            .iter()
                            .zip(&electron_shell.coefficients[index])
                        {
                            let exponent = exponent.parse::<f64>()?;
                            let coefficient = coefficient.parse::<f64>()?;

                            let norm = utils::gaussian_norm(exponent, angular);

                            exponents.push(exponent);
                            coefficients.push(norm * coefficient);
                        }

                        element_basis.push(ContractedGaussian {
                            coefficients,
                            exponents,
                            angular: [i, j, k].map(|i| i as u32),
                        });
                    }
                }
            }

            atomic_mapping.insert(element, element_basis);
        }

        Ok(Self(atomic_mapping))
    }
}

// generate all (i, j, k) such that i + j + k = angular
fn generate_angular_vectors(angular_magnitude: i32) -> Vec<(i32, i32, i32)> {
    let mut angular_vectors = Vec::with_capacity(8);

    for (i, j, k) in itertools::iproduct!(
        0..=angular_magnitude,
        0..=angular_magnitude,
        0..=angular_magnitude
    ) {
        if i + j + k == angular_magnitude {
            angular_vectors.push((i, j, k));
        }
    }

    angular_vectors
}
