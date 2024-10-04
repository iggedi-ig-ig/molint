use std::{collections::HashMap, fs::File, path::Path};

use serde::Deserialize;

use crate::{periodic_table::ElementType, system::Atom};

/// Data associated with a contracted gaussian. stored as a struct of lists.
// TODO(perf): SmallVec (or even complete stack storage?)
#[derive(Clone, Debug, PartialEq, PartialOrd)]
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

impl BasisSet {
    pub fn load(path: impl AsRef<Path>) -> anyhow::Result<Self> {
        let basis_set: BseBasisSet = serde_json::from_reader(File::open(path)?)?;
        basis_set.try_into()
    }

    pub fn atomic_basis(&self, atom: &Atom) -> &[ContractedGaussian] {
        let element_type = ElementType::from_ordinal(atom.ordinal)
            .unwrap_or_else(|| panic!("failed to convert ordinal {} to ElementType", atom.ordinal));
        &self.0[&element_type]
    }
}

#[derive(Deserialize, Debug)]
struct BseBasisSet {
    elements: HashMap<ElementType, BseElectronicConfiguration>,
}

#[derive(Deserialize, Debug)]
struct BseElectronicConfiguration {
    electron_shells: Vec<BseElectronShell>,
}

#[derive(Deserialize, Debug)]
#[allow(unused)]
struct BseElectronShell {
    function_type: String,
    angular_momentum: Vec<i32>,
    exponents: Vec<String>,
    coefficients: Vec<Vec<String>>,
}

impl TryFrom<BseBasisSet> for BasisSet {
    // TODO: use a "better" type for error
    type Error = anyhow::Error;

    fn try_from(value: BseBasisSet) -> Result<Self, Self::Error> {
        let mut atomic_mapping = HashMap::with_capacity(value.elements.len());

        // TODO: this is pretty deeply nested, this can definitely be improved somehow
        for (element, configuration) in value.elements {
            let mut element_basis = Vec::new();

            for electron_shell in &configuration.electron_shells {
                if electron_shell.function_type != "gto" {
                    log::warn!("skipping function type {}", electron_shell.function_type);
                    continue;
                }

                for (index, &angular_magnitude) in
                    electron_shell.angular_momentum.iter().enumerate()
                {
                    let angular_vectors = generate_angular_vectors(angular_magnitude);

                    for angular @ (i, j, k) in angular_vectors {
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

                            let norm = gaussian_norm(exponent, angular);

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

/// Normalization constant of a [GaussianPrimitive] with the given parameters
fn gaussian_norm(exponent: f64, angular: (i32, i32, i32)) -> f64 {
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
