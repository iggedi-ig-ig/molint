use std::{collections::HashMap, fs::File, path::Path};

use serde::Deserialize;
use smallvec::SmallVec;

use crate::{periodic_table::ElementType, system::Atom};

/// Data associated with a contracted gaussian. stored as a struct of lists.
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
    // TODO(style): use a "better" type for error
    type Error = anyhow::Error;

    fn try_from(value: BseBasisSet) -> Result<Self, Self::Error> {
        let mut atomic_mapping = HashMap::with_capacity(value.elements.len());

        // TODO(style): this is pretty deeply nested, this can definitely be improved somehow
        for (element, configuration) in value.elements {
            let mut element_basis = Vec::new();

            for electron_shell in &configuration.electron_shells {
                match electron_shell.function_type.as_str() {
                    "gto" | "gto_cartesian" => {
                        // these are equivalent and can be computed using the current integral
                        // implementation. Angular momentum is represented as cartesian
                        // polynomials.
                    }
                    "gto_spherical" => {
                        // In contrast to "gto" and "gto_cartesian", basis functions are
                        // (supposed to be) represented as spherical harmonics. The current
                        // implementation, however, does not support integrating over these.

                        // TODO(completeness): implement conversion function from "gto_spherical"
                        //  to "gto_cartesian"
                        log::warn!(
                            r#"skipping unsupported basis function type "gto_spherical" on element {element:?}"#
                        );
                        continue;
                    }
                    function_type => {
                        log::warn!("skipping unknown basis function type {function_type} on element {element:?}");
                        continue;
                    }
                }

                for (index, &angular_magnitude) in
                    electron_shell.angular_momentum.iter().enumerate()
                {
                    let angular_vectors = generate_angular_vectors(angular_magnitude);

                    for angular in angular_vectors {
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
                            coefficients: coefficients.into(),
                            exponents: exponents.into(),
                            angular,
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
fn generate_angular_vectors(angular_magnitude: i32) -> Vec<[i32; 3]> {
    let mut angular_vectors = Vec::with_capacity(8);

    for (i, j, k) in itertools::iproduct!(
        0..=angular_magnitude,
        0..=angular_magnitude,
        0..=angular_magnitude
    ) {
        if i + j + k == angular_magnitude {
            angular_vectors.push([i, j, k]);
        }
    }

    angular_vectors
}

/// Normalization constant of a [GaussianPrimitive] with the given parameters
fn gaussian_norm(exponent: f64, [i, j, k]: [i32; 3]) -> f64 {
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
