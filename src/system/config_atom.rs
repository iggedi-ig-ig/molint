use nalgebra::Point3;
use serde::{Deserialize, Serialize};

use super::Atom;

/// A helper type for serialization of [Atom]s
#[derive(Debug, Serialize, Deserialize)]
pub(super) struct ConfigAtom {
    element: String,
    position: [f64; 3],
}

impl TryFrom<ConfigAtom> for Atom {
    type Error = anyhow::Error;

    fn try_from(value: ConfigAtom) -> Result<Self, Self::Error> {
        Ok(Self {
            ordinal: value.element.parse()?,
            position: Point3::from(value.position),
        })
    }
}
