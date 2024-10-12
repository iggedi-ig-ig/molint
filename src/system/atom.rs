use nalgebra::Point3;

#[derive(Copy, Clone, Debug)]
pub struct Atom {
    /// also the charge of the nucleus
    pub ordinal: usize,
    pub position: Point3<f64>,
}
