use crate::utils::FourVec;

pub mod breit_wheeler;
pub mod photon_scatter;

/// Polarization types for cross sections
#[derive(Debug, Clone, Copy)]
pub enum Polarization {
    Unpolarized,
    Linear,
    Circular,
}

/// Trait to define the cross section interface for processes
pub trait CrossSection {
    fn total(&self, s: f64) -> f64;
    fn differential(&self, s: f64, theta: f64) -> f64;
    fn sample_theta(&self, s: f64, rng: &mut impl rand::Rng) -> f64;
    fn sample_outgoing(&self, s: f64, rng: &mut impl rand::Rng) -> FourVec;
}
