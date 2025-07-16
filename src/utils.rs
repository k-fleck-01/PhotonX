//! src/utils.rs
//!
//! Utility functions and constants for the simulation.

use nalgebra::{Vector3, Vector4};

// Define type aliases for vectors
pub type ThreeVec = Vector3<f64>;
pub type FourVec = Vector4<f64>;

pub trait MinkowskiMetric {
    /// Calculates norms of 4-vectors using the Minkowski metric.
    fn minkowski_dot(&self, other: &Self) -> f64;
    fn minkowski_norm(&self) -> f64;
    fn minkowski_sq_norm(&self) -> f64;
}

impl MinkowskiMetric for FourVec {
    /// Assumes Minkowski metric with signature (+---).
    fn minkowski_dot(&self, other: &Self) -> f64 {
        self[0] * other[0] - self[1] * other[1] - self[2] * other[2] - self[3] * other[3]
    }

    fn minkowski_norm(&self) -> f64 {
        self.minkowski_sq_norm().abs().sqrt()
    }

    fn minkowski_sq_norm(&self) -> f64 {
        self.minkowski_dot(self)
    }
}

pub mod constants {
    pub const PI: f64 = std::f64::consts::PI;
    pub const PI_2: f64 = std::f64::consts::PI / 2.0;
    pub const TWOPI: f64 = 2.0 * std::f64::consts::PI;

    pub const MC2: f64 = physical_constants::ELECTRON_MASS_ENERGY_EQUIVALENT_IN_MEV; // in MeV
    pub const RE: f64 = physical_constants::CLASSICAL_ELECTRON_RADIUS; // in m
}
