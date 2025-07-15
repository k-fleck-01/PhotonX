//! src/utils.rs
//!
//! Utility functions and constants for the simulation.

use nalgebra::{Vector3, Vector4};

// Define type aliases for vectors
pub type ThreeVec = Vector3<f64>;
pub type FourVec = Vector4<f64>;

pub mod constants {
    pub const PI: f64 = std::f64::consts::PI;
    pub const PI_2: f64 = std::f64::consts::PI / 2.0;
    pub const TWOPI: f64 = 2.0 * std::f64::consts::PI;

    pub const MC2: f64 = physical_constants::ELECTRON_MASS_ENERGY_EQUIVALENT_IN_MEV; // in MeV
    pub const RE: f64 = physical_constants::CLASSICAL_ELECTRON_RADIUS; // in m
}
