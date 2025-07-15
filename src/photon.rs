//! src/photon.rs
//!
//! Represents a single or macro photon from the driver beam,
//! along with useful methods for interaction.

use crate::utils::{FourVec, ThreeVec};

pub struct Photon {
    pub position: FourVec, // four-position (ct, x, y, z)
    pub momentum: FourVec, // four-momentum (E, c px, c py, c pz)
    pub weight: f64,       // weight of the photon
}

impl Photon {
    pub fn new(position: FourVec, momentum: FourVec, weight: f64) -> Self {
        Photon {
            position,
            momentum,
            weight,
        }
    }

    pub fn mandelstam_s(self, other: &Photon) -> f64 {
        // Mandelstam variable s = (p1 + p2)^2
        let p1 = self.momentum;
        let p2 = other.momentum;
        (p1 + p2).norm_squared()
    }
}

pub struct LorentzBooster {}
impl LorentzBooster {
    pub fn new() -> Self {
        LorentzBooster {}
    }

    pub fn boost(&self, photon: &Photon, velocity: ThreeVec) -> Photon {
        // Apply Lorentz boost to the photon
        let gamma = 1.0 / (1.0 - velocity.norm_squared().sqrt()).sqrt();
        let boosted_momentum = photon.momentum * gamma; // Simplified for demonstration
        Photon {
            position: photon.position,
            momentum: boosted_momentum,
            weight: photon.weight,
        }
    }
}
