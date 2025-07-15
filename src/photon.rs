//! src/photon.rs
//!
//! Represents a single or macro photon from the driver beam,
//! along with useful methods for interaction.

use crate::utils::{FourVec, ThreeVec};
use nalgebra::{vector, zero};

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

    pub fn mandelstam_s(&self, other: &Photon) -> f64 {
        // Mandelstam variable s = (p1 + p2)^2
        let p_sum = &self.momentum + &other.momentum;
        p_sum.norm_squared()
    }
}

pub struct LorentzBooster {
    pub beta: ThreeVec, // velocity vector (v/c)
}

impl LorentzBooster {
    pub fn new() -> Self {
        LorentzBooster {
            beta: vector![0.0, 0.0, 1.0],
        }
    }

    pub fn calculate_boost_factor(&mut self, p: &FourVec) {
        // Calculate the boost factor from the total four-momentum in the
        // lab frame.
        self.beta = vector![p[1] / p[0], p[2] / p[0], p[3] / p[0]];
    }

    pub fn boost(&self, v: &FourVec) -> FourVec {
        let gamma = 1.0 / (1.0 - self.beta.norm_squared()).sqrt();
        let spatial = vector![v[1], v[2], v[3]];
        let beta_dot_x = self.beta.dot(&spatial);
        let t_prime = gamma * (v[0] - beta_dot_x);
        let spatial_prime = &spatial
            + ((gamma - 1.0) * beta_dot_x / self.beta.norm_squared() - gamma * v[0]) * &self.beta;

        FourVec::new(
            t_prime,
            spatial_prime[0],
            spatial_prime[1],
            spatial_prime[2],
        )
    }

    pub fn inverse_boost(&self, v: &FourVec) -> FourVec {
        let inv = LorentzBooster {
            beta: -self.beta.clone(),
        };
        inv.boost(v)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_boost() {
        let p1 = FourVec::new(10.0, 0.0, 0.0, 10.0);
        let p2 = FourVec::new(20.0, 0.0, 0.0, -20.0);
        let total = p1 + p2;

        let mut booster = LorentzBooster::new();
        booster.calculate_boost_factor(&total);

        // In ZMF, both particles should have the same energy
        // and opposite momenta
        let boosted_p1 = booster.boost(&p1);
        let boosted_p2 = booster.boost(&p2);

        assert!(
            (boosted_p1[0] - boosted_p2[0]).abs() < 1e-6,
            "Energies should match"
        );
        assert!(
            (boosted_p1[1] + boosted_p2[1]).abs() < 1e-6,
            "Momenta should cancel out"
        );
        assert!(
            (boosted_p1[2] + boosted_p2[2]).abs() < 1e-6,
            "Momenta should cancel out"
        );
        assert!(
            (boosted_p1[3] + boosted_p2[3]).abs() < 1e-6,
            "Momenta should cancel out"
        );
    }
}
