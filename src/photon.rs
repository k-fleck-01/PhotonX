//! src/photon.rs
//!
//! Represents a single or macro photon from the driver beam,
//! along with useful methods for interaction.

use crate::utils::{FourVec, MinkowskiMetric, ThreeVec};
use nalgebra::vector;

#[derive(Debug, Clone)]
pub struct Photon {
    pub position: FourVec, // four-position (ct, x, y, z)
    pub momentum: FourVec, // four-momentum (E, c px, c py, c pz)
    pub weight: f64,       // weight of the photon
}

impl Photon {
    pub fn new(position: FourVec, momentum: FourVec, weight: f64) -> Self {
        // Check that photon momentum is null
        let k2 = momentum.minkowski_sq_norm();
        if k2.abs() > 1.0e-6 {
            panic!("Photon momentum is not null: k^2 = {}", k2);
        }

        Photon {
            position,
            momentum,
            weight,
        }
    }

    pub fn mandelstam_s(&self, other: &Photon) -> f64 {
        // Mandelstam variable s = (p1 + p2)^2
        let p_sum = &self.momentum + &other.momentum;
        p_sum.minkowski_sq_norm()
    }

    pub fn direction(&self) -> ThreeVec {
        // Return a unit vector in the direction of photon's three-momentum.
        vector![self.momentum[1], self.momentum[2], self.momentum[3]].normalize()
    }

    pub fn angle_with(&self, other: &Photon) -> f64 {
        // Calculate the angle between two photons using their directions
        let dir1 = self.direction();
        let dir2 = other.direction();
        dir1.dot(&dir2).acos()
    }
}

#[derive(Debug, Clone, Default)]
pub struct LorentzBooster {
    /// Allows for Lorentz boosting to the ZMF.
    pub beta: ThreeVec, // velocity vector (v/c)
}

impl LorentzBooster {
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

    pub fn boost_to_zmf(&mut self, p1: &Photon, p2: &Photon) -> (Photon, Photon) {
        // Boost two photons to the zero momentum frame (ZMF)
        let total = &p1.momentum + &p2.momentum;
        self.calculate_boost_factor(&total);

        let boosted_p1 = self.boost(&p1.momentum);
        let boosted_x1 = self.boost(&p1.position);

        let boosted_p2 = self.boost(&p2.momentum);
        let boosted_x2 = self.boost(&p2.position);

        (
            Photon::new(boosted_x1, boosted_p1, p1.weight),
            Photon::new(boosted_x2, boosted_p2, p2.weight),
        )
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

        let mut booster = LorentzBooster::default();
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

    #[test]
    fn test_photon() {
        let photon1 = Photon::new(
            FourVec::new(1.0, 0.0, 0.0, 0.0),
            FourVec::new(1.0, 0.0, 0.0, 1.0),
            1.0,
        );
        let photon2 = Photon::new(
            FourVec::new(1.0, 0.0, 0.0, 0.0),
            FourVec::new(1.0, 0.0, 0.0, -1.0),
            1.0,
        );

        let s = photon1.mandelstam_s(&photon2);
        assert!(s - 4.0 < 1e-6, "Mandelstam s should be 4.0");
    }

    #[test]
    #[should_panic(expected = "Photon momentum is not null")]
    fn test_photon_invalid() {
        // Check panic on non-null momentum
        Photon::new(
            FourVec::new(1.0, 0.0, 0.0, 0.0),
            FourVec::new(1.0, 0.0, 0.0, 2.0), // Non-null momentum
            1.0,
        );
    }

    #[test]
    fn test_angle() {
        let photon1 = Photon::new(
            FourVec::new(1.0, 0.0, 0.0, 0.0),
            FourVec::new(1.0, 0.0, 0.0, 1.0),
            1.0,
        );

        let photon2 = Photon::new(
            FourVec::new(1.0, 0.0, 0.0, 0.0),
            FourVec::new(1.0, 0.0, 0.0, -1.0),
            1.0,
        );

        let theta = photon1.angle_with(&photon2);
        assert!(
            (theta - std::f64::consts::PI).abs() < 1e-6,
            "Angle should be PI radians"
        );

        let photon3 = Photon::new(
            FourVec::new(1.0, 0.0, 0.0, 0.0),
            FourVec::new(1.0, 1.0, 0.0, 0.0),
            1.0,
        );

        let theta = photon1.angle_with(&photon3);

        assert!(
            (theta - std::f64::consts::FRAC_PI_2).abs() < 1e-6,
            "Angle should be PI/2 radians"
        );
    }
}
