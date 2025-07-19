//! cross_sections/breit_wheeler.rs
//!
//! Cross section for the Breit-Wheeler process.

use super::{CrossSection, Polarization};
use crate::utils::constants::{MC2, RE};
use crate::utils::{FourVec, MinkowskiMetric};
use std::f64;

pub struct BreitWheelerProcess {
    pub pmode: Polarization,
}

impl CrossSection for BreitWheelerProcess {
    fn total(&self, s: f64) -> f64 {
        match self.pmode {
            Polarization::Unpolarized => BreitWheelerProcess::total_unpolarized(s),
            Polarization::Circular => todo!(),
            _ => todo!(),
        }
    }

    fn differential(&self, s: f64, theta: f64) -> f64 {
        match self.pmode {
            Polarization::Unpolarized => BreitWheelerProcess::differential_unpolarized(s, theta),
            Polarization::Circular => unimplemented!(),
            Polarization::Linear => unimplemented!(),
        }
    }

    fn sample_theta(&self, s: f64, rng: &mut impl rand::Rng) -> f64 {
        /// Samples the polar angle theta in the ZMF
        let max_sigma = self.differential(s, f64::consts::FRAC_2_PI);

        loop {
            let theta = rng.gen_range(0.0..f64::consts::PI);
            let y = rng.gen_range(0.0..max_sigma);

            if y < self.differential(s, theta) {
                return theta;
            }
        }
    }

    fn sample_outgoing(&self, s: f64, rng: &mut impl rand::Rng) -> crate::utils::FourVec {
        /// Generates the four-momentum of one of the outgoing leptons in the ZMF.
        let energy = 0.5 * s.sqrt();
        let theta = self.sample_theta(s, rng);
        let phi = rng.gen_range(0.0..(2.0 * f64::consts::PI));

        let pmag = (energy.powi(2) - MC2.powi(2)).sqrt();
        let px = pmag * theta.sin() * phi.cos();
        let py = pmag * theta.sin() * phi.sin();
        let pz = pmag * theta.cos();
        FourVec::new(energy, px, py, pz)
    }
}

impl BreitWheelerProcess {
    pub fn new(pmode: Polarization) -> Self {
        Self { pmode }
    }

    /// Total unpolarized cross section for the Breit-Wheeler process.
    fn total_unpolarized(s: f64) -> f64 {
        // Check above threshold
        if s < 4.0 {
            return 0.0;
        }

        let beta = (1.0 - 4.0 / s).sqrt();
        let larg = (1.0 + beta) / (1.0 - beta);
        let log = larg.ln();

        let term1 = (3.0 - beta.powi(4)) * log;
        let term2 = 2.0 * beta * (2.0 - beta.powi(2));
        f64::consts::FRAC_PI_2 * RE * RE * (term1 - term2) * (1.0 - beta.powi(2))
    }

    /// Differential unpolarized cross section for the Breit-Wheeler process,
    /// do/dO, in the ZMF.
    fn differential_unpolarized(s: f64, theta: f64) -> f64 {
        if s < 4.0 {
            return 0.0;
        }

        let beta = (1.0 - 4.0 / s).sqrt();
        let sin2t = f64::sin(theta).powi(2);
        let cos2t = 1.0 - sin2t;

        let term1 = 1.0 + 2.0 * beta * sin2t - beta.powi(4) - (beta.powi(2) * sin2t).powi(2);
        let term2 = (1.0 - beta.powi(2) * cos2t).powi(2);
        (RE * RE * beta / s) * (term1 / term2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_total_sigma() {
        let pol = Polarization::Unpolarized;
        let bw_process = BreitWheelerProcess::new(pol);

        // Check below threshold
        let mut s = 2.0;
        let mut sigma = bw_process.total(s);
        assert_eq!(sigma, 0.0);

        // Check above threshold
        s = 5.0;
        sigma = bw_process.total(s);
        assert!(sigma > 0.0);
    }

    #[test]
    fn test_sample() {
        let pol = Polarization::Unpolarized;
        let bw_process = BreitWheelerProcess::new(pol);
        let mut rng = rand::thread_rng();

        let s = 5.0; // Above threshold
        let p_out = bw_process.sample_outgoing(s, &mut rng);

        // p_out^2 = m^2
        let p2 = p_out.minkowski_sq_norm();
        assert!(
            (p2 - MC2.powi(2)).abs() < 1e-6,
            "Outgoing particle does not have correct mass"
        );
    }
}
