//! cross_sections/breit_wheeler.rs
//!
//! Cross section for the Breit-Wheeler process.

use crate::utils::constants::{PI_2, RE};
use std::f64;

/// Total cross section for the Breit-Wheeler process.
pub fn total_sigma(s: f64) -> f64 {
    // Check above threshold
    if s < 4.0 {
        return 0.0;
    }

    let beta = (1.0 - 4.0 / s).sqrt();
    let larg = (1.0 + beta) / (1.0 - beta);
    let log = larg.ln();

    let term1 = (3.0 - beta.powi(4)) * log;
    let term2 = 2.0 * beta * (2.0 - beta.powi(2));
    let sigma = PI_2 * RE * RE * (term1 - term2) * (1.0 - beta.powi(2));
    sigma
}

/// Differential cross section for the Breit-Wheeler process,
/// do/dO, in the ZMF.
pub fn diff_sigma(s: f64, theta: f64) -> f64 {
    if s < 4.0 {
        return 0.0;
    }

    let beta = (1.0 - 4.0 / s).sqrt();
    let sin2t = f64::sin(theta).powi(2);
    let cos2t = 1.0 - sin2t;

    let term1 = 1.0 + 2.0 * beta * sin2t - beta.powi(4) - (beta.powi(2) * sin2t).powi(2);
    let term2 = (1.0 - beta.powi(2) * cos2t).powi(2);
    let dsigma = (RE * RE * beta / s) * (term1 / term2);
    dsigma
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_total_sigma() {
        // Check below threshold
        let mut s = 2.0;
        let mut sigma = total_sigma(s);
        assert_eq!(sigma, 0.0);

        // Check above threshold
        s = 5.0;
        sigma = total_sigma(s);
        assert!(sigma > 0.0);
    }
}
