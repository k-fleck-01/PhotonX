//! src/interaction.rs
//!
//! Controls the interaction of photons, querying if an interaction occurs
//! and generating any necessary secondaries.
//!

use crate::photon::Photon;

pub struct InteractionGenerator {}

impl InteractionGenerator {
    pub fn new() -> Self {
        InteractionGenerator {}
    }

    pub fn does_interact(&self, photon1: &Photon, photon2: &Photon) -> bool {
        // TODO: Add interaction logic
    }
}
