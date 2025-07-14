### Makefile for PhotonX

### Default target
.PHONY: all
all: build

### Rust run and build
.PHONY: build
build:
	cargo build --release

.PHONY: run
run:
	cargo run --release

### Tabulate photon scattering cross sections
.PHONY: tabulate
tabulate:
	python3 scripts/tabulate_cross_sections.py

### Clean Rust build artifacts
.PHONY: clean
clean:
	cargo clean

### Clean generated data files
.PHONY: clean-data
clean-data:
	rm -f data/*.h5


