// NOTE: trybuild didn't able to fail with compile error, same for rust-analyzer, why??
use ccmat::{lattice_angstrom, sites_frac_coord, CrystalBuilder};

fn main() {
    let lattice = lattice_angstrom![
        a = (1.0, 0.0, 0.0),
        b = (0.0, 1.0, 0.0),
        c = (0.0, 0.0, 1.0),
    ];
    let sites = sites_frac_coord![
        (0.0, 0.0, 0.0), 8;
        (0.0, 0.0, 1.5), 8; // !error here, 1.5 out of range
    ];
    let _ = CrystalBuilder::new()
        .with_lattice(&lattice)
        .with_frac_sites(&sites)
        .build()
        .unwrap();
}
