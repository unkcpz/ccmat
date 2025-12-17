use ccmat::{lattice_angstrom, sites_frac_coord, CrystalBuilder};

fn main() {
    let lattice = lattice_angstrom![
        a = (1.0, 0.0, 0.0),
        b = (0.0, 1.0, 0.0),
        c = (0.0, 0.0, 1.0),
    ];
    let sites = sites_frac_coord![
        (0.0, 0.0, 0.0), 8;
        (0.0, 0.0, 0.5), 8;
    ];
    let crystal = CrystalBuilder::new()
        .with_lattice(&lattice)
        .with_sites(sites)
        .build()
        .unwrap();

    dbg!(crystal);
}
