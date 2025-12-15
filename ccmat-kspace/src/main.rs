use ccmat_core::{atomic_number, lattice_angstrom, sites_frac_coord, BravaisClass, CrystalBuilder};
use ccmat_kspace::find_path;
use ccmat_symmetry::analyze_symmetry;

fn main() {
    let lattice = lattice_angstrom![(4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0),];

    let sites = sites_frac_coord![
        (0.0, 0.0, 0.0), atomic_number!(C);
        (0.5, 0.5, 0.5), atomic_number!(C);
        (0.0, 0.0, 0.1), atomic_number!(O);
        (0.5, 0.5, 0.6), atomic_number!(O);
    ];

    let s = CrystalBuilder::new()
        .with_lattice(&lattice)
        .with_sites(&sites)
        .build()
        .unwrap();

    let syminfo = analyze_symmetry(&s, 1e-5).unwrap();
    assert_eq!(syminfo.bravais_class(), BravaisClass::tI);

    let _ = find_path(&s, 1e-5, 1e-7);
}
