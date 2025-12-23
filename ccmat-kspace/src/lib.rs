mod path;

use ccmat_core::{
    math::{approx_f64, Matrix3, TransformationMatrix, Vector3},
    matrix_3x3, BravaisClass, Crystal, CrystalBuilder, FracCoord, SiteFraction,
};
use ccmat_symmetry::{analyze_symmetry, moyo_wrapper::NiggliReduce, SymmetryInfo};
use tracing::warn;

use crate::path::{KpathEval, KpathInfo};

#[allow(clippy::too_many_lines)]
fn find_p_matrix(syminfo: &SymmetryInfo) -> (Matrix3, Matrix3) {
    let bravais_class = syminfo.bravais_class();
    let spg_number = syminfo.spg_number();

    let (p, inv_p) = match bravais_class {
        BravaisClass::cP
        | BravaisClass::tP
        | BravaisClass::hP
        | BravaisClass::oP
        | BravaisClass::mP
        // for aP, should have already obtained the primitive cell
        | BravaisClass::aP => (
            matrix_3x3![
                1, 0, 0;
                0, 1, 0;
                0, 0, 1;
            ],
            matrix_3x3![
                1, 0, 0;
                0, 1, 0;
                0, 0, 1;
            ],
        ),
        BravaisClass::cF | BravaisClass::oF => (
            (1.0 / 2.0)
                * matrix_3x3![
                    0, 1, 1;
                    1, 0, 1;
                    1, 1, 0;
                ],
            matrix_3x3![
                -1,  1,  1;
                 1, -1,  1;
                 1,  1, -1;
            ],
        ),
        BravaisClass::cI | BravaisClass::tI | BravaisClass::oI => (
            (1.0 / 2.0)
                * matrix_3x3![
                    -1,  1,  1;
                     1, -1,  1;
                     1,  1, -1;
                ],
            matrix_3x3![
                0, 1, 1;
                1, 0, 1;
                1, 1, 0;
            ],
        ),
        BravaisClass::hR => (
            (1.0 / 3.0)
                * matrix_3x3![
                    2, -1, -1;
                    1,  1, -2;
                    1,  1,  1;
                ],
            matrix_3x3![
                 1,  0,  1;
                -1,  1,  1;
                 0, -1,  1;
            ],
        ),
        BravaisClass::oS => match spg_number {
            // oA
            x if (38..=41).contains(&x) => (
                (1.0 / 2.0)
                    * matrix_3x3![
                        0,  0,  2;
                        1,  1,  0;
                       -1,  1,  0;
                    ],
                matrix_3x3![
                    0,  1, -1;
                    0,  1,  1;
                    1,  0,  0;
                ],
            ),
            // oC
            x if (20..=21).contains(&x) || (35..=37).contains(&x) || (63..=68).contains(&x) => (
                (1.0 / 2.0)
                    * matrix_3x3![
                        1,  1,  0;
                       -1,  1,  0;
                        0,  0,  2;
                    ],
                matrix_3x3![
                    1, -1,  0;
                    1,  1,  0;
                    0,  0,  1;
                ],
            ),
            _ => unreachable!("oS bravais lattice spacegroup number in wrong range"),
        },
        BravaisClass::mC => (
            (1.0 / 2.0)
                * matrix_3x3![
                    1, -1,  0;
                    1,  1,  0;
                    0,  0,  2;
                ],
            matrix_3x3![
                 1,  1,  0;
                -1,  1,  0;
                 0,  0,  1;
            ],
        ),
    };

    (p, inv_p)
}

fn find_primitive_hpkot(
    standardize_structure: &Crystal,
    syminfo: &SymmetryInfo,
    symprec: f64,
) -> Result<(Crystal, Matrix3, Vec<usize>), Box<dyn std::error::Error + Send + Sync>> {
    let (tp, inv_tp) = find_p_matrix(syminfo);

    let lattice_priv = standardize_structure.lattice().change_basis_by(&tp);

    let (positions, species) = (
        standardize_structure.positions(),
        standardize_structure.species(),
    );

    // The matrix are from above defined function safe to cast
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let nvolume = inv_tp.det() as usize;

    let find_position = |sites: &[SiteFraction], p: Vector3<FracCoord>| -> Option<usize> {
        sites
            .iter()
            // this is the Iter::position, return the first index of iter search,
            // a bit confuse under the context of structure's position
            .position(|site| {
                let ps = site.position();
                let p = p.map(|i| f64::from(i) - f64::from(i).floor());
                let ps = ps.map(|i| f64::from(i) - f64::from(i).floor());
                approx_f64(p[0], ps[0], symprec)
                    && approx_f64(p[1], ps[1], symprec)
                    && approx_f64(p[2], ps[2], symprec)
            })
    };

    let mut sites: Vec<SiteFraction> = Vec::with_capacity(positions.len() / nvolume);
    let mut mapping: Vec<usize> = Vec::with_capacity(positions.len());
    // brute forcely filter sites with duplicate position
    for (position, specie) in positions.iter().zip(species.iter()) {
        let new_position = position.change_basis_by(&tp)?;
        if let Some(idx) = find_position(&sites, new_position) {
            // duplicate site
            mapping.push(idx);
        } else {
            let atomic_number = specie.atomic_number();
            mapping.push(sites.len());
            sites.push(SiteFraction::new(new_position, atomic_number));
        }
    }

    #[cfg(debug_assertions)]
    assert_eq!(mapping.len(), positions.len());

    let crystal = CrystalBuilder::new()
        .with_lattice(&lattice_priv)
        .with_frac_sites(sites)
        .build()?;
    Ok((crystal, tp, mapping))
}

#[allow(non_camel_case_types)]
#[derive(Debug)]
pub(crate) enum ExtBravaisClass {
    // Triclinic
    #[allow(dead_code)]
    aP1, // reserved for aP2 + aP3, ref: hpkot paper (Table 94).
    aP2,
    aP3,
    // Monoclinic
    mP1,
    mC1,
    mC2,
    mC3,
    // Orthorhombic
    oP1,
    oA1,
    oA2,
    oC1,
    oC2,
    oF1,
    oF2,
    oF3,
    oI1,
    oI2,
    oI3,
    // Tetragonal
    tP1,
    tI1,
    tI2,
    // Rhombohedral
    hR1,
    hR2,
    // Hexagonal
    hP1,
    hP2,
    // Cubic
    cP1,
    cP2,
    cF1,
    cF2,
    cI1,
}

/// # Errors
/// ??
///
/// # Panics
/// ??
#[allow(clippy::too_many_lines)]
pub fn find_path(
    crystal: &Crystal,
    symprec: f64,
    threshold: f64,
) -> Result<(&'static KpathInfo, KpathEval, Crystal), Box<dyn std::error::Error + Send + Sync>> {
    let syminfo = analyze_symmetry(crystal, symprec)?;
    let structure_std = syminfo.standardize_structure();
    let spg_number = syminfo.spg_number();

    let (structure_priv, _, _) = find_primitive_hpkot(&structure_std, &syminfo, symprec)?;
    let lattice_params = structure_priv.lattice().lattice_params();

    // lattice parameters are for standard conventional cell
    let (a, b, c, alpha, beta, gamma) = structure_std.lattice().lattice_params();
    let a: f64 = a.into();
    let b: f64 = b.into();
    let c: f64 = c.into();
    let _: f64 = alpha.into();
    let beta: f64 = beta.into();
    let _: f64 = gamma.into();

    let ext_bravais = match syminfo.bravais_class() {
        BravaisClass::aP => {
            // get the niggli reduced reciprocal lattice from standard lattice and back to real
            // space.
            // XXX: to get a niggli_reduce this quite cumbersome with type casting... how to
            // improve?? I should find a way that LatticeReciprocal can call niggli_reduce but
            // without the need of ccmat_core depend on moyo. I should make niggli_reduce into a trait.
            let rlatt_niggli_reduced = structure_std.lattice().reciprocal().niggli_reduce()?;

            let (ka, kb, kc, kalpha, kbeta, kgamma) = rlatt_niggli_reduced.lattice_params();

            let ka: f64 = ka.into();
            let kb: f64 = kb.into();
            let kc: f64 = kc.into();
            let kalpha: f64 = kalpha.into();
            let kbeta: f64 = kbeta.into();
            let kgamma: f64 = kgamma.into();

            let mut matrix_mapping: [(f64, TransformationMatrix); 3] = [
                (
                    f64::abs(kb * kc * f64::cos(kalpha)),
                    // a'=b, b'=c, c'=a
                    matrix_3x3![
                        0 0 1;
                        1 0 0;
                        0 1 0;
                    ],
                ),
                (
                    f64::abs(kc * ka * f64::cos(kbeta)),
                    // a'=c, b'=a, c'=b
                    matrix_3x3![
                        0 1 0;
                        0 0 1;
                        1 0 0;
                    ],
                ),
                (
                    f64::abs(ka * kb * f64::cos(kgamma)),
                    // a'=a, b'=b, c'=c
                    matrix_3x3![
                        1 0 0;
                        0 1 0;
                        0 0 1;
                    ],
                ),
            ];

            matrix_mapping.sort_by(|x, y| {
                x.0.partial_cmp(&y.0)
                    .expect("f64::NaN appears in matrix mapping")
            });
            let mt = std::mem::take(&mut matrix_mapping[0].1);
            let lattice = rlatt_niggli_reduced.reciprocal();
            let lattice = lattice.change_basis_by(&mt);
            let klattice = lattice.reciprocal();

            // Make them all-acute or all-obtuse with the additional conditions
            // explained in HPKOT
            // Note: cos > 0 => angle < 90deg

            // TODO: naming?? just kalpha??
            let (_, _, _, kalpha3, kbeta3, kgamma3) = klattice.lattice_params();
            let (kalpha3, kbeta3, kgamma3): (f64, f64, f64) =
                (kalpha3.into(), kbeta3.into(), kgamma3.into());

            if f64::cos(kalpha3).abs() < threshold {
                warn!("aP lattice, but k_alpha3 ~ 90 degrees");
            }

            if f64::cos(kbeta3).abs() < threshold {
                warn!("aP lattice, but k_beta3 ~ 90 degrees");
            }

            if f64::cos(kgamma3).abs() < threshold {
                warn!("aP lattice, but k_gamma3 ~ 90 degrees");
            }

            let s = (
                f64::cos(kalpha3) > 0.0,
                f64::cos(kbeta3) > 0.0,
                f64::cos(kgamma3) > 0.0,
            );

            let m3 = match s {
                // 1a || 1b
                (true, true, true) | (false, false, false) => {
                    matrix_3x3![
                        1, 0, 0;
                        0, 1, 0;
                        0, 0, 1;
                    ]
                }
                // 2a || 2b
                (true, false, false) | (false, true, true) => {
                    matrix_3x3![
                        1,  0, 0;
                        0, -1, 0;
                        0,  0, 1;
                    ]
                }
                // 3a || 3b
                (false, true, false) | (true, false, true) => {
                    matrix_3x3![
                        -1,  0,  0;
                         0,  1,  0;
                         0,  0, -1;
                    ]
                }
                // 4a || 4b
                (false, false, true) | (true, true, false) => {
                    matrix_3x3![
                        -1,  0, 0;
                         0, -1, 0;
                         0,  0, 1;
                    ]
                }
            };

            let klattice = lattice.change_basis_by(&m3).reciprocal();
            let (_, _, _, kalpha3, kbeta3, kgamma3) = klattice.lattice_params();
            let (kalpha3, kbeta3, kgamma3): (f64, f64, f64) =
                (kalpha3.into(), kbeta3.into(), kgamma3.into());

            let s = (
                f64::cos(kalpha3) > 0.0,
                f64::cos(kbeta3) > 0.0,
                f64::cos(kgamma3) > 0.0,
            );

            match s {
                // all-acute
                (true, true, true) => ExtBravaisClass::aP2,
                // all-obtuse
                (false, false, false) => ExtBravaisClass::aP3,
                _ => unreachable!("Unexpected aP triclinic lattice"),
            }
        }
        BravaisClass::mP => ExtBravaisClass::mP1,
        BravaisClass::mC => {
            let cosbeta = f64::cos(beta);
            if f64::abs(b - a * f64::sqrt(1.0 - cosbeta * cosbeta)) < threshold {
                warn!("mC lattice, but b ~ a*sin(beta)");
            }

            if b < a * f64::sqrt(1.0 - cosbeta * cosbeta) {
                ExtBravaisClass::mC1
            } else {
                if f64::abs(-a * cosbeta / c + (a * a) * (1.0 - cosbeta * cosbeta) / (b * b) - 1.0)
                    < threshold
                {
                    warn!("mC lattice, but -a*cos(beta)/c + a^2*sin(beta)^2/b^2 ~ 1");
                }

                if -a * cosbeta / c + (a * a) * (1.0 - cosbeta * cosbeta) / (b * b) < 1.0 {
                    ExtBravaisClass::mC2
                } else {
                    ExtBravaisClass::mC3
                }
            }
        }
        BravaisClass::oP => ExtBravaisClass::oP1,
        BravaisClass::oS => match spg_number {
            // oA
            x if (38..=41).contains(&x) => {
                if f64::abs(b - c) < threshold {
                    warn!("oA lattice, but b ~ c");
                }
                if b < c {
                    ExtBravaisClass::oA1
                } else {
                    ExtBravaisClass::oA2
                }
            }
            // oC
            x if (20..=21).contains(&x) || (35..=37).contains(&x) || (63..=68).contains(&x) => {
                if f64::abs(b - a) < threshold {
                    warn!("oC lattice, but a ~ b");
                }
                if a < b {
                    ExtBravaisClass::oC1
                } else {
                    ExtBravaisClass::oC2
                }
            }
            _ => unreachable!("oS bravais lattice spacegroup number in wrong range"),
        },
        BravaisClass::oF => {
            if f64::abs(1.0 / (a * a) - (1.0 / (b * b) + 1.0 / (c * c))) < threshold {
                warn!("oF lattice, but 1/a^2 ~ 1/b^2 + 1/c^2");
            }
            if f64::abs(1.0 / (c * c) - (1.0 / (a * a) + 1.0 / (b * b))) < threshold {
                warn!("oF lattice, but 1/c^2 ~ 1/a^2 + 1/b^2");
            }
            if 1.0 / (a * a) > 1.0 / (b * b) + 1.0 / (c * c) {
                ExtBravaisClass::oF1
            } else if 1.0 / (c * c) > 1.0 / (a * a) + 1.0 / (b * b) {
                ExtBravaisClass::oF2
            } else {
                ExtBravaisClass::oF3
            }
        }
        BravaisClass::oI => {
            #[derive(Debug)]
            enum Face {
                A, // oI2
                B, // oI3
                C, // oI1
            }
            impl std::fmt::Display for Face {
                fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
                    match self {
                        Face::A => write!(f, "A"),
                        Face::B => write!(f, "B"),
                        Face::C => write!(f, "C"),
                    }
                }
            }

            let mut vec = [(a, Face::A), (b, Face::B), (c, Face::C)];
            vec.sort_by(|x, y| {
                y.0.partial_cmp(&x.0)
                    .expect("lattice length compare impossible to be NaN")
            });

            if f64::abs(vec[0].0 - vec[1].0) < threshold {
                warn!(
                    "oI lattice, but the two longest vectors {} and {} have almost the same length",
                    vec[0].1, vec[1].1,
                );
            }
            match vec[1].1 {
                Face::A => ExtBravaisClass::oI2,
                Face::B => ExtBravaisClass::oI3,
                Face::C => ExtBravaisClass::oI1,
            }
        }
        BravaisClass::tP => ExtBravaisClass::tP1,
        BravaisClass::tI => {
            if (a - c).abs() < threshold {
                warn!("tI lattice, but a ~ c");
            }

            if c < a {
                ExtBravaisClass::tI1
            } else {
                ExtBravaisClass::tI2
            }
        }
        BravaisClass::hR => {
            if f64::abs(f64::sqrt(3.0) * a - f64::sqrt(2.0) * c) < threshold {
                warn!("hR lattice, but sqrt(3)a almost equal to sqrt(2)c");
            }
            if f64::sqrt(3.0) * a < f64::sqrt(2.0) * c {
                ExtBravaisClass::hR1
            } else {
                ExtBravaisClass::hR2
            }
        }
        BravaisClass::hP => {
            // 143..=163 without 150, 152, 154..=156.
            if [
                143, 144, 145, 146, 147, 148, 149, 151, 153, 157, 159, 160, 161, 162, 163,
            ]
            .contains(&spg_number)
            {
                ExtBravaisClass::hP1
            } else {
                ExtBravaisClass::hP2
            }
        }
        BravaisClass::cP => match spg_number {
            x if (195..=206).contains(&x) => ExtBravaisClass::cP1,
            x if (207..=230).contains(&x) => ExtBravaisClass::cP2,
            _ => unreachable!("cP bravais lattice spacegroup number in wrong range"),
        },
        BravaisClass::cF => match spg_number {
            x if (195..=206).contains(&x) => ExtBravaisClass::cF1,
            x if (207..=230).contains(&x) => ExtBravaisClass::cF2,
            _ => unreachable!("cF bravais lattice spacegroup number in wrong range"),
        },
        BravaisClass::cI => ExtBravaisClass::cI1,
    };

    let path_info = path::lookup(&ext_bravais);
    let path_eval = path::eval(path_info, lattice_params)?;

    Ok((path_info, path_eval, structure_priv))
}

#[allow(
    non_snake_case,
    clippy::unreadable_literal,
    clippy::excessive_precision
)]
#[cfg(test)]
mod tests {
    use ccmat_core::{atomic_number, lattice_angstrom, sites_frac_coord, CrystalBuilder};
    use ccmat_symmetry::analyze_symmetry;
    use tracing_test::traced_test;

    use crate::find_path;
    use crate::find_primitive_hpkot;
    use crate::BravaisClass;

    macro_rules! assert_eq_approx_vec3 {
        ($a:expr, $b:expr) => {
            assert_eq_approx_vec3!($a, $b, 1e-12)
        };
        ($a:expr, $b:expr, $tol:expr) => {
            let (a, b) = ($a, $b);
            for i in 0..3 {
                if (a[i] - b[i]).abs() > $tol {
                    panic!(
                        "assertion failed: `{:?} ≈ {:?}`, index {}, diff `{}`, tol `{}`",
                        a,
                        b,
                        i,
                        (a[i] - b[i]).abs(),
                        $tol
                    );
                }
            }
        };
    }

    // this test is only for the purpose to align with the python test in seekpath.
    // The test there is ill-defined, since the function assume the input structure is standardized
    // but the test sample didn't hold that assumption.
    #[test]
    fn find_primitive_hpkot_bcc() {
        let lattice = lattice_angstrom![(4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0),];

        let sites = sites_frac_coord![
            (0.0, 0.0, 0.0), atomic_number!(C);   // C
            (0.5, 0.5, 0.5), atomic_number!(C);   // C
            (0.0, 0.25, 0.0), atomic_number!(O);  // O
            (0.5, 0.75, 0.5), atomic_number!(O);  // O
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-5).unwrap();

        // NOTE: the function assume the s is standardized, but it is not here.
        let (s_priv, _, mapping) = find_primitive_hpkot(&s, &syminfo, 1e-6).unwrap();

        assert_eq!(mapping, [0, 0, 1, 1]);
        assert_eq_approx_vec3!(s_priv.lattice().a().map(f64::from), [-2.0, 2.0, 2.0]);
        assert_eq_approx_vec3!(s_priv.lattice().b().map(f64::from), [2.0, -2.0, 2.0]);
        assert_eq_approx_vec3!(s_priv.lattice().c().map(f64::from), [2.0, 2.0, -2.0]);

        assert_eq_approx_vec3!(s_priv.positions()[0].map(f64::from), [0.0, 0.0, 0.0]);
        assert_eq_approx_vec3!(s_priv.positions()[1].map(f64::from), [0.25, 0.0, 0.25]);
    }

    // same as test above, simply to align with seekpath test
    #[test]
    fn find_primitive_hpkot_oA() {
        let lattice = lattice_angstrom![(9.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0),];

        let sites = sites_frac_coord![
            (0.0, 0.5, 0.46903476), atomic_number!(C);
            (0.0, 0.5, 0.15103982), atomic_number!(O);
            (0.0, 0.0, 0.65103982), atomic_number!(O);
            (0.5, 0.5, 0.87367305), atomic_number!(O);
            (0.0, 0.0, 0.96903476), atomic_number!(C);
            (0.5, 0.0, 0.37367305), atomic_number!(O);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-5).unwrap();

        // NOTE: the function assumes the structure is standardized — same as earlier test
        let (mut s_priv, _, mapping) = find_primitive_hpkot(&s, &syminfo, 1e-6).unwrap();

        // expected primitive cell
        assert_eq_approx_vec3!(s_priv.lattice().a().map(f64::from), [0.0, 1.5, -1.5]);
        assert_eq_approx_vec3!(s_priv.lattice().b().map(f64::from), [0.0, 1.5, 1.5]);
        assert_eq_approx_vec3!(s_priv.lattice().c().map(f64::from), [9.0, 0.0, 0.0]);

        assert_eq!(mapping, [0, 1, 1, 2, 0, 2]);

        let expected_positions = [
            [0.03096524, 0.96903476, 0.0],
            [0.34896018, 0.65103982, 0.0],
            [-0.37367305, 1.37367305, 0.5],
        ];

        for (i, pos) in s_priv.positions().iter().enumerate() {
            assert_eq_approx_vec3!(pos.map(f64::from), expected_positions[i]);
        }

        // for the coordinates wrapping into [0.0, 1.0)
        s_priv.wrap_frac_positions();

        let expected_positions = [
            [0.03096524, 0.96903476, 0.0],
            [0.34896018, 0.65103982, 0.0],
            [0.62632695, 0.37367305, 0.5],
        ];

        for (i, pos) in s_priv.positions().iter().enumerate() {
            assert_eq_approx_vec3!(pos.map(f64::from), expected_positions[i]);
        }
    }

    // different symprec result in getting differnt bravais lattice
    #[test]
    fn symprec_bravais_lattice() {
        let lattice = lattice_angstrom![(4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.00001),];

        let sites = sites_frac_coord![
            (0.0, 0.0, 0.0), atomic_number!(C);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-8).unwrap();
        assert_eq!(syminfo.bravais_class(), BravaisClass::tP);

        let syminfo = analyze_symmetry(&s, 1e-3).unwrap();
        assert_eq!(syminfo.bravais_class(), BravaisClass::cP);
    }

    #[traced_test]
    #[test]
    fn tI_warn() {
        let lattice = lattice_angstrom![(4.0, 0.0, 0.0), (0.0, 4.0, 0.0), (0.0, 0.0, 4.0),];

        let sites = sites_frac_coord![
            (0.0, 0.0, 0.0), atomic_number!(C);
            (0.5, 0.5, 0.5), atomic_number!(C);
            (0.0, 0.0, 0.1), atomic_number!(O);
            (0.5, 0.5, 0.6), atomic_number!(O);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-5).unwrap();
        assert_eq!(syminfo.bravais_class(), BravaisClass::tI);

        let _ = find_path(&s, 1e-5, 1e-7);

        assert!(logs_contain("tI lattice, but a ~ c"));
    }

    #[ignore = "https://github.com/spglib/moyo/issues/191"]
    #[traced_test]
    #[test]
    fn oF1_warn() {
        let lattice = lattice_angstrom![
            (f64::sqrt(1.0 / (1.0 / 16.0 + 1.0 / 25.0)), 0.0, 0.0),
            (0.0, 4.0, 0.0),
            (0.0, 0.0, 5.0),
        ];

        let sites = sites_frac_coord![
            (0.0, 0.0, 0.0), atomic_number!(C);
            (0.0, 0.5, 0.5), atomic_number!(C);
            (0.5, 0.0, 0.5), atomic_number!(C);
            (0.5, 0.5, 0.0), atomic_number!(C);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-5).unwrap();
        assert_eq!(syminfo.bravais_class(), BravaisClass::oF);

        let _ = find_path(&s, 1e-5, 1e-7);

        assert!(logs_contain("oF lattice, but 1/c^2 ~ 1/a^2 + 1/b^2"));
    }

    #[traced_test]
    #[test]
    fn oF2_warn() {
        let lattice = lattice_angstrom![
            (10.0, 0.0, 0.0),
            (0.0, 21.0, 0.0),
            (0.0, 0.0, f64::sqrt(1.0 / (1.0 / 100.0 + 1.0 / 441.0))),
        ];

        let sites = sites_frac_coord![
            (0.1729328200000002, 0.5632488700000001, 0.9531259500000002), atomic_number!(C);
            (0.8270671799999998, 0.4367511299999999, 0.9531259500000002), atomic_number!(C);
            (0.0770671799999998, 0.3132488700000001, 0.7031259500000002), atomic_number!(C);
            (0.9229328200000002, 0.6867511299999999, 0.7031259500000002), atomic_number!(C);
            (0.1729328200000002, 0.0632488700000001, 0.4531259500000002), atomic_number!(C);
            (0.8270671799999998, 0.9367511299999998, 0.4531259500000002), atomic_number!(C);
            (0.0770671799999998, 0.8132488700000001, 0.2031259500000002), atomic_number!(C);
            (0.9229328200000002, 0.1867511299999999, 0.2031259500000002), atomic_number!(C);
            (0.6729328200000002, 0.5632488700000001, 0.4531259500000002), atomic_number!(C);
            (0.3270671799999998, 0.4367511299999999, 0.4531259500000002), atomic_number!(C);
            (0.5770671799999998, 0.3132488700000001, 0.2031259500000002), atomic_number!(C);
            (0.4229328200000002, 0.6867511299999999, 0.2031259500000002), atomic_number!(C);
            (0.6729328200000002, 0.0632488700000001, 0.9531259500000002), atomic_number!(C);
            (0.3270671799999998, 0.9367511299999998, 0.9531259500000002), atomic_number!(C);
            (0.5770671799999998, 0.8132488700000001, 0.7031259500000002), atomic_number!(C);
            (0.4229328200000002, 0.1867511299999999, 0.7031259500000002), atomic_number!(C);
            (0.0, 0.5, 0.4701481000000003), atomic_number!(O);
            (0.75, 0.75, 0.2201481000000003), atomic_number!(O);
            (0.0, 0.0, 0.9701481000000002), atomic_number!(O);
            (0.75, 0.25, 0.7201481000000003), atomic_number!(O);
            (0.5, 0.5, 0.9701481000000002), atomic_number!(O);
            (0.25, 0.75, 0.7201481000000003), atomic_number!(O);
            (0.5, 0.0, 0.4701481000000003), atomic_number!(O);
            (0.25, 0.25, 0.2201481000000003), atomic_number!(O);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-5).unwrap();
        assert_eq!(syminfo.bravais_class(), BravaisClass::oF);

        let _ = find_path(&s, 1e-5, 1e-7);

        assert!(logs_contain("oF lattice, but 1/c^2 ~ 1/a^2 + 1/b^2"));
    }

    #[traced_test]
    #[test]
    fn oI_bc_warn() {
        let lattice = lattice_angstrom![(4.0, 0.0, 0.0), (0.0, 5.0, 0.0), (0.0, 0.0, 5.0),];

        let sites = sites_frac_coord![
            (0.0, 0.0, 0.0), atomic_number!(C);
            (0.5, 0.5, 0.5), atomic_number!(C);
            (0.0, 0.0, 0.1), atomic_number!(O);
            (0.5, 0.5, 0.6), atomic_number!(O);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-5).unwrap();
        assert_eq!(syminfo.bravais_class(), BravaisClass::oI);

        let _ = find_path(&s, 1e-5, 1e-7);

        assert!(logs_contain(
            "oI lattice, but the two longest vectors B and C have almost the same length"
        ));
    }

    #[traced_test]
    #[test]
    fn oC_warn() {
        let lattice = lattice_angstrom![(3.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 5.0),];

        let sites = sites_frac_coord![
            (0.5000000000000000, 0.1136209299999999, 0.7500967299999999), atomic_number!(C);
            (0.5000000000000000, 0.8863790700000000, 0.2500967299999999), atomic_number!(C);
            (0.0000000000000000, 0.6136209300000000, 0.7500967299999999), atomic_number!(C);
            (0.0000000000000000, 0.3863790700000001, 0.2500967299999999), atomic_number!(C);
            (0.0000000000000000, 0.8444605049999999, 0.7659032699999999), atomic_number!(O);
            (0.0000000000000000, 0.1555394950000001, 0.2659032699999999), atomic_number!(O);
            (0.5000000000000000, 0.3444605049999999, 0.7659032699999999), atomic_number!(O);
            (0.5000000000000000, 0.6555394950000001, 0.2659032699999999), atomic_number!(O);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-5).unwrap();
        assert_eq!(syminfo.bravais_class(), BravaisClass::oS);
        assert_eq!(syminfo.spg_number(), 36);

        let _ = find_path(&s, 1e-5, 1e-7);

        assert!(logs_contain("oC lattice, but a ~ b"));
    }

    #[traced_test]
    #[test]
    fn oA_warn() {
        let lattice = lattice_angstrom![(9.0, 0.0, 0.0), (0.0, 3.0, 0.0), (0.0, 0.0, 3.0),];

        let sites = sites_frac_coord![
            (0.0000000000000000, 0.0000000000000000, 0.0309652399999998), atomic_number!(C);
            (0.0000000000000000, 0.5000000000000000, 0.5309652399999998), atomic_number!(C);
            (0.0000000000000000, 0.5000000000000000, 0.8489601849999999), atomic_number!(O);
            (0.5000000000000000, 0.5000000000000000, 0.1263269549999999), atomic_number!(O);
            (0.0000000000000000, 0.0000000000000000, 0.3489601849999999), atomic_number!(O);
            (0.5000000000000000, 0.0000000000000000, 0.6263269549999999), atomic_number!(O);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&s, 1e-5).unwrap();
        assert_eq!(syminfo.bravais_class(), BravaisClass::oS);
        assert_eq!(syminfo.spg_number(), 38);

        let _ = find_path(&s, 1e-5, 1e-7);

        assert!(logs_contain("oA lattice, but b ~ c"));
    }

    // TODO: test on mC and oP for warning messages as well as planed in seekpath
}
