/*
extxyz next generation specification:

The specification backward compatible with the extxyz specification but become with more formalized format in writing specification to make it non-ambiguse in representing the structure.

It is call "xyz", so always 3D by default, if it is used for 2d, the z is all 0.

The positions are always in cartesian coordinates in angstrom, since original extxyz lack of ability to pre-parse.
Meanwhile, it was originally for molecule, it is straightforward to assume the cartesian coordinates.

in the ng, there are default properties if not provided, the default properties try to backward compatible with libAtom/extxyz.
If assigned in the key, those settings can be override.
- lattice unit (angstrom by default), only if Lattice in the "Properties"
- pbc (true:true:true) and the shape of Properties is pos:3
- species:1:pos:3 -> if there is no such settings in the "Properties" this is the default setting. "libAtom/extxyz" use this already.
- the legacy libAtom/extxyz is quit unusable, the regular xyz file will result in fail parsing.
- the parser are suppose to aggresively parse any valid pair of key=value and leave the rest "comment". if no key=value pairs the whole line is a comment. otherwise raise errors.
*/
use std::io::BufRead;

use ccmat_core::{
    atomic_number_from_symbol, math::Vector3, Angstrom, Crystal, CrystalBuilder, Lattice, Molecule,
    MoleculeBuilder, SiteCartesian,
};

#[derive(Debug)]
pub enum ParseError {
    WrongParser,
    ParseFailed {
        source: Box<dyn std::error::Error + Send + Sync>,
    },
    ConvertFailed {
        message: String,
    },
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::WrongParser => write!(f, "wrong parser"),
            ParseError::ParseFailed { source } => write!(f, "parsing failed from - {source}"),
            ParseError::ConvertFailed { message } => write!(
                f,
                "unable to covert and construct to the structure after parsing: {message}"
            ),
        }
    }
}

#[derive(Debug)]
pub enum Structure {
    Crystal(Crystal),
    Molecule(Molecule),
}

impl TryFrom<extxyz::Frame> for Structure {
    type Error = ParseError;

    fn try_from(value: extxyz::Frame) -> Result<Self, Self::Error> {
        let natoms = value.natoms();
        let info_map = value.info();
        let arrs_map = value.arrs();

        let species: Vec<u8> = if let Some(v) = arrs_map.get("species") {
            if let extxyz::Value::VecText(v, n) = v {
                debug_assert_eq!(*n, natoms);

                v.iter()
                    .map(|t| atomic_number_from_symbol(t))
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|err| ParseError::ConvertFailed {
                        message: format!("{err}"),
                    })?
            } else {
                Err(ParseError::ConvertFailed {
                    message: "species must be a vec of Text".to_string(),
                })?
            }
        } else {
            Err(ParseError::ConvertFailed {
                message: "'species' not parsed in the properties, wrong input format".to_string(),
            })?
        };

        let positions: Vec<Vector3<Angstrom>> = if let Some(m) = arrs_map.get("pos") {
            if let extxyz::Value::MatrixFloat(m, shape) = m {
                debug_assert_eq!(shape.0, natoms);

                // only try to support 3d structure at the moment, make sense the crate is ext"xyz"
                debug_assert_eq!(shape.1, 3);

                m.iter()
                    .map(|pos| {
                        let mut vp = [Angstrom::from(0.0); 3];
                        for (i, &p) in pos.iter().enumerate() {
                            vp[i] = Angstrom::from(*p);
                        }
                        Vector3(vp)
                    })
                    .collect()
            } else {
                Err(ParseError::ConvertFailed {
                    message: "species must be a matrix of float".to_string(),
                })?
            }
        } else {
            Err(ParseError::ConvertFailed {
                message: "'pos' not parsed in the properties, wrong input format".to_string(),
            })?
        };

        let sites_cart: Vec<SiteCartesian> = positions
            .into_iter()
            .zip(species)
            .map(|(pos, atomic_number)| SiteCartesian::new(pos, atomic_number))
            .collect();

        if let Some(v) = info_map.get("Lattice") {
            // a crystal
            if let extxyz::Value::MatrixFloat(v, (nrows, ncols)) = v {
                if *nrows != 3 || *ncols != 3 {
                    return Err(ParseError::ConvertFailed {
                        message: "Lattice must has shape (3, 3)".to_string(),
                    });
                }
                let mut latt = [[0.0; 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        latt[i][j] = *v[i][j];
                    }
                }
                let latt = Lattice::from_angstroms(latt);
                let crystal = CrystalBuilder::new()
                    .with_lattice(&latt)
                    .with_cart_sites(sites_cart)
                    .build()
                    .map_err(|err| ParseError::ConvertFailed {
                        message: format!("cannot build the crystal, {err}"),
                    })?;

                Ok(Structure::Crystal(crystal))
            } else {
                Err(ParseError::ConvertFailed {
                    message: "Lattice must be a matrix of floats".to_string(),
                })
            }
        } else {
            // molecule
            let mol = MoleculeBuilder::new()
                .with_sites(sites_cart)
                .build()
                .map_err(|err| ParseError::ConvertFailed {
                    message: format!("cannot build the molecule, {err}"),
                })?;

            Ok(Structure::Molecule(mol))
        }
    }
}

impl std::error::Error for ParseError {}

///
/// # Errors
/// ???
pub fn parse<R>(r: &mut R, ext: &str) -> Result<Structure, ParseError>
where
    R: BufRead,
{
    match ext {
        "xyz" => {
            let s: Structure = extxyz::read_frame(r)
                .map_err(|err| ParseError::ParseFailed {
                    source: Box::new(err),
                })?
                .try_into()?;
            Ok(s)
        }
        "cif" => {
            todo!()
        }
        _ => unimplemented!(),
    }
}

#[allow(clippy::similar_names)]
#[cfg(test)]
mod tests {
    use std::io::Cursor;

    use ccmat_core::{lattice_angstrom, FracCoord};

    use super::*;

    macro_rules! assert_eq_approx {
        ($a:expr, $b:expr) => {{
            assert_eq_approx!($a, $b, 1e-12)
        }};
        ($a:expr, $b:expr, $tol:expr) => {{
            let (left, right) = ($a, $b);
            if (left - right).abs() > $tol {
                panic!(
                    "assertion failed: `{} ≈ {}`, diff:  `{}`, tol: `{}`",
                    left,
                    right,
                    (left - right).abs(),
                    $tol
                );
            }
        }};
    }

    #[test]
    fn test_parse_crystal_default() {
        let ext = "xyz";
        let content = include_str!("../assets/sio2.xyz");

        let mut rd = Cursor::new(content.as_bytes());
        let Structure::Crystal(c) = parse(&mut rd, ext).unwrap() else {
            panic!()
        };

        let latt_got = c.lattice();
        let latt_expect = lattice_angstrom![
            a = (5.0, 0.0, 0.0),
            b = (0.0, 5.0, 0.0),
            c = (0.0, 0.0, 5.0),
        ];

        for i in 0..2 {
            assert_eq_approx!(f64::from(latt_got.a()[i]), f64::from(latt_expect.a()[i]));
            assert_eq_approx!(f64::from(latt_got.b()[i]), f64::from(latt_expect.b()[i]));
            assert_eq_approx!(f64::from(latt_got.c()[i]), f64::from(latt_expect.c()[i]));
        }

        // cartesian convert to frac
        let pos_got_oxygen = c.positions()[2];
        let pos_expect = Vector3([
            FracCoord::from(0.25),
            FracCoord::from(0.25),
            FracCoord::from(0.25),
        ]);

        for i in 0..3 {
            assert_eq_approx!(f64::from(pos_got_oxygen[i]), f64::from(pos_expect[i]));
        }

        assert_eq!(c.species()[2].atomic_number(), 8);
    }

    #[test]
    fn test_parse_mol_default() {
        let ext = "xyz";
        let content = include_str!("../assets/benz.xyz");

        let mut rd = Cursor::new(content.as_bytes());
        let Structure::Molecule(c) = parse(&mut rd, ext).unwrap() else {
            panic!()
        };

        // cartesian convert to frac
        let pos_got_c5 = c.positions()[5];
        let pos_expect_c5 = Vector3([
            Angstrom::from(-1.209_657),
            Angstrom::from(0.698_396),
            Angstrom::from(0.0),
        ]);

        let pos_got_h5 = c.positions()[11];
        let pos_expect_h5 = Vector3([
            Angstrom::from(-2.156_659),
            Angstrom::from(1.245_145),
            Angstrom::from(0.0),
        ]);

        for i in 0..3 {
            assert_eq_approx!(f64::from(pos_got_c5[i]), f64::from(pos_expect_c5[i]));
            assert_eq_approx!(f64::from(pos_got_h5[i]), f64::from(pos_expect_h5[i]));
        }

        assert_eq!(c.species()[2].atomic_number(), 6);
        assert_eq!(c.species()[10].atomic_number(), 1);
    }
}
