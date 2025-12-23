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
        _ => unimplemented!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_crystal_default() {}

    #[test]
    fn test_parse_mol_default() {}
}
