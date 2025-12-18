use std::io::BufRead;

use ccmat_core::{math::Vector3, Angstrom, Crystal, CrystalBuilder, Lattice, Molecule};

#[derive(Debug)]
enum ParseError {
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

enum Structure {
    Crystal(Crystal),
    Molecule(Molecule),
}

impl TryFrom<extxyz::Frame> for Structure {
    type Error = ParseError;

    fn try_from(value: extxyz::Frame) -> Result<Self, Self::Error> {
        let natoms = value.natoms();
        let info = value.info();
        if let Some(v) = info.get("Lattice") {
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
                let a: Vector3<Angstrom> = Vector3(latt[0].map(|v| v.into()));
                let b: Vector3<Angstrom> = Vector3(latt[1].map(|v| v.into()));
                let c: Vector3<Angstrom> = Vector3(latt[2].map(|v| v.into()));
                let latt = Lattice::new(a, b, c);
                // CrystalBuilder::new()
                //     .with_lattice(&latt)
                //     .with_sites(sites)
                //     .build()
                //     .map_err(|err| ParseError::ConvertFailed {
                //         message: format!("cannot build the crystal, {err}"),
                //     })?;
                todo!()
            } else {
                Err(ParseError::ConvertFailed {
                    message: "Lattice must be a matrix of floats".to_string(),
                })
            }
        } else {
            // molecule
            todo!()
        }
    }
}

impl std::error::Error for ParseError {}

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
