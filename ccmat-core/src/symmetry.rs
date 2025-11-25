/* Interop with moyo for symmetry utils.
 *
 * Here I interact with crate::moyo, and keep actual moyo behind it, use it from user's view
 * on how I think the moyo APIs should be better reshaped.
 * The goal is to have `moyo_wrapper` drop in replaced by `moyo`.
 *
 */

// TODO: thiserror

use std::borrow::Cow;

use crate::math::{Matrix3, TransformationMatrix, Vector3};
use crate::moyo_wrapper;
use crate::structure::Centering;
use crate::BravaisClass;
use crate::{lattice_angstrom, Crystal, CrystalBuilder, FracCoord, Site};

/// delegation of `moyo_wrapper` to ccmat API users.
pub struct SymmetryInfo {
    inner: moyo_wrapper::SymmetryInfo,
}

impl From<moyo_wrapper::BravaisClass> for BravaisClass {
    fn from(bv: moyo_wrapper::BravaisClass) -> Self {
        match bv {
            moyo_wrapper::BravaisClass::aP => BravaisClass::aP,
            moyo_wrapper::BravaisClass::mP => BravaisClass::mP,
            moyo_wrapper::BravaisClass::mC => BravaisClass::mC,
            moyo_wrapper::BravaisClass::oP => BravaisClass::oP,
            moyo_wrapper::BravaisClass::oS => BravaisClass::oS,
            moyo_wrapper::BravaisClass::oF => BravaisClass::oF,
            moyo_wrapper::BravaisClass::oI => BravaisClass::oI,
            moyo_wrapper::BravaisClass::tP => BravaisClass::tP,
            moyo_wrapper::BravaisClass::tI => BravaisClass::tI,
            moyo_wrapper::BravaisClass::hR => BravaisClass::hR,
            moyo_wrapper::BravaisClass::hP => BravaisClass::hP,
            moyo_wrapper::BravaisClass::cP => BravaisClass::cP,
            moyo_wrapper::BravaisClass::cF => BravaisClass::cF,
            moyo_wrapper::BravaisClass::cI => BravaisClass::cI,
        }
    }
}

impl From<moyo_wrapper::Centering> for Centering {
    fn from(c: moyo_wrapper::Centering) -> Self {
        match c {
            moyo_wrapper::Centering::P => Centering::P,
            moyo_wrapper::Centering::A => Centering::A,
            moyo_wrapper::Centering::B => Centering::B,
            moyo_wrapper::Centering::C => Centering::C,
            moyo_wrapper::Centering::I => Centering::I,
            moyo_wrapper::Centering::R => Centering::R,
            moyo_wrapper::Centering::F => Centering::F,
        }
    }
}

impl SymmetryInfo {
    /// Space group number (1-230)
    #[must_use]
    pub fn spg_number(&self) -> u32 {
        self.inner.spg_number()
    }

    /// Hall symbol number (1-530).
    #[must_use]
    pub fn hall_number(&self) -> u32 {
        self.inner.hall_number()
    }

    /// Bravais class
    #[must_use]
    pub fn bravais_class(&self) -> BravaisClass {
        self.inner.bravais_class().into()
    }

    /// Hall symbol
    #[must_use]
    pub fn hall_symbol(&self) -> Cow<'_, str> {
        self.inner.hall_symbol()
    }

    /// Spage group symbol (aka "Hermann–Mauguin (International) symbol")
    #[must_use]
    pub fn spacegroup_symbol(&self) -> Cow<'_, str> {
        self.inner.spagegroup_symbol()
    }

    /// Check if contain inversion symmetry.
    #[must_use]
    pub fn has_inversion(&self) -> bool {
        self.inner.has_inversion()
    }

    /// Crystal in standard structure
    #[must_use]
    pub fn standardize_structure(&self) -> Crystal {
        self.inner.std_cell().into()
    }

    pub fn std_rotation(&self) -> Matrix3 {
        Matrix3(self.inner.std_rotation())
    }
}

/// analyze symmetry
///
/// # Errors
/// moyo not able to analyze the symmetry of the structure.
pub fn analyze_symmetry(
    crystal: &Crystal,
    symprec: f64,
) -> Result<SymmetryInfo, Box<dyn std::error::Error + Send + Sync>> {
    let cell = crystal.into();
    let inner = moyo_wrapper::analyze_symmetry(&cell, symprec)?;
    let sym_info = SymmetryInfo { inner };
    Ok(sym_info)
}

type Basis = [Vector3<f64>; 3];

/// further wrap ``moyo::math::niggili::niggli_reduce`` (however not exposed) into function where the types ccmat confortable to work with.
pub(crate) fn niggli_reduce(
    basis: Basis,
) -> Result<(Basis, TransformationMatrix), Box<dyn std::error::Error + Send + Sync>> {
    let basis = basis.map(|v| *v);
    match moyo_wrapper::niggli_reduce(basis) {
        Ok(result) => {
            let basis: [Vector3<f64>; 3] = result.0.map(Vector3);
            let mt = result
                .1
                .map(|v| [f64::from(v[0]), f64::from(v[1]), f64::from(v[2])]);
            Ok((basis, Matrix3(mt)))
        }
        Err(err) => Err(format!("niggli reduction failed (moyo as symmery engine): {err}").into()),
    }
}

impl From<&Crystal> for moyo_wrapper::Cell {
    fn from(s: &Crystal) -> Self {
        let a = s.lattice().a().map(f64::from);
        let b = s.lattice().b().map(f64::from);
        let c = s.lattice().c().map(f64::from);
        let lattice = moyo_wrapper::__macro::lattice!(a, b, c);

        // TODO: moyo need an api or macro to create positions.
        let positions = s
            .positions()
            .iter()
            .map(|p| [f64::from(p[0]), f64::from(p[1]), f64::from(p[2])])
            .collect();

        let numbers = s
            .species()
            .iter()
            .map(|s| s.atomic_number().into())
            .collect();

        moyo_wrapper::CellBuilder::new()
            .with_lattice(lattice)
            .with_positions(positions)
            .with_numbers(numbers)
            .build()
    }
}

// If willing to give the ownership.
impl From<Crystal> for moyo_wrapper::Cell {
    fn from(s: Crystal) -> Self {
        (&s).into()
    }
}

// The reference version not provide, since the moyo::Cell only used internally thus
// assume no need to hold its ownership.
impl From<moyo_wrapper::Cell> for Crystal {
    fn from(cell: moyo_wrapper::Cell) -> Self {
        let a = cell.lattice().basis().0;
        let b = cell.lattice().basis().1;
        let c = cell.lattice().basis().2;

        let lattice = lattice_angstrom![
            a = (a[0], a[1], a[2]),
            b = (b[0], b[1], b[2]),
            c = (c[0], c[1], c[2]),
        ];

        let positions = cell.positions();
        let numbers = cell.numbers();

        // Length are guranteed to be the same because the moyo::Cell is constructed from ccmat
        // use `into_iter` to move and avoid allocation.
        let sites: Vec<Site> = positions
            .iter()
            .zip(numbers)
            .map(|(pos, num)| {
                let pos = pos.map(FracCoord::from);
                let pos = Vector3(pos);
                let num: u8 = (*num).try_into().expect("atomic number not in 0..128");
                Site::new(pos, num)
            })
            .collect();

        // build_uncheck since moyo::Cell is the intermediate structure from Crystal and guranteed
        // to be valid, otherwise it is a bug.
        CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_sites(&sites)
            .build_uncheck()
    }
}

#[cfg(test)]
mod tests {
    use crate::{atomic_number, lattice_angstrom, sites_frac_coord, CrystalBuilder};

    use super::*;

    #[test]
    fn spacegroup() {
        const X_4F: f64 = 0.3046;

        let lattice = lattice_angstrom![
            a = (4.603, 0.0, 0.0),
            b = (0.0, 4.603, 0.0),
            c = (0.0, 0.0, 4.603),
        ];
        let sites = sites_frac_coord![
            (0.0, 0.0, 0.0), atomic_number!(Ti);               // Ti(2a)
            (0.5, 0.5, 0.5), atomic_number!(Ti);               // Ti(2a)
            (X_4F, X_4F, 0.0), atomic_number!(O);              // O(4f)
            (1.0 - X_4F, 1.0 - X_4F, 0.0), atomic_number!(O);  // O(4f)
            (-X_4F + 0.5, X_4F + 0.5, 0.5), atomic_number!(O); // O(4f)
            (X_4F + 0.5, -X_4F + 0.5, 0.5), atomic_number!(O); // O(4f)
        ];
        let crystal = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_sites(&sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&crystal, 1e-4).unwrap();
        assert_eq!(syminfo.spg_number(), 136);
    }
}
