/* Interop with moyo for symmetry utils.
 *
 * Here I interact with crate::moyo, and keep actual moyo behind it, use it from user's view
 * on how I think the moyo APIs should be better reshaped.
 * The goal is to have `moyo_wrapper` drop in replaced by `moyo`.
 *
 */

// TODO: thiserror

use std::borrow::Cow;

use ccmat_core::math::{Matrix3, Vector3};
use ccmat_core::{
    lattice_angstrom, Basis, Crystal, CrystalBuilder, FracCoord, HasBasis, SiteFraction,
};
use std::cmp;

use crate::moyo_wrapper;

/// delegation of `moyo_wrapper` to ccmat API users.
pub struct SymmetryInfo {
    inner: moyo_wrapper::SymmetryInfo,
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Centering {
    P, // Primitive
    A, // A-face centered
    B, // B-face centered
    C, // C-face centered
    I, // Body centered
    R, // Rhombohedral (obverse setting)
    F, // Face centered
}

#[allow(non_camel_case_types)]
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BravaisClass {
    // Triclinic
    aP,
    // Monoclinic
    mP,
    mC,
    // Orthorhombic
    oP,
    oS,
    oF,
    oI,
    // Tetragonal
    tP,
    tI,
    // Rhombohedral
    hR,
    // Hexagonal
    hP,
    // Cubic
    cP,
    cF,
    cI,
}

impl std::fmt::Display for BravaisClass {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            // Triclinic
            BravaisClass::aP => "aP",
            // Monoclinic
            BravaisClass::mP => "mP",
            BravaisClass::mC => "mC",
            // Orthorhombic
            BravaisClass::oP => "oP",
            BravaisClass::oS => "oS",
            BravaisClass::oF => "oF",
            BravaisClass::oI => "oI",
            // Tetragonal
            BravaisClass::tP => "tP",
            BravaisClass::tI => "tI",
            // Rhombohedral
            BravaisClass::hR => "hR",
            // Hexagonal
            BravaisClass::hP => "hP",
            // Cubic
            BravaisClass::cP => "cP",
            BravaisClass::cF => "cF",
            BravaisClass::cI => "cI",
        };
        write!(f, "{s}")
    }
}

impl From<moyo_wrapper::BravaisClass> for BravaisClass {
    fn from(bv: moyo_wrapper::BravaisClass) -> Self {
        match bv.raw() {
            moyo::data::BravaisClass::aP => BravaisClass::aP,
            moyo::data::BravaisClass::mP => BravaisClass::mP,
            moyo::data::BravaisClass::mC => BravaisClass::mC,
            moyo::data::BravaisClass::oP => BravaisClass::oP,
            moyo::data::BravaisClass::oS => BravaisClass::oS,
            moyo::data::BravaisClass::oF => BravaisClass::oF,
            moyo::data::BravaisClass::oI => BravaisClass::oI,
            moyo::data::BravaisClass::tP => BravaisClass::tP,
            moyo::data::BravaisClass::tI => BravaisClass::tI,
            moyo::data::BravaisClass::hR => BravaisClass::hR,
            moyo::data::BravaisClass::hP => BravaisClass::hP,
            moyo::data::BravaisClass::cP => BravaisClass::cP,
            moyo::data::BravaisClass::cF => BravaisClass::cF,
            moyo::data::BravaisClass::cI => BravaisClass::cI,
        }
    }
}

impl From<moyo_wrapper::Centering> for Centering {
    fn from(c: moyo_wrapper::Centering) -> Self {
        match c.raw() {
            moyo::data::Centering::P => Centering::P,
            moyo::data::Centering::A => Centering::A,
            moyo::data::Centering::B => Centering::B,
            moyo::data::Centering::C => Centering::C,
            moyo::data::Centering::I => Centering::I,
            moyo::data::Centering::R => Centering::R,
            moyo::data::Centering::F => Centering::F,
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
        let cell = self.inner.std_cell();
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
        let sites: Vec<SiteFraction> = positions
            .iter()
            .zip(numbers)
            .map(|(pos, num)| {
                let pos = pos.map(FracCoord::from);
                let pos = Vector3(pos);
                let num: u8 = (*num).try_into().expect("atomic number not in 0..128");
                SiteFraction::new(pos, num)
            })
            .collect();

        // build_uncheck since moyo::Cell is the intermediate structure from Crystal and guranteed
        // to be valid, otherwise it is a bug.
        CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build_uncheck()
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
    s: &Crystal,
    symprec: f64,
) -> Result<SymmetryInfo, Box<dyn std::error::Error + Send + Sync>> {
    let a = s.lattice().a().map(f64::from);
    let b = s.lattice().b().map(f64::from);
    let c = s.lattice().c().map(f64::from);
    let lattice = moyo_wrapper::__macro::lattice!(a, b, c);

    // TODO: moyo need an api or macro to create positions.
    let positions = s
        .positions_fraction()
        .iter()
        .map(|p| [f64::from(p[0]), f64::from(p[1]), f64::from(p[2])])
        .collect();

    let numbers = s
        .species()
        .iter()
        .map(|s| s.atomic_number().into())
        .collect();

    let cell = moyo_wrapper::CellBuilder::new()
        .with_lattice(lattice)
        .with_positions(positions)
        .with_numbers(numbers)
        .build();
    let inner = moyo_wrapper::analyze_symmetry(&cell, symprec)?;
    let sym_info = SymmetryInfo { inner };
    Ok(sym_info)
}

// This can become a symmetry extent trait
// TODO: this trait can be moved to moyo, as the proposed way of better interface.
pub trait SymmetryExt {
    /// 3x3 lattice in angstrom
    fn lattice(&self) -> [[f64; 3]; 3];
    /// positions in fraction
    fn positions(&self) -> Vec<[f64; 3]>;
    /// atomic number
    fn numbers(&self) -> Vec<i32>;

    /// # Errors
    /// TODO: ???
    fn is_supercell(&self, symprec: f64) -> Result<bool, Box<dyn std::error::Error + Send + Sync>> {
        let m = self.lattice();
        let (a, b, c) = (m[0], m[1], m[2]);
        let lattice = moyo_wrapper::Lattice::new(a, b, c);
        let positions = self.positions();
        let numbers = self.numbers();

        let natoms = numbers.len();
        let cell = moyo_wrapper::CellBuilder::new()
            .with_lattice(lattice)
            .with_positions(positions)
            .with_numbers(numbers)
            .build();
        let sym_info = moyo_wrapper::analyze_symmetry(&cell, symprec)?;
        let cell_priv = sym_info.prim_std_cell();
        match natoms.cmp(&cell_priv.numbers().len()) {
            cmp::Ordering::Less => unreachable!("primitive cell cannot have more atoms."),
            cmp::Ordering::Equal => Ok(false),
            cmp::Ordering::Greater => Ok(true),
        }
    }

    // TODO: this can become the generic interface of moyo symmetry analysis call
}

impl SymmetryExt for Crystal {
    fn lattice(&self) -> [[f64; 3]; 3] {
        let a: [f64; 3] = self.lattice().a().map(f64::from);
        let b: [f64; 3] = self.lattice().b().map(f64::from);
        let c: [f64; 3] = self.lattice().c().map(f64::from);
        [a, b, c]
    }

    fn positions(&self) -> Vec<[f64; 3]> {
        self.positions_fraction()
            .iter()
            .map(|p| p.map(f64::from))
            .collect()
    }

    fn numbers(&self) -> Vec<i32> {
        self.species()
            .iter()
            .map(|s| i32::from(s.atomic_number()))
            .collect()
    }
}

pub trait NiggliReduceExt: Sized {
    fn niggli_reduce(&self) -> Result<Self, Box<dyn std::error::Error + Send + Sync>>;
}

impl<T> NiggliReduceExt for T
where
    T: HasBasis + From<Basis>,
{
    fn niggli_reduce(&self) -> Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        let basis = self.basis();
        let basis = basis.map(|v| *v);
        let (basis, _) = moyo_wrapper::niggli_reduce(basis)?;
        let basis: [Vector3<f64>; 3] = basis.map(Vector3);
        Ok(Self::from(basis))
    }
}

#[cfg(test)]
mod tests {
    use ccmat_core::{atomic_number, lattice_angstrom, sites_frac_coord, CrystalBuilder};

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
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let syminfo = analyze_symmetry(&crystal, 1e-4).unwrap();
        assert_eq!(syminfo.spg_number(), 136);
    }
}
