/* Experiments (wishlist) on moyo APIs
 *
 * Everything in this moyo_wrapper mod are expected to go to official moyo crate.
 * In experiment, this mod should only depend on official moyo crate.
 *
 * I suggest to make most of fields private and have field access APIs. This can benifit for
 * keeping the current moyo inner implementation untouched, and leave the door for future
 * improvment without dramaticall breaking the API.
 *
 * The API design I had is based on my perspective on the usage of moyo is for getting all symmetry
 * related information by one pass. Therefore I make the assumption that the inner structure is
 * immutable in the process. If there are new structure expected in the output such as primitive or
 * standardized cell, moyo make a new allocation for those instead of mutate the original input
 * structure.
 * Based on this assumption, the function to return the part of the cell (positions, numbers,
 * lattice) are all return the view (the reference) of the original data.
 *
 */
use std::borrow::Cow;

// NOTE: I expecte
// - moyo::base::MoyoError -> moyo::Error,
// - data::{arithmetic_crystal_class_entry, hall_symbol_entry} are private methods.
// - MoyoDataSet -> SymmetryInfo as I show and impl below
use moyo::{
    self,
    data::{arithmetic_crystal_class_entry, hall_symbol_entry},
    MoyoDataset,
};

pub(crate) struct BravaisClass {
    pub(crate) inner: moyo::data::BravaisClass,
}

pub(crate) struct Centering {
    pub(crate) inner: moyo::data::Centering,
}

pub(crate) mod __macro {
    macro_rules! __vec3_angstrom {
        ([$x:expr, $y:expr, $z:expr]) => {
            [$x, $y, $z]
        };
        (($x:expr, $y:expr, $z:expr)) => {
            [$x, $y, $z]
        };
    }

    macro_rules! lattice {
        (
            $a:tt,
            $b:tt,
            $c:tt $(,)?
        ) => {{
            $crate::moyo_wrapper::Lattice::new($a, $b, $c)
        }};
    }
    pub(crate) use lattice;

    /// macro to create a positions list
    ///
    /// # Examples
    ///
    /// ```ignore
    /// use ccmat_core::sites_frac_coord;
    ///
    /// let _ = positions![
    ///     0.0, 0.0, 0.0;
    ///     0.0, 0.0, 0.5;
    /// ];
    /// ```
    #[allow(unused_macros)]
    macro_rules! positions {
        () => {
            // empty cell
            Vec::new()
        };
        ( $(
            $x:expr,$y:expr,$z:expr,
          );+ $(;)?
        ) => {{
            let positions = vec![
                $(
                    position!($x, $y, $z)
                ),+
            ];
            positions
        }};
    }
}

pub(crate) struct LatticeSet;
pub(crate) struct LatticeNotSet;
pub(crate) struct PositionsSet;
pub(crate) struct PositionsNotSet;
pub(crate) struct NumbersSet;
pub(crate) struct NumbersNotSet;

pub(crate) struct Lattice {
    inner: moyo::base::Lattice,
}

impl Lattice {
    #[must_use]
    pub fn new(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> Self {
        Lattice {
            inner: moyo::base::Lattice::from_basis([a, b, c]),
        }
    }

    /// the final API to access basis I expect is:
    ///
    /// let a = cell.lattice().basis().0;
    /// let b = cell.lattice().basis().1;
    /// let c = cell.lattice().basis().2;
    ///
    /// I want to avoid any specific knownledge on the return type such as linear algebra data
    /// structure from `nalgebra`.
    pub(crate) fn basis(&self) -> ([f64; 3], [f64; 3], [f64; 3]) {
        let a = self.inner.basis.column(0);
        let b = self.inner.basis.column(1);
        let c = self.inner.basis.column(2);

        ([a[0], a[1], a[2]], [b[0], b[1], b[2]], [c[0], c[1], c[2]])
    }
}

/// Wrapper of `moyo::base::Cell` to explore idiomatic API design for moyo.
pub(crate) struct Cell {
    inner: moyo::base::Cell,
}

impl Cell {
    /// Ideally the return type can be a view a reference of cell lattice.
    ///
    /// However, I have to create the wrapper type, otherwise it is tricky to do memory mapping from
    /// inner lattice to the wrapper lattice.
    /// After integrate in moyo, this can directly return a reference.
    pub(crate) fn lattice(&self) -> Lattice {
        Lattice {
            inner: self.inner.lattice.clone(),
        }
    }

    pub(crate) fn positions(&self) -> Vec<[f64; 3]> {
        self.inner
            .positions
            .iter()
            .map(|p| [p[0], p[1], p[2]])
            .collect()
    }

    pub(crate) fn numbers(&self) -> &[i32] {
        &self.inner.numbers
    }
}

/// I adapt what I did for the crystal builder pattern in ccmat. The builder pattern with compile
/// time state check can avoid mistake for data initialization.
pub(crate) struct CellBuilder<LatticeSetState, PositionsSetState, NumbersSetState> {
    cell: Cell,
    _lattice: std::marker::PhantomData<LatticeSetState>,
    _positions: std::marker::PhantomData<PositionsSetState>,
    _numbers: std::marker::PhantomData<NumbersSetState>,
}

impl Default for CellBuilder<LatticeNotSet, PositionsNotSet, NumbersNotSet> {
    fn default() -> Self {
        let init_moyo_cell = moyo::base::Cell {
            // simply allocate in memory for futher init
            lattice: moyo::base::Lattice::from_basis([
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0],
            ]),
            positions: vec![],
            numbers: vec![],
        };
        Self {
            cell: Cell {
                inner: init_moyo_cell,
            },
            _lattice: std::marker::PhantomData,
            _positions: std::marker::PhantomData,
            _numbers: std::marker::PhantomData,
        }
    }
}

impl CellBuilder<LatticeNotSet, PositionsNotSet, NumbersNotSet> {
    #[must_use]
    pub(crate) fn new() -> Self {
        Self::default()
    }
}

impl<P, N> CellBuilder<LatticeNotSet, P, N> {
    /// Passing to initialize the lattice.
    /// transfer the ownership, API user need to explicitly clone in order to use lattice
    /// independent from the construct `Cell`.
    #[must_use]
    pub(crate) fn with_lattice(self, lattice: Lattice) -> CellBuilder<LatticeSet, P, N> {
        let inner = moyo::base::Cell {
            lattice: lattice.inner,
            ..self.cell.inner
        };
        CellBuilder {
            cell: Cell { inner },
            _lattice: std::marker::PhantomData,
            _positions: std::marker::PhantomData,
            _numbers: std::marker::PhantomData,
        }
    }
}

impl<L, N> CellBuilder<L, PositionsNotSet, N> {
    /// Passing to initialize the positions.
    /// transfer the ownership, API user need to explicitly clone in order to use positions vec
    /// independent from the construct `Cell`.
    #[must_use]
    pub(crate) fn with_positions(
        self,
        positions: Vec<[f64; 3]>,
    ) -> CellBuilder<L, PositionsSet, N> {
        let positions = positions
            .into_iter()
            .map(|p| nalgebra::Vector3::new(p[0], p[1], p[2]))
            .collect();
        let inner = moyo::base::Cell {
            positions,
            ..self.cell.inner
        };
        CellBuilder {
            cell: Cell { inner },
            _lattice: std::marker::PhantomData,
            _positions: std::marker::PhantomData,
            _numbers: std::marker::PhantomData,
        }
    }
}

impl<L, P> CellBuilder<L, P, NumbersNotSet> {
    /// Passing to initialize the numbers.
    /// transfer the ownership, API user need to explicitly clone in order to use numbers vec
    /// independent from the construct `Cell`.
    #[must_use]
    pub(crate) fn with_numbers(self, numbers: Vec<i32>) -> CellBuilder<L, P, NumbersSet> {
        let inner = moyo::base::Cell {
            numbers,
            ..self.cell.inner
        };
        CellBuilder {
            cell: Cell { inner },
            _lattice: std::marker::PhantomData,
            _positions: std::marker::PhantomData,
            _numbers: std::marker::PhantomData,
        }
    }
}

impl CellBuilder<LatticeSet, PositionsSet, NumbersSet> {
    /// build the final cell
    ///
    /// I assume moyo for performance doesn't do validation on the structure. To eliminate the ill
    /// defined structure for instance where the length of positions and numbers are not equal.
    pub(crate) fn build(self) -> Cell {
        self.cell
    }
}

#[derive(Debug)]
pub struct MoyoError {
    inner: moyo::base::MoyoError,
}

impl std::fmt::Display for MoyoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:?}", self.inner)
    }
}

impl std::error::Error for MoyoError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        Some(&self.inner)
    }
}

impl From<moyo::base::MoyoError> for MoyoError {
    fn from(value: moyo::base::MoyoError) -> Self {
        Self { inner: value }
    }
}

pub trait NiggliReduce: Sized {
    fn niggli_reduce(&self) -> Result<Self, MoyoError>;
}

/// Wrapper of `MoyoDataset` with handy (in my opinion) APIs.
pub(crate) struct SymmetryInfo {
    inner: MoyoDataset,
}

impl SymmetryInfo {
    /// Space group number (1-230)
    ///
    /// # Panics
    /// Panic if space group number return from moyo is negative, i32 -> u32 fail, should be a bug
    /// then.
    #[must_use]
    pub(crate) fn spg_number(&self) -> u32 {
        self.inner
            .number
            .try_into()
            .expect("spage group number not in 1..=230")
    }

    /// Hall symbol number (1-530).
    ///
    /// # Panics
    /// Panic if hall number return from moyo is negative, i32 -> u32 fail, should be a bug
    /// then.
    #[must_use]
    pub(crate) fn hall_number(&self) -> u32 {
        self.inner
            .hall_number
            .try_into()
            .expect("spage group number not in 1..=230")
    }

    /// Bravais class
    ///
    /// # Panics
    /// When moyo failed to get the hall symbol from the `hall_number`, shouldn't happened in ccmat
    /// since the `hall_number` is computed by moyo, if panic it is a bug.
    #[must_use]
    pub(crate) fn bravais_class(&self) -> BravaisClass {
        let hall_number = self.inner.hall_number;
        let hall_symbol =
            hall_symbol_entry(hall_number).expect("unable to get hall symbol from hall_number");
        let arithmetic_entry =
            arithmetic_crystal_class_entry(hall_symbol.arithmetic_number).unwrap();

        BravaisClass {
            inner: arithmetic_entry.bravais_class,
        }
    }

    /// Centering
    ///
    /// # Panics
    /// When moyo failed to get the hall symbol from the `hall_number`, shouldn't happened in ccmat
    /// since the `hall_number` is computed by moyo, if panic it is a bug.
    #[allow(dead_code)]
    #[must_use]
    pub(crate) fn centring(&self) -> Centering {
        let hall_number = self.inner.hall_number;
        let hall_symbol =
            hall_symbol_entry(hall_number).expect("unable to get hall symbol from hall_number");
        Centering {
            inner: hall_symbol.centering,
        }
    }

    /// Hall symbol
    ///
    /// # Panics
    /// When moyo failed to get the hall symbol from the `hall_number`, shouldn't happened in ccmat
    /// since the `hall_number` is computed by moyo, if panic it is a bug.
    #[must_use]
    pub(crate) fn hall_symbol(&self) -> Cow<'_, str> {
        let hall_number = self.inner.hall_number;
        let hall_symbol =
            hall_symbol_entry(hall_number).expect("unable to get hall symbol from hall_number");
        Cow::Borrowed(hall_symbol.hall_symbol)
    }

    /// Check if contain inversion symmetry.
    #[must_use]
    pub(crate) fn has_inversion(&self) -> bool {
        let hall_symbol = self.hall_symbol();
        hall_symbol.starts_with('-')
    }

    /// Standard Cell
    pub(crate) fn std_cell(&self) -> Cell {
        Cell {
            // clone because std_cell can be not the original cell.
            inner: self.inner.std_cell.clone(),
        }
    }
}

/// analyze symmetry
///
/// # Errors
/// ???
// TODO: thiserror
pub(crate) fn analyze_symmetry(
    cell: &Cell,
    symprec: f64,
) -> Result<SymmetryInfo, Box<dyn std::error::Error + Send + Sync>> {
    let inner = MoyoDataset::with_default(&cell.inner, symprec)?;
    let sym_info = SymmetryInfo { inner };
    Ok(sym_info)
}

type TransformMatrix = [[i32; 3]; 3];
type Basis = [[f64; 3]; 3];

// I expect there are three methods: niggli_reduce, niggli_reduce_unchecked, and is_niggli_reduced.
// In ccmat, I may take the approach of current moyo by only providing them as the methods of
// Lattice. (TBD) Meanwhile, in ccmat, I only provide niggli_reduce with always validating the
// result.
// The reason is I recon moyo is more basic lib to provide the extensive operations.
//
// I expect the basis has the same type as input for `Lattice::from_basis()`.
#[allow(dead_code)]
pub(crate) fn niggli_reduce_unchecked(basis: Basis) -> (Basis, TransformMatrix) {
    // NOTE: there is implementation at ``moyo::math::niggli::niggli_reduce`` but not exposed
    // I have to use Lattice::niggli_reduce`` as work around.
    let lattice = moyo::base::Lattice::from_basis(basis);
    let (lattice_converted, mat_trans) = lattice.unchecked_niggli_reduce();
    (lattice_converted.basis.into(), mat_trans.into())
}

pub(crate) fn niggli_reduce(basis: Basis) -> Result<(Basis, TransformMatrix), MoyoError> {
    let lattice = moyo::base::Lattice::from_basis(basis);
    let (lattice_converted, mat_trans) = lattice.niggli_reduce().map_err(MoyoError::from)?;

    Ok((lattice_converted.basis.into(), mat_trans.into()))
}

#[allow(dead_code)]
pub(crate) fn is_niggli_reduced(basis: Basis) -> bool {
    let lattice = moyo::base::Lattice::from_basis(basis);
    lattice.is_niggli_reduced()
}
