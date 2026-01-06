/*
 * structure.rs contains of basic structure utils.
 * notes:
 * Operation safety is guranteed by the type.
 * Use Angstrom as the major internal and default API unit to be consistent with xx/xx.
 * Internally use ``FracCoord`` to represent the position to make sure fractional
 * coords don’t change if lattice changes shape or scale.
 *
 * - Use 'lattice'
 *
 * Compile time build errors include:
 * - when fractional coordinates x not satisfy 0 <= x < 1.0.
 * - multi-set of lattice and sites.
 * - not set lattice or sites.
 *
 * Following errors or runtime validation.
 * - exact duplicate sites (? this might be suitable as compile time error, but how?)
 * - lattice vectors on the same plane
 *
 */

use std::ops::Add;

use crate::math::{TransformationMatrix, Vector3};
use crate::symbol_from_atomic_number;

// TODO: naming convention for vars, check IUCr or cif specification
// Give a table to compare in between different popular tools.

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Angstrom(pub f64);

impl From<Angstrom> for f64 {
    fn from(value: Angstrom) -> Self {
        value.0
    }
}

impl From<f64> for Angstrom {
    fn from(value: f64) -> Self {
        Angstrom(value)
    }
}

impl Add<Angstrom> for Angstrom {
    type Output = Angstrom;

    fn add(self, rhs: Angstrom) -> Self::Output {
        Angstrom::from(f64::from(self) + f64::from(rhs))
    }
}

/// Unit the inverse Angstrom 1/Å for vectors in reciprocal space.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct InvAngstrom(pub f64);

impl Add<InvAngstrom> for InvAngstrom {
    type Output = InvAngstrom;

    fn add(self, rhs: InvAngstrom) -> Self::Output {
        InvAngstrom::from(f64::from(self) + f64::from(rhs))
    }
}

impl From<InvAngstrom> for f64 {
    fn from(value: InvAngstrom) -> Self {
        value.0
    }
}

impl From<f64> for InvAngstrom {
    fn from(value: f64) -> Self {
        InvAngstrom(value)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bohr(pub f64);

impl From<Bohr> for f64 {
    fn from(value: Bohr) -> Self {
        value.0
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FracCoord(pub f64);

impl std::fmt::Display for FracCoord {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{:15.9}", self.0)
    }
}

impl From<FracCoord> for f64 {
    fn from(value: FracCoord) -> Self {
        value.0
    }
}

impl From<f64> for FracCoord {
    fn from(value: f64) -> Self {
        FracCoord(value)
    }
}

impl Vector3<FracCoord> {
    /// m is the matrix transform the coordinates, use inv(m) to transform the vector.
    ///
    /// # Errors
    ///
    /// error if the det of the transformation matrix is 0, non-invertible or singular.
    pub fn change_basis_by(
        &self,
        m: &TransformationMatrix,
    ) -> Result<Self, Box<dyn std::error::Error + Send + Sync>> {
        let x: f64 = self[0].into();
        let y: f64 = self[1].into();
        let z: f64 = self[2].into();

        // m inv
        let m = m
            .inv()
            .ok_or("singular transformation matrix, det(m) ~ 0")?;

        // (x', y', z') = (x, y, z) * m_inv;
        let new_x = x * m[0][0] + y * m[0][1] + z * m[0][2];
        let new_y = x * m[1][0] + y * m[1][1] + z * m[1][2];
        let new_z = x * m[2][0] + y * m[2][1] + z * m[2][2];

        Ok(Vector3([new_x.into(), new_y.into(), new_z.into()]))
    }

    #[must_use]
    pub fn into_cartesian(&self, latt: &Lattice) -> Vector3<Angstrom> {
        let (x, y, z): (f64, f64, f64) = (self[0].into(), self[1].into(), self[2].into());
        x * latt.a() + y * latt.b() + z * latt.c()
    }
}

impl Vector3<Angstrom> {
    #[must_use]
    pub fn into_fraction(&self, latt: &Lattice) -> Vector3<FracCoord> {
        // use reciprocal to avoid redo the inverse etc. might be less intuitive
        let recip = latt.reciprocal();

        // Angstrom * InvAngstrom -> FracCoord
        let coord: Vector3<f64> = Vector3([self[0].into(), self[1].into(), self[2].into()]);
        let a_star = Vector3(recip.a_star().map(f64::from));
        let b_star = Vector3(recip.b_star().map(f64::from));
        let c_star = Vector3(recip.c_star().map(f64::from));

        let factor = 1.0 / (2.0 * std::f64::consts::PI);
        let x = factor * dot(&coord, &a_star);
        let y = factor * dot(&coord, &b_star);
        let z = factor * dot(&coord, &c_star);

        Vector3([FracCoord::from(x), FracCoord::from(y), FracCoord::from(z)])
    }
}

/// `frac!` macro to create `FracCoord` and validate the value is in between [0.0, 1.0)
/// at compile time.
#[macro_export]
macro_rules! frac {
    ($x:expr) => {{
        let frac_coord = $crate::FracCoord($x);
        const {
            assert!(
                (0.0 <= $x && $x < 1.0),
                "invalid fractional coordinate: must satisfy 0.0 <= x < 1.0"
            );
        }
        frac_coord
    }};
}

/// `angstrom!` macro to create `Angstrom`.
#[macro_export]
macro_rules! angstrom {
    ($x:expr) => {{
        let cart_coord = $crate::Angstrom($x);
        cart_coord
    }};
}

/// macro to set the sites (in fraction coordinate)
///
/// # Examples
///
/// ```
/// use ccmat_core::sites_frac_coord;
///
/// let _ = sites_frac_coord![
///     (0.0, 0.0, 0.0), 8;
///     (0.0, 0.0, 0.5), 8;
/// ];
/// ```
#[macro_export]
macro_rules! sites_frac_coord {
    () => {
        Vec::new()
    };
    ( $(
        ($x:expr,$y:expr,$z:expr), $kind:expr
      );+ $(;)?
    ) => {{
        let sites = vec![
            $(
                $crate::SiteFraction::new(
                    $crate::math::Vector3::<$crate::FracCoord>([
                        $crate::frac!($x),
                        $crate::frac!($y),
                        $crate::frac!($z),
                    ]),
                    $kind,
                )
            ),+
        ];
        sites
    }};
}

/// macro to set the sites (in cartesian coordinate in the unit of angstrom)
///
/// # Examples
///
/// ```
/// use ccmat_core::sites_cart_coord;
///
/// let _ = sites_cart_coord![
///     (0.0, 0.0, 0.0), 8;
///     (0.0, 0.0, 3.12), 8;
/// ];
/// ```
#[macro_export]
macro_rules! sites_cart_coord {
    () => {
        Vec::new()
    };
    ( $(
        ($x:expr,$y:expr,$z:expr), $kind:expr
      );+ $(;)?
    ) => {{
        let sites = vec![
            $(
                $crate::SiteCartesian::new(
                    $crate::math::Vector3::<$crate::Angstrom>([
                        $crate::angstrom!($x),
                        $crate::angstrom!($y),
                        $crate::angstrom!($z),
                    ]),
                    $kind,
                )
            ),+
        ];
        sites
    }};
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

pub type Basis = [Vector3<f64>; 3];

/// Lattice
/// inner data structure of the struct are private.
/// TODO: this can derive Copy
#[derive(Debug, Clone)]
pub struct Lattice {
    a: Vector3<Angstrom>,
    b: Vector3<Angstrom>,
    c: Vector3<Angstrom>,
}

/// f64 wrapper for radians
#[derive(Debug, Copy, Clone)]
pub struct Rad(f64);

impl From<Rad> for f64 {
    fn from(value: Rad) -> Self {
        value.0
    }
}

impl From<f64> for Rad {
    fn from(value: f64) -> Self {
        Rad(value)
    }
}

// TODO: I should have a proc-macro for impl all such from f64 traits

/// f64 wrapper for value with unit of volume (Angstrom ^ 2)
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Volume(f64);

impl From<Volume> for f64 {
    fn from(value: Volume) -> Self {
        value.0
    }
}

impl From<f64> for Volume {
    fn from(value: f64) -> Self {
        Volume(value)
    }
}

/// dot product
fn dot(v: &Vector3<f64>, u: &Vector3<f64>) -> f64 {
    v[0] * u[0] + v[1] * u[1] + v[2] * u[2]
}

/// cross product
fn cross(u: &Vector3<f64>, v: &Vector3<f64>) -> Vector3<f64> {
    Vector3::<f64>([
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ])
}

// TODO: Lattice and LatticeReciprocal can be generalized through Unit.
// for instance it is not has trait HasBasis to cover both, need to be more generic
pub trait HasBasis {
    fn basis(&self) -> Basis;
}

impl HasBasis for Lattice {
    fn basis(&self) -> Basis {
        [self.a().into(), self.b().into(), self.c().into()]
    }
}

impl HasBasis for LatticeReciprocal {
    fn basis(&self) -> Basis {
        [
            self.a_star().into(),
            self.b_star().into(),
            self.c_star().into(),
        ]
    }
}

impl From<Basis> for Lattice {
    fn from(bs: Basis) -> Self {
        let a: Vector3<Angstrom> = bs[0].into();
        let b: Vector3<Angstrom> = bs[1].into();
        let c: Vector3<Angstrom> = bs[2].into();
        Self { a, b, c }
    }
}

impl From<Basis> for LatticeReciprocal {
    fn from(bs: Basis) -> Self {
        let a: Vector3<InvAngstrom> = bs[0].into();
        let b: Vector3<InvAngstrom> = bs[1].into();
        let c: Vector3<InvAngstrom> = bs[2].into();
        Self { a, b, c }
    }
}

impl Lattice {
    // TODO: how to use type system to validate the row/column definition?
    #[must_use]
    pub fn new(a: Vector3<Angstrom>, b: Vector3<Angstrom>, c: Vector3<Angstrom>) -> Self {
        Lattice { a, b, c }
    }

    /// Constructs the lattice from Cartesian lattice vectors expressed in angstroms.
    ///
    /// The input is a 3×3 array where each row represents a lattice vector
    /// `(a, b, c)` in Cartesian coordinates, with components given in angstroms.
    /// Each component is converted into the internal `Angstrom` unit type.
    ///
    /// # Parameters
    ///
    /// - `latt`: A 3×3 array of lattice vectors in angstroms.
    ///
    /// # Example
    ///
    /// ```
    /// use ccmat_core::Lattice;
    ///
    /// let latt = [
    ///     [5.43, 0.0, 0.0],
    ///     [0.0, 5.43, 0.0],
    ///     [0.0, 0.0, 5.43],
    /// ];
    ///
    /// let latt = Lattice::from_angstroms(latt);
    /// ```
    #[must_use]
    pub fn from_angstroms(latt: [[f64; 3]; 3]) -> Self {
        let a: Vector3<Angstrom> = Vector3(latt[0].map(Angstrom::from));
        let b: Vector3<Angstrom> = Vector3(latt[1].map(Angstrom::from));
        let c: Vector3<Angstrom> = Vector3(latt[2].map(Angstrom::from));
        Self { a, b, c }
    }

    #[must_use]
    pub fn a(&self) -> Vector3<Angstrom> {
        self.a
    }

    #[must_use]
    pub fn b(&self) -> Vector3<Angstrom> {
        self.b
    }

    #[must_use]
    pub fn c(&self) -> Vector3<Angstrom> {
        self.c
    }

    pub fn lattice_params(&self) -> (Angstrom, Angstrom, Angstrom, Rad, Rad, Rad) {
        let va = self.a.map(f64::from);
        let length_a = f64::sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);

        let vb = self.b.map(f64::from);
        let length_b = f64::sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);

        let vc = self.c.map(f64::from);
        let length_c = f64::sqrt(vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);

        let cos_alpha = (vb[0] * vc[0] + vb[1] * vc[1] + vb[2] * vc[2]) / (length_b * length_c);
        let cos_beta = (va[0] * vc[0] + va[1] * vc[1] + va[2] * vc[2]) / (length_a * length_c);
        let cos_gamma = (va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2]) / (length_a * length_b);

        (
            length_a.into(),
            length_b.into(),
            length_c.into(),
            cos_alpha.acos().into(),
            cos_beta.acos().into(),
            cos_gamma.acos().into(),
        )
    }

    #[must_use]
    pub fn length_a(&self) -> Angstrom {
        self.lattice_params().0
    }

    #[must_use]
    pub fn length_b(&self) -> Angstrom {
        self.lattice_params().1
    }

    #[must_use]
    pub fn length_c(&self) -> Angstrom {
        self.lattice_params().2
    }

    #[must_use]
    pub fn rad_alpha(&self) -> Rad {
        self.lattice_params().3
    }

    #[must_use]
    pub fn rad_beta(&self) -> Rad {
        self.lattice_params().4
    }

    #[must_use]
    pub fn rad_gamma(&self) -> Rad {
        self.lattice_params().5
    }

    #[must_use]
    pub fn volume(&self) -> Volume {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());

        // a⋅(b×c)
        Volume(dot(&a, &cross(&b, &c)))
    }

    #[must_use]
    pub fn reciprocal(&self) -> LatticeReciprocal {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());
        let volume: f64 = Volume(dot(&a, &cross(&b, &c))).into();
        let a_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&b, &c);
        let b_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&c, &a);
        let c_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&a, &b);

        let a_star = Vector3(a_star.map(InvAngstrom::from));
        let b_star = Vector3(b_star.map(InvAngstrom::from));
        let c_star = Vector3(c_star.map(InvAngstrom::from));

        LatticeReciprocal::new(a_star, b_star, c_star)
    }

    /// Lattice is represented in the new basis
    ///
    /// (a', b', c') = (a, b, c) * { m00 m01 m02 }
    ///                            { m10 m11 m12 }
    ///                            { m20 m21 m22 }
    ///
    /// a' = m00 * a + m10 * b + m20 * c
    /// b' = m01 * a + m11 * b + m21 * c
    /// c' = m02 * a + m12 * b + m22 * c
    #[must_use]
    pub fn change_basis_by(&self, m: &TransformationMatrix) -> Self {
        let (a, b, c): (Vector3<f64>, Vector3<f64>, Vector3<f64>) =
            (self.a.into(), self.b.into(), self.c.into());
        let ap = m[0][0] * a + m[1][0] * b + m[2][0] * c;
        let bp = m[0][1] * a + m[1][1] * b + m[2][1] * c;
        let cp = m[0][2] * a + m[1][2] * b + m[2][2] * c;

        Self::new(ap.into(), bp.into(), cp.into())
    }
}

// TODO: add lattice_bohr!()
// TODO: impl Display to Lattice for pretty print.

/// Create a [`Lattice`] using vectors expressed in **Ångström** units.
///
/// This macro constructs a [`Lattice`] from three lattice vectors (`a`, `b`, and `c`)
/// expressed as tuples or arrays of three floating-point numbers.
///
/// - Each component is converted to ``Vector3<Angstrom>`` automatically.
/// - Both `(x, y, z)` and `[x, y, z]` tuple/array syntax are supported for vector.
/// - Trailing commas are optional.
///
/// It supports both **named** and **positional** forms:
///
/// - **Named form** (explicit `a=`, `b=`, `c=`):
///   ```
///   use ccmat_core::lattice_angstrom;
///
///   let latt = lattice_angstrom!(
///       a = (1.0, 0.0, 0.0),
///       b = (0.0, 1.0, 0.0),
///       c = (0.0, 0.0, 1.0),
///   );
///   ```
///
/// - **Positional form** (omit names, ordered as `a`, `b`, `c`):
///   ```
///   use ccmat_core::{lattice_angstrom, Lattice, Angstrom};
///
///   let latt = lattice_angstrom!(
///       (1.0, 0.0, 0.0),
///       (0.0, 1.0, 0.0),
///       (0.0, 0.0, 1.0),
///   );
///   ```
///
/// # Errors
/// - None at compile time; this macro expands directly to constructor calls.
///
/// # Example
///
/// It is also okay to use '[]' instead of '()' for each vector.
///
/// ```
/// use ccmat_core::{lattice_angstrom, Lattice, Angstrom};
///
/// let latt = lattice_angstrom!(
///     a = [2.5, 0.0, 0.0],
///     b = [0.0, 2.5, 0.0],
///     c = [0.0, 0.0, 2.5],
/// );
/// println!("{:?}", latt);
/// ```
#[macro_export]
macro_rules! lattice_angstrom {
    (
        a = $a:tt,
        b = $b:tt,
        c = $c:tt $(,)?
    ) => {{
        macro_rules! __vec3_angstrom {
            ([$x:expr, $y:expr, $z:expr]) => {
                $crate::math::Vector3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
            (($x:expr, $y:expr, $z:expr)) => {
                $crate::math::Vector3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
        }

        let lattice = $crate::Lattice::new(
            __vec3_angstrom!($a),
            __vec3_angstrom!($b),
            __vec3_angstrom!($c),
        );
        lattice
    }};
    (
        $a:tt,
        $b:tt,
        $c:tt $(,)?
    ) => {{
        macro_rules! __vec3_angstrom {
            ([$x:expr, $y:expr, $z:expr]) => {
                $crate::math::Vector3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
            (($x:expr, $y:expr, $z:expr)) => {
                $crate::math::Vector3::<$crate::Angstrom>([
                    $crate::Angstrom($x),
                    $crate::Angstrom($y),
                    $crate::Angstrom($z),
                ])
            };
        }
        let lattice = $crate::Lattice::new(
            __vec3_angstrom!($a),
            __vec3_angstrom!($b),
            __vec3_angstrom!($c),
        );
        lattice
    }};
}

pub struct LatticeReciprocal {
    // internal use a not a_star, but the API is a_star to make it very explicit.
    a: Vector3<InvAngstrom>,
    b: Vector3<InvAngstrom>,
    c: Vector3<InvAngstrom>,
}

impl LatticeReciprocal {
    #[must_use]
    pub fn new(
        a_star: Vector3<InvAngstrom>,
        b_star: Vector3<InvAngstrom>,
        c_star: Vector3<InvAngstrom>,
    ) -> Self {
        LatticeReciprocal {
            a: a_star,
            b: b_star,
            c: c_star,
        }
    }

    #[must_use]
    pub fn reciprocal(&self) -> Lattice {
        let (a, b, c) = (self.a.into(), self.b.into(), self.c.into());
        let volume: f64 = Volume(dot(&a, &cross(&b, &c))).into();
        let a_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&b, &c);
        let b_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&c, &a);
        let c_star = 1.0 / volume * (2.0 * std::f64::consts::PI) * cross(&a, &b);

        let a_star = Vector3(a_star.map(Angstrom::from));
        let b_star = Vector3(b_star.map(Angstrom::from));
        let c_star = Vector3(c_star.map(Angstrom::from));

        Lattice::new(a_star, b_star, c_star)
    }

    pub fn lattice_params(&self) -> (InvAngstrom, InvAngstrom, InvAngstrom, Rad, Rad, Rad) {
        let va = self.a.map(f64::from);
        let length_a = f64::sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);

        let vb = self.b.map(f64::from);
        let length_b = f64::sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);

        let vc = self.c.map(f64::from);
        let length_c = f64::sqrt(vc[0] * vc[0] + vc[1] * vc[1] + vc[2] * vc[2]);

        let cos_alpha = (vb[0] * vc[0] + vb[1] * vc[1] + vb[2] * vc[2]) / (length_b * length_c);
        let cos_beta = (va[0] * vc[0] + va[1] * vc[1] + va[2] * vc[2]) / (length_a * length_c);
        let cos_gamma = (va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2]) / (length_a * length_b);

        (
            length_a.into(),
            length_b.into(),
            length_c.into(),
            f64::acos(cos_alpha).into(),
            f64::acos(cos_beta).into(),
            f64::acos(cos_gamma).into(),
        )
    }

    /// Lattice is represented in the new basis
    ///
    /// (a', b', c') = (a, b, c) * { m00 m01 m02 }
    ///                            { m10 m11 m12 }
    ///                            { m20 m21 m22 }
    ///
    /// a' = m00 * a + m10 * b + m20 * c
    /// b' = m01 * a + m11 * b + m21 * c
    /// c' = m02 * a + m12 * b + m22 * c
    #[must_use]
    pub fn change_basis_by(&self, m: &TransformationMatrix) -> Self {
        let (a, b, c): (Vector3<f64>, Vector3<f64>, Vector3<f64>) =
            (self.a.into(), self.b.into(), self.c.into());
        let ap = m[0][0] * a + m[1][0] * b + m[2][0] * c;
        let bp = m[0][1] * a + m[1][1] * b + m[2][1] * c;
        let cp = m[0][2] * a + m[1][2] * b + m[2][2] * c;
        Self::new(ap.into(), bp.into(), cp.into())
    }

    #[must_use]
    pub fn a_star(&self) -> Vector3<InvAngstrom> {
        self.a
    }

    #[must_use]
    pub fn b_star(&self) -> Vector3<InvAngstrom> {
        self.b
    }

    #[must_use]
    pub fn c_star(&self) -> Vector3<InvAngstrom> {
        self.c
    }
}

#[derive(Debug)]
pub struct MoleculeValidateError {
    message: String,
}

impl std::fmt::Display for MoleculeValidateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for MoleculeValidateError {}

#[derive(Debug)]
pub struct CrystalValidateError {
    message: String,
}

impl std::fmt::Display for CrystalValidateError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}", self.message)
    }
}

impl std::error::Error for CrystalValidateError {}

pub struct LatticeSet;
pub struct LatticeNotSet;
pub struct SitesSet;
pub struct SitesNotSet;

/// use builder pattern so the validation is runtime
///
/// # Example
/// ```
/// use ccmat_core::*;
///
/// let lattice = lattice_angstrom![
///     a = (1.0, 0.0, 0.0),
///     b = (0.0, 1.0, 0.0),
///     c = (0.0, 0.0, 1.0),
/// ];
/// let sites = vec![];
/// let crystal = CrystalBuilder::new()
///     .with_lattice(&lattice)
///     .with_frac_sites(sites)
///     .build()
///     .unwrap();
/// ```
#[derive(Debug)]
pub struct CrystalBuilder<LatticeSetState, SiteSetState> {
    crystal: Crystal,
    _lattice: std::marker::PhantomData<LatticeSetState>,
    _sites: std::marker::PhantomData<SiteSetState>,
}

impl Default for CrystalBuilder<LatticeNotSet, SitesNotSet> {
    fn default() -> Self {
        Self {
            crystal: Crystal {
                lattice: lattice_angstrom!([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0],),
                positions: vec![],
                species: vec![],
            },
            _lattice: std::marker::PhantomData,
            _sites: std::marker::PhantomData,
        }
    }
}

impl CrystalBuilder<LatticeNotSet, SitesNotSet> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

impl<S> CrystalBuilder<LatticeNotSet, S> {
    // TODO: should Lattice pass as ref?
    //
    /// Set the lattice for the crystal, where the lattice can be defined using `lattice_angstrom!`
    /// macro.
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmat_core::*;
    /// let lattice = lattice_angstrom![
    ///     a = (1.0, 0.0, 0.0),
    ///     b = (0.0, 1.0, 0.0),
    ///     c = (0.0, 0.0, 1.0),
    /// ];
    /// ```
    #[must_use]
    pub fn with_lattice(self, lattice: &Lattice) -> CrystalBuilder<LatticeSet, S> {
        CrystalBuilder {
            crystal: Crystal {
                // TODO: transfer the ownership instead clone?
                lattice: lattice.clone(),
                ..self.crystal
            },
            _lattice: std::marker::PhantomData,
            _sites: std::marker::PhantomData,
        }
    }
}

impl<L> CrystalBuilder<L, SitesNotSet> {
    /// Set the sites (in fraction coordinate) for the crystal, where sites can into an iterator of sites:
    ///
    /// # Examples
    /// ```
    /// use ccmat_core::*;
    /// let sites = sites_frac_coord![
    ///     (0.6666666666666667, 0.3333333333333333, 0.8333333333333333), atomic_number!(In);
    ///     (0.3333333333333333, 0.6666666666666666, 0.1666666666666666), atomic_number!(In);
    ///     (0.0000000000000000, 0.0000000000000000, 0.4999999999999999), atomic_number!(In);
    ///     (0.0000000000000000, 0.0000000000000000, 0.0000000000000000), atomic_number!(Hg);
    ///     (0.6666666666666666, 0.3333333333333333, 0.3333333333333333), atomic_number!(Hg);
    ///     (0.3333333333333333, 0.6666666666666666, 0.6666666666666666), atomic_number!(Hg);
    /// ];
    /// ```
    #[must_use]
    pub fn with_frac_sites<I>(self, sites: I) -> CrystalBuilder<L, SitesSet>
    where
        I: IntoIterator<Item = SiteFraction>,
    {
        let (positions, species) = sites
            .into_iter()
            .map(|atom| (atom.position, atom.specie.clone()))
            .collect();

        CrystalBuilder {
            crystal: Crystal {
                positions,
                species,
                ..self.crystal
            },
            _lattice: std::marker::PhantomData,
            _sites: std::marker::PhantomData,
        }
    }
}

impl CrystalBuilder<LatticeSet, SitesNotSet> {
    /// Different from `with_frac_sites`, the `with_cart_sites` can only be called after `with_lattice`
    /// to set lattice first.
    ///
    /// Set the sites (in cartesian coordinate) for the crystal, where sites can into an iterator of sites:
    ///
    /// # Examples
    /// ```
    /// use ccmat_core::*;
    /// let sites = sites_frac_coord![
    ///     (0.6666666666666667, 0.3333333333333333, 0.8333333333333333), atomic_number!(In);
    ///     (0.3333333333333333, 0.6666666666666666, 0.1666666666666666), atomic_number!(In);
    ///     (0.0000000000000000, 0.0000000000000000, 0.4999999999999999), atomic_number!(In);
    ///     (0.0000000000000000, 0.0000000000000000, 0.0000000000000000), atomic_number!(Hg);
    ///     (0.6666666666666666, 0.3333333333333333, 0.3333333333333333), atomic_number!(Hg);
    ///     (0.3333333333333333, 0.6666666666666666, 0.6666666666666666), atomic_number!(Hg);
    /// ];
    /// ```
    #[must_use]
    pub fn with_cart_sites<I>(self, sites: I) -> CrystalBuilder<LatticeSet, SitesSet>
    where
        I: IntoIterator<Item = SiteCartesian>,
    {
        let latt = self.crystal.lattice();
        let (positions, species) = sites
            .into_iter()
            .map(|site| {
                let pos_cart = site.position;
                let pos = pos_cart.into_fraction(&latt);
                (pos, site.specie)
            })
            .collect();

        CrystalBuilder {
            crystal: Crystal {
                positions,
                species,
                ..self.crystal
            },
            _lattice: std::marker::PhantomData,
            _sites: std::marker::PhantomData,
        }
    }
}

impl CrystalBuilder<LatticeSet, SitesSet> {
    fn validate(&self) -> Result<(), CrystalValidateError> {
        // TODO: call validate
        if !self.crystal.positions.len() == self.crystal.species.len() {
            return Err(CrystalValidateError {
                message: "crystal valid failed".to_string(),
            });
        }
        Ok(())
    }

    // build without runtime validation this is for proc macro which valid in compile time.
    // It is also used inside crate where the crystal is known to be valid.
    #[must_use]
    pub fn build_uncheck(self) -> Crystal {
        debug_assert!(self.crystal.positions.len() == self.crystal.species.len());

        self.crystal
    }

    /// build and validate the it is a valid crystal.
    /// At the moment only validate the size(positions) == size(numbers)
    ///
    /// # Errors
    /// ??
    pub fn build(self) -> Result<Crystal, CrystalValidateError> {
        self.validate()?;

        let crystal = self.build_uncheck();

        Ok(crystal)
    }
}

// TODO: partial occupation on sites
#[derive(Debug, Clone)]
pub struct Specie {
    atomic_number: u8,
}

impl Specie {
    fn new(atomic_number: u8) -> Self {
        Specie { atomic_number }
    }

    #[must_use]
    pub fn atomic_number(&self) -> u8 {
        self.atomic_number
    }

    #[must_use]
    pub fn symbol(&self) -> String {
        symbol_from_atomic_number(self.atomic_number()).expect("not a valid atomic number")
    }
}

// TODO: would be helpful to have a macro for initialize site
// something like: site_frac!([0.0, 0.5, 0.5], Si)

/// A crystallographic site.
///
/// A `SiteFraction` represents an atomic position within a crystal unit cell,
/// expressed in **fractional coordinates**, together with its chemical
/// identity (`Specie`).
///
/// The `position` is given in fractional coordinates relative to the
/// lattice vectors of the crystal. Each component typically lies in
/// the range `[0, 1)`, but this is not enforced.
///
/// # Examples
///
/// ```
/// use ccmat_core::*;
/// use ccmat_core::math::Vector3;
///
/// let site = SiteFraction::new(
///     Vector3([frac!(0.0), frac!(0.5), frac!(0.5)]),
///     atomic_number!(Si),
/// );
///
/// let pos = site.position();
/// ```
#[derive(Debug)]
pub struct SiteFraction {
    position: Vector3<FracCoord>,
    specie: Specie,
}

impl SiteFraction {
    /// Creates a new crystallographic site.
    ///
    /// Returns a `SiteFraction` with the given position and species.
    #[must_use]
    pub fn new(position: Vector3<FracCoord>, atomic_number: u8) -> Self {
        SiteFraction {
            position,
            specie: Specie::new(atomic_number),
        }
    }

    /// Returns the fractional position of the site.
    ///
    /// The returned vector is expressed in fractional coordinates with respect to the crystal lattice.
    #[must_use]
    pub fn position(&self) -> Vector3<FracCoord> {
        self.position
    }
}

#[derive(Debug)]
pub struct SiteCartesian {
    position: Vector3<Angstrom>,
    specie: Specie,
}

impl SiteCartesian {
    /// Creates a new crystallographic site.
    ///
    /// Returns a `SiteCartesian` with the given position and species.
    #[must_use]
    pub fn new(position: Vector3<Angstrom>, atomic_number: u8) -> Self {
        SiteCartesian {
            position,
            specie: Specie::new(atomic_number),
        }
    }

    /// Returns the cartesian position of the site.
    ///
    /// The returned vector is expressed in cartesian coordinates (in Angstrom) with respect to the crystal lattice.
    #[must_use]
    pub fn position(&self) -> Vector3<Angstrom> {
        self.position
    }
}

#[derive(Debug, Clone)]
pub struct Molecule {
    positions: Vec<Vector3<Angstrom>>,
    species: Vec<Specie>,
}

#[derive(Debug)]
pub struct MoleculeBuilder<SiteSetState> {
    molecule: Molecule,
    _sites: std::marker::PhantomData<SiteSetState>,
}

impl MoleculeBuilder<SitesNotSet> {
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }
}

impl Default for MoleculeBuilder<SitesNotSet> {
    fn default() -> Self {
        Self {
            molecule: Molecule {
                positions: vec![],
                species: vec![],
            },
            _sites: std::marker::PhantomData,
        }
    }
}

impl MoleculeBuilder<SitesNotSet> {
    /// Set the sites (in cartesian coordinate) for the molecue, where sites can into an iterator of sites:
    ///
    /// # Examples
    /// ```
    /// use ccmat_core::*;
    ///
    /// let sites = sites_cart_coord![
    ///     (2.0, 3.568, 0.67), atomic_number!(In);
    ///     (0.78, 0.6677, 0.0989), atomic_number!(Hg);
    /// ];
    /// ```
    #[must_use]
    #[allow(clippy::unused_self)]
    pub fn with_sites<I>(self, sites: I) -> MoleculeBuilder<SitesSet>
    where
        I: IntoIterator<Item = SiteCartesian>,
    {
        let (positions, species) = sites
            .into_iter()
            .map(|atom| (atom.position, atom.specie.clone()))
            .collect();

        MoleculeBuilder {
            molecule: Molecule { positions, species },
            _sites: std::marker::PhantomData,
        }
    }
}

impl MoleculeBuilder<SitesSet> {
    fn validate(&self) -> Result<(), MoleculeValidateError> {
        // TODO: call validate
        if !self.molecule.positions.len() == self.molecule.species.len() {
            return Err(MoleculeValidateError {
                message: "molecule valid failed".to_string(),
            });
        }
        Ok(())
    }

    // build without runtime validation this is for proc macro which valid in compile time.
    // It is also used inside crate where the crystal is known to be valid.
    #[must_use]
    pub fn build_uncheck(self) -> Molecule {
        debug_assert!(self.molecule.positions.len() == self.molecule.species.len());

        self.molecule
    }

    /// build and validate the it is a valid molecule.
    /// At the moment only validate the size(positions) == size(numbers)
    ///
    /// # Errors
    /// ??
    pub fn build(self) -> Result<Molecule, MoleculeValidateError> {
        self.validate()?;

        let mol = self.build_uncheck();

        Ok(mol)
    }
}

impl Molecule {
    /// vec of positions in Cartesian coordinate
    // XXX: I may need a view_positions which return ref, not take ownership
    // XXX: is this a sound API that this is not sync with Crystal's `positions` which return frac
    // coordinates?
    // I am leaning to change crystal's positions to return also Angstrom and rename position ->
    // positions_fraction.
    #[must_use]
    pub fn positions(&self) -> Vec<Vector3<Angstrom>> {
        self.positions.clone()
    }

    #[must_use]
    pub fn species(&self) -> &[Specie] {
        &self.species
    }
}

// Crystal is the public API so should be align with the real world convention.
// I did not expose the data structure for crystal directly but the builder.
// Internally fileds data structures are private to keep API stable.
// Now I try to align it similar to moyo's Cell data structure.
//
//
// TODO: not yet generic, but for 3D only at the moment, generalize when doing 2D and 1D
// in prototype, this struct include as much information as possible. may need to be generic.
// use rust's type system to check the problem of frac coordinates and direct coordinates.
#[derive(Debug, Clone)]
pub struct Crystal {
    lattice: Lattice,
    positions: Vec<Vector3<FracCoord>>,
    species: Vec<Specie>,
}

impl Crystal {
    // XXX: is this method redundant??
    #[must_use]
    pub fn builder() -> CrystalBuilder<LatticeNotSet, SitesNotSet> {
        CrystalBuilder::new()
    }

    #[must_use]
    pub fn lattice(&self) -> Lattice {
        // TODO: is Cow possible?
        self.lattice.clone()
    }

    #[must_use]
    pub fn lattice_reciprocal(&self) -> LatticeReciprocal {
        self.lattice().reciprocal()
    }

    /// vec of positions in fractional coordinate
    #[must_use]
    pub fn positions_fraction(&self) -> Vec<Vector3<FracCoord>> {
        // TODO: avoid clone in readonly?
        self.positions.clone()
    }

    /// vec of positions in Cartesian coordinate
    #[must_use]
    pub fn positions(&self) -> Vec<Vector3<Angstrom>> {
        // TODO: avoid clone in readonly?
        self.positions
            .iter()
            .map(|p| p.into_cartesian(&self.lattice()))
            .collect()
    }

    #[must_use]
    pub fn species(&self) -> &[Specie] {
        &self.species
    }

    #[must_use]
    pub fn volume(&self) -> Volume {
        self.lattice.volume()
    }

    /// For positions in fraction coordinates, wrap the positions in the cell.
    /// Because there is periodic in the crystal system.
    pub fn wrap_frac_positions(&mut self) {
        for p in &mut self.positions {
            let p_ = p.map(|i| FracCoord::from(f64::from(i) - f64::from(i).floor()));
            *p = Vector3(p_);
        }
    }
}

// pub fn find_primitive_spglib(
//     crystal: &Crystal,
// ) -> Result<(Crystal, PMatrix, InvPMatrix), Box<dyn std::error::Error + Sync + Send>> {
//     todo!()
// }

#[cfg(test)]
mod tests {
    use crate::{atomic_number, matrix_3x3};

    use super::*;

    // TODO: dup (search the same macro name in the repo), move to test_utils?
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
    fn build_crystal_compile_error() {
        let t = trybuild::TestCases::new();
        t.compile_fail("tests/build_crystal/^fail_*.rs");
    }

    #[test]
    fn macro_lattice_angstrom() {
        let _ = lattice_angstrom![
            a = (1.0, 0.0, 0.0),
            b = (0.0, 1.0, 0.0),
            c = (0.0, 0.0, 1.0),
        ];
        let _ = lattice_angstrom![(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0),];
        let _ = lattice_angstrom![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0],];
        // trailing comma ','
        let _ = lattice_angstrom![[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    }

    #[test]
    fn macro_sites_frac() {
        let _: Vec<SiteFraction> = sites_frac_coord![];
        let _ = sites_frac_coord![
            (0.0, 0.0, 0.0), 8;
            (0.0, 0.0, 0.5), 8;
        ];
    }

    #[test]
    fn reciprocal() {
        let lattice = lattice_angstrom![
            // no orthogonal cell
            a = (2.0, 0.5, 0.0),
            b = (0.0, 3.0, 0.5),
            c = (0.5, 0.0, 4.0),
        ];

        let latt2 = lattice.reciprocal().reciprocal();

        // a
        assert_eq_approx!(f64::from(latt2.a[0]), 2.0);
        assert_eq_approx!(f64::from(latt2.a[1]), 0.5);
        assert_eq_approx!(f64::from(latt2.a[2]), 0.0);

        // b
        assert_eq_approx!(f64::from(latt2.b[2]), 0.5);

        // c
        assert_eq_approx!(f64::from(latt2.c[1]), 0.0);
    }

    #[test]
    fn crystal_volume() {
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

        assert_eq!(crystal.volume(), Volume(4.603 * 4.603 * 4.603));
    }

    #[test]
    fn frac_to_cart_round_trip() {
        let latt = lattice_angstrom![
            // no orthogonal cell
            a = (2.0, 0.5, 0.0),
            b = (0.0, 3.0, 0.5),
            c = (0.5, 0.0, 4.0),
        ];

        let pos = Vector3([
            FracCoord::from(0.3),
            FracCoord::from(0.1),
            FracCoord::from(0.28),
        ]);
        let back_pos = pos.into_cartesian(&latt).into_fraction(&latt);
        for i in 0..3 {
            assert_eq_approx!(f64::from(pos[i]), f64::from(back_pos[i]));
        }
    }

    #[test]
    fn latt_change_basis_by() {
        let lattice = lattice_angstrom![
            // no orthogonal cell
            a = (2.0, 0.5, 0.0),
            b = (0.0, 3.0, 0.5),
            c = (0.5, 0.0, 4.0),
        ];

        let tmatrix = matrix_3x3![
            1 0 0;
            0 1 0;
            0 0 1;
        ];

        let latt = lattice.change_basis_by(&tmatrix);
        // a
        assert_eq_approx!(f64::from(latt.a[0]), 2.0);
        assert_eq_approx!(f64::from(latt.a[1]), 0.5);
        assert_eq_approx!(f64::from(latt.a[2]), 0.0);

        // b
        assert_eq_approx!(f64::from(latt.b[2]), 0.5);

        // c
        assert_eq_approx!(f64::from(latt.c[1]), 0.0);

        // rotate in xy plane
        let tmatrix = matrix_3x3![
            sqrt(3.)/2.,    -1./2.,     0;
            1./2.,     sqrt(3.)/2.,     0;
            0,                   0,     1;
        ];

        let latt = lattice.change_basis_by(&tmatrix);

        //  1.73205  1.93301  0.25
        // -1.0      2.34808  0.433013
        //  0.5      0.0      4.0

        // a
        assert_eq_approx!(f64::from(latt.a[0]), 1.73205, 1e-4);
        assert_eq_approx!(f64::from(latt.a[1]), 1.93301, 1e-4);
        assert_eq_approx!(f64::from(latt.a[2]), 0.25, 1e-4);

        // b
        assert_eq_approx!(f64::from(latt.b[2]), 0.433_013, 1e-4);

        // c
        assert_eq_approx!(f64::from(latt.c[1]), 0.0);
    }

    #[test]
    fn position_change_basis_by() {
        let p = Vector3::<FracCoord>([
            FracCoord::from(1.0),
            FracCoord::from(2.0),
            FracCoord::from(3.0),
        ]);
        let m = matrix_3x3![
            0 1 0;
            1 0 0;
            0 0 1;
        ];

        let pt = p.change_basis_by(&m).unwrap();
        assert_eq_approx!(f64::from(pt[0]), 2.0);
        assert_eq_approx!(f64::from(pt[1]), 1.0);
        assert_eq_approx!(f64::from(pt[2]), 3.0);
    }
}
