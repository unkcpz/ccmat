use std::ops::{Add, Deref, Index, Mul};

use crate::structure::{Angstrom, InvAngstrom};
pub use ccmat_macros::matrix_3x3;

#[must_use]
pub fn approx_f64(a: f64, b: f64, tol: f64) -> bool {
    f64::abs(a - b) < tol
}

#[derive(Default, Debug)]
pub struct Matrix3(pub [[f64; 3]; 3]);

impl Matrix3 {
    #[must_use]
    pub fn det(&self) -> f64 {
        // det=a(ei−fh)−b(di−fg)+c(dh−eg)
        let m = self.0;
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }

    // transpose
    #[must_use]
    pub fn t(&self) -> Self {
        let m = self.0;
        let inner = [
            [m[0][0], m[1][0], m[2][0]],
            [m[0][1], m[1][1], m[2][1]],
            [m[0][2], m[1][2], m[2][2]],
        ];
        Matrix3(inner)
    }

    #[allow(clippy::many_single_char_names)]
    #[must_use]
    pub fn inv(&self) -> Option<Self> {
        let m = &self.0;

        let det = self.det();

        if det.abs() < 1e-12 {
            return None; // non-invertible or singular
        }

        let a = m[0][0];
        let b = m[0][1];
        let c = m[0][2];
        let d = m[1][0];
        let e = m[1][1];
        let f = m[1][2];
        let g = m[2][0];
        let h = m[2][1];
        let i = m[2][2];

        let inv_det = 1.0 / det;

        let adj = [
            [(e * i - f * h), -(b * i - c * h), (b * f - c * e)],
            [-(d * i - f * g), (a * i - c * g), -(a * f - c * d)],
            [(d * h - e * g), -(a * h - b * g), (a * e - b * d)],
        ];

        let mut inv = [[0.0; 3]; 3];

        for row in 0..3 {
            for col in 0..3 {
                inv[row][col] = adj[row][col] * inv_det;
            }
        }

        Some(Self(inv))
    }
}

impl Index<usize> for Matrix3 {
    type Output = [f64; 3];

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Mul<Matrix3> for f64 {
    type Output = Matrix3;

    /// Scalar multiply for a `Matrix3`
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmat_core::matrix_3x3;
    ///
    /// let mat3 = matrix_3x3![
    ///     1, 0, 0;
    ///     0, 1, 0;
    ///     1, 0, 1;
    /// ];
    /// let mat3 = 1./2. * mat3;
    /// assert_eq!(mat3[0][0], 1./2.);
    /// assert_eq!(mat3[2][0], 1./2.);
    /// ```
    fn mul(self, rhs: Matrix3) -> Self::Output {
        let v = rhs.0.map(|x| x.map(|x| x * self));
        Matrix3(v)
    }
}

impl Mul<Vector3<f64>> for &Matrix3 {
    type Output = Vector3<f64>;

    /// Left transformation on a vector
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmat_core::matrix_3x3;
    /// use ccmat_core::math::Vector3;
    ///
    /// let rot_mat = matrix_3x3![
    ///      cos(PI/3.), -sin(PI/3.), 0;
    ///      sin(PI/3.), cos(PI/3.), 0;
    ///      0,  0,  1;
    /// ];
    ///
    /// let v: Vector3<f64> = Vector3([1., 0., 0.]);
    /// let new_v = &rot_mat * v;
    ///
    /// assert!(new_v[0] - 0.5 < 1e-9);
    /// assert!(new_v[1] - f64::sqrt(3.)/2. < 1e-9);
    /// assert!(new_v[2] - 0. < 1e-9);
    /// ```
    fn mul(self, rhs: Vector3<f64>) -> Self::Output {
        let m = self;
        let p0 = m[0][0] * rhs[0] + m[0][1] * rhs[1] + m[0][2] * rhs[2];
        let p1 = m[1][0] * rhs[0] + m[1][1] * rhs[1] + m[1][2] * rhs[2];
        let p2 = m[2][0] * rhs[0] + m[2][1] * rhs[1] + m[2][2] * rhs[2];

        Vector3([p0, p1, p2])
    }
}

pub type TransformationMatrix = Matrix3;
pub type RotationMatrix = Matrix3;

#[macro_export]
macro_rules! matrix_3x3 {
    ( $($tokens:tt)* ) => {{
        let inner = $crate::math::matrix_3x3!($($tokens)*);
        $crate::math::Matrix3(inner)
    }};
}

// TODO: add handy idx accessing method
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Vector3<T>(pub [T; 3]);

impl<T> Add<Vector3<T>> for Vector3<T>
where
    T: Add<Output = T> + Copy,
{
    type Output = Self;

    fn add(self, rhs: Vector3<T>) -> Self::Output {
        Vector3([
            (*self)[0] + (*rhs)[0],
            (*self)[1] + (*rhs)[1],
            (*self)[2] + (*rhs)[2],
        ])
    }
}

impl<T> Deref for Vector3<T> {
    type Target = [T; 3];

    fn deref(&self) -> &Self::Target {
        &self.0
    }
}

impl<T> Index<usize> for Vector3<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl Mul<Vector3<f64>> for f64 {
    type Output = Vector3<f64>;

    /// Scalar multiply (on left) for a `Vector3<f64>`
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmat_core::math::Vector3;
    ///
    /// let vec3 = Vector3::<f64>([2.0, 2.0, 4.0]);
    /// assert_eq!(0.1 * vec3, Vector3::<f64>([0.2, 0.2, 0.4]));
    /// ```
    fn mul(self, rhs: Vector3<f64>) -> Self::Output {
        let v = rhs.map(|x| x * self);
        Vector3::<f64>(v)
    }
}

impl Mul<f64> for Vector3<f64> {
    type Output = Vector3<f64>;

    /// Scalar multiply (on right) for a `Vector3<f64>`
    ///
    /// # Examples
    ///
    /// ```
    /// use ccmat_core::math::Vector3;
    ///
    /// let vec3 = Vector3::<f64>([2.0, 2.0, 4.0]);
    /// assert_eq!(vec3 * 0.1, Vector3::<f64>([0.2, 0.2, 0.4]));
    /// ```
    fn mul(self, rhs: f64) -> Self::Output {
        let v = self.map(|x| x * rhs);
        Vector3::<f64>(v)
    }
}

impl From<Vector3<Angstrom>> for Vector3<f64> {
    fn from(v: Vector3<Angstrom>) -> Self {
        Vector3::<f64>(v.map(f64::from))
    }
}

impl From<Vector3<f64>> for Vector3<Angstrom> {
    fn from(v: Vector3<f64>) -> Self {
        Vector3::<Angstrom>(v.map(Angstrom::from))
    }
}

impl From<Vector3<InvAngstrom>> for Vector3<f64> {
    fn from(v: Vector3<InvAngstrom>) -> Self {
        Vector3::<f64>(v.map(f64::from))
    }
}

impl From<Vector3<f64>> for Vector3<InvAngstrom> {
    fn from(v: Vector3<f64>) -> Self {
        Vector3::<InvAngstrom>(v.map(InvAngstrom::from))
    }
}

#[cfg(test)]
mod tests {

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
    fn matrix_3x3() {
        let mat = matrix_3x3![
            cos(PI/3.)   2           3;
            4            5     1.0+6.1;
            7            8   sqrt(9.0);
        ];

        assert_eq_approx!(mat[0][0], 0.5);
        assert_eq_approx!(mat[1][2], 7.1);
    }
}
