#![allow(non_upper_case_globals)]
use evalexpr::{
    eval_float_with_context, ContextWithMutableVariables, DefaultNumericTypes, HashMapContext,
    Value,
};

use ccmat_core::{Angstrom, FracCoord, Rad};

use crate::ExtBravaisClass;

#[allow(non_camel_case_types)]
#[derive(Debug, Hash, Eq, PartialEq)]
pub enum HighSymmetryPoint {
    A,
    A_0,
    B,
    B_0,
    B_2,
    C,
    C_0,
    C_2,
    C_4,
    D,
    DELTA_0,
    D_0,
    D_2,
    E,
    E_0,
    E_2,
    E_4,
    F,
    F_0,
    F_2,
    F_4,
    G,
    #[allow(clippy::upper_case_acronyms)]
    GAMMA,
    G_0,
    G_2,
    G_4,
    G_6,
    H,
    H_0,
    H_2,
    H_4,
    H_6,
    I,
    I_2,
    J_0,
    K,
    K_2,
    K_4,
    L,
    LAMBDA_0,
    L_0,
    L_2,
    L_4,
    M,
    M_0,
    M_2,
    M_4,
    M_6,
    M_8,
    N,
    N_2,
    N_4,
    N_6,
    P,
    P_0,
    P_2,
    Q_0,
    R,
    R_0,
    R_2,
    S,
    SIGMA_0,
    S_0,
    S_2,
    S_4,
    S_6,
    T,
    T_2,
    U,
    U_0,
    U_2,
    V,
    V_0,
    V_2,
    W,
    W_2,
    X,
    X_1,
    Y,
    Y_0,
    Y_2,
    Y_4,
    Z,
    Z_0,
    Z_2,
}

// XXX pub? in order to be tested against original python impl
#[must_use]
pub(crate) fn lookup(ext_bravais: &ExtBravaisClass) -> &'static KpathInfo {
    match ext_bravais {
        // Triclinic
        ExtBravaisClass::aP1 => unimplemented!("aP1 is reserved for aP2+aP3, check hpkot paper"),
        ExtBravaisClass::aP2 => &INFO_aP2,
        ExtBravaisClass::aP3 => &INFO_aP3,
        // Monoclinic
        ExtBravaisClass::mP1 => &INFO_mP1,
        ExtBravaisClass::mC1 => &INFO_mC1,
        ExtBravaisClass::mC2 => &INFO_mC2,
        ExtBravaisClass::mC3 => &INFO_mC3,
        // Orthorhombic
        ExtBravaisClass::oP1 => &INFO_oP1,
        ExtBravaisClass::oA1 => &INFO_oA1,
        ExtBravaisClass::oA2 => &INFO_oA2,
        ExtBravaisClass::oC1 => &INFO_oC1,
        ExtBravaisClass::oC2 => &INFO_oC2,
        ExtBravaisClass::oF1 => &INFO_oF1,
        ExtBravaisClass::oF2 => &INFO_oF2,
        ExtBravaisClass::oF3 => &INFO_oF3,
        ExtBravaisClass::oI1 => &INFO_oI1,
        ExtBravaisClass::oI2 => &INFO_oI2,
        ExtBravaisClass::oI3 => &INFO_oI3,
        // Tetragonal
        ExtBravaisClass::tP1 => &INFO_tP1,
        ExtBravaisClass::tI1 => &INFO_tI1,
        ExtBravaisClass::tI2 => &INFO_tI2,
        // Rhombohedral
        ExtBravaisClass::hR1 => &INFO_hR1,
        ExtBravaisClass::hR2 => &INFO_hR2,
        // Hexagonal
        ExtBravaisClass::hP1 => &INFO_hP1,
        ExtBravaisClass::hP2 => &INFO_hP2,
        // Cubic
        ExtBravaisClass::cP1 => &INFO_cP1,
        ExtBravaisClass::cP2 => &INFO_cP2,
        ExtBravaisClass::cF1 => &INFO_cF1,
        ExtBravaisClass::cF2 => &INFO_cF2,
        ExtBravaisClass::cI1 => &INFO_cI1,
    }
}

#[derive(Debug)]
struct FracCoordExpr(&'static str);

macro_rules! var {
    ($ex:expr) => {
        FracCoordExpr(stringify!($ex))
    };
}

macro_rules! expr {
    ($ex:expr) => {{
        // this is to turn the evaluation always perform as float type.
        let s = concat!("1.0 * ", stringify!($ex));
        FracCoordExpr(s)
    }};
}

/// ```ignore
/// static KPARAM_hR1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
///     D: a * a / 4 / c / c,
///     Y: 5 / 6 - 2 *D,
///     N: 1 / 3 + D,
/// );
/// ```
///
/// become
///
/// ```ignore
/// static KPARAM_hR1: &[(FracCoordExpr, FracCoordExpr)] = &[
///     (var!(D), expr!(a * a / 4 / c / c)),
///     (var!(Y), expr!(5 / 6 - 2 * D)),
///     (var!(N), expr!(1 / 3 + D)),
/// ];
/// ```
macro_rules! kparams {
    () => {{ &[] }};
    ($($lhs:ident : $rhs:expr ,)+$(,)?) => {{
        &[
            $((var!($lhs), expr!($rhs)),)+
        ]
    }};
}

/// ```ignore
///
/// static PATH_hR1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
///     GAMMA T,
///     T H_2,
///     H_0 L,
///     L GAMMA,
///     GAMMA S_0,
///     S_2 F,
///     F GAMMA,
/// ];
/// ```
///
/// become
///
/// ```ignore
/// static PATH_hR1: &[(HighSymmetryPoint, HighSymmetryPoint)] = &[
///     (HighSymmetryPoint::GAMMA, HighSymmetryPoint::T),
///     (HighSymmetryPoint::T, HighSymmetryPoint::H_2),
///     (HighSymmetryPoint::H_0, HighSymmetryPoint::L),
///     (HighSymmetryPoint::L, HighSymmetryPoint::GAMMA),
///     (HighSymmetryPoint::GAMMA, HighSymmetryPoint::S_0),
///     (HighSymmetryPoint::S_2, HighSymmetryPoint::F),
///     (HighSymmetryPoint::F, HighSymmetryPoint::GAMMA),
/// ];
/// ```
macro_rules! kpath {
    ($($src:ident $dst:ident , )+$(,)?) => {{
        &[
            $( (HighSymmetryPoint::$src, HighSymmetryPoint::$dst), )+
        ]
    }};
}

macro_rules! kpoints {
    ( $( $symp:ident : $x:expr , $y:expr , $z:expr ,)+$(,)? ) => {{
        &[
            $( (HighSymmetryPoint::$symp, (expr!($x), expr!($y), expr!($z))), )+
        ]
    }};
}

#[derive(Debug)]
pub struct KpathInfo {
    path: &'static [(HighSymmetryPoint, HighSymmetryPoint)],
    points: &'static [(
        HighSymmetryPoint,
        (FracCoordExpr, FracCoordExpr, FracCoordExpr),
    )],
    kparam: &'static [(FracCoordExpr, FracCoordExpr)],
}

#[derive(Debug)]
pub struct KpathEval {
    path: &'static [(HighSymmetryPoint, HighSymmetryPoint)],
    points: Vec<(
        &'static HighSymmetryPoint,
        (FracCoord, FracCoord, FracCoord),
    )>,
}

impl std::fmt::Display for KpathEval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "---- path ----")?;
        for p in self.path {
            writeln!(f, "{:?} -> {:?}", p.0, p.1)?;
        }
        writeln!(f, "---- ---- ----")?;

        writeln!(f, "---- high symmetry points ----")?;
        for pt in &self.points {
            writeln!(f, "{:?}: {}, {}, {}", pt.0, pt.1 .0, pt.1 .1, pt.1 .2)?;
        }
        writeln!(f, "---- -------------------- ----")?;
        Ok(())
    }
}

pub(crate) fn eval(
    path_info: &KpathInfo,
    lattice_params: (Angstrom, Angstrom, Angstrom, Rad, Rad, Rad),
) -> Result<KpathEval, Box<dyn std::error::Error + Send + Sync>> {
    let (a, b, c, alpha, beta, gamma) = lattice_params;

    let mut context = HashMapContext::<DefaultNumericTypes>::new();
    context.set_value("a".into(), Value::from_float(a.into()))?;
    context.set_value("b".into(), Value::from_float(b.into()))?;
    context.set_value("c".into(), Value::from_float(c.into()))?;
    context.set_value("alpha".into(), Value::from_float(alpha.into()))?;
    context.set_value("beta".into(), Value::from_float(beta.into()))?;
    context.set_value("gamma".into(), Value::from_float(gamma.into()))?;

    for (var, expr) in path_info.kparam {
        let res = eval_float_with_context(expr.0, &context)?;
        context.set_value(var.0.into(), Value::from_float(res))?;
    }

    let mut points = vec![];
    for (p, (x, y, z)) in path_info.points {
        let eval_x = eval_float_with_context(x.0, &context)?;
        let eval_y = eval_float_with_context(y.0, &context)?;
        let eval_z = eval_float_with_context(z.0, &context)?;

        points.push((p, (FracCoord(eval_x), FracCoord(eval_y), FracCoord(eval_z))));
    }

    let kpath_eval = KpathEval {
        path: path_info.path,
        points,
    };
    Ok(kpath_eval)
}

// --------------------------------- start --------------------------------------

// ---- aP2 ----
static PATH_aP2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    Y GAMMA,
    GAMMA Z,
    R GAMMA,
    GAMMA T,
    U GAMMA,
    GAMMA V,
];

static POINTS_aP2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Z: 0, 0, 1/2,
    Y: 0, 1/2, 0,
    X: 1/2, 0, 0,
    V: 1/2, 1/2, 0,
    U: 1/2, 0, 1/2,
    T: 0, 1/2, 1/2,
    R: 1/2, 1/2, 1/2,
];

static KPARAM_aP2: &[(FracCoordExpr, FracCoordExpr)] = kparams![];

pub static INFO_aP2: KpathInfo = KpathInfo {
    path: PATH_aP2,
    points: POINTS_aP2,
    kparam: KPARAM_aP2,
};

// ---- aP3 ----
static PATH_aP3: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    Y GAMMA,
    GAMMA Z,
    R_2 GAMMA,
    GAMMA T_2,
    U_2 GAMMA,
    GAMMA V_2,
];

static POINTS_aP3: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Z: 0, 0, 1/2,
    Y: 0, 1/2, 0,
    Y_2: 0, -1/2, 0,
    X: 1/2, 0, 0,
    V_2: 1/2, -1/2, 0,
    U_2: -1/2, 0, 1/2,
    T_2: 0, -1/2, 1/2,
    R_2: -1/2, -1/2, 1/2,
];

static KPARAM_aP3: &[(FracCoordExpr, FracCoordExpr)] = kparams![];

pub static INFO_aP3: KpathInfo = KpathInfo {
    path: PATH_aP3,
    points: POINTS_aP3,
    kparam: KPARAM_aP3,
};

// ---- mP1 ----
static PATH_mP1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA Z,
    Z D,
    D B,
    B GAMMA,
    GAMMA A,
    A E,
    E Z,
    Z C_2,
    C_2 Y_2,
    Y_2 GAMMA,
];

static POINTS_mP1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Z: 0, 1/2, 0,
    B: 0, 0, 1/2,
    B_2: 0, 0, -1/2,
    Y: 1/2, 0, 0,
    Y_2: -1/2, 0, 0,
    C: 1/2, 1/2, 0,
    C_2: -1/2, 1/2, 0,
    D: 0, 1/2, 1/2,
    D_2: 0, 1/2, -1/2,
    A: -1/2, 0, 1/2,
    E: -1/2, 1/2, 1/2,
    H: -Y, 0, 1-N,
    H_2: -1+Y, 0, N,
    H_4: -Y, 0, -N,
    M: -Y, 1/2, 1-N,
    M_2: -1+Y, 1/2, N,
    M_4: -Y, 1/2, -N,
];

static KPARAM_mP1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    Y: (1 + a / c * cosbeta) / 2 / sinbeta / sinbeta,
    N: 1/2 + Y * c * cosbeta / a,
);

pub static INFO_mP1: KpathInfo = KpathInfo {
    path: PATH_mP1,
    points: POINTS_mP1,
    kparam: KPARAM_mP1,
};

// ---- mC1 ----
static PATH_mC1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA C,
    C_2 Y_2,
    Y_2 GAMMA,
    GAMMA M_2,
    M_2 D,
    D_2 A,
    A GAMMA,
    L_2 GAMMA,
    GAMMA V_2,
];

static POINTS_mC1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Y_2: -1/2, 1/2, 0,
    Y_4: 1/2, -1/2, 0,
    A: 0, 0, 1/2,
    M_2: -1/2, 1/2, 1/2,
    V: 1/2, 0, 0,
    V_2: 0, 1/2, 0,
    L_2: 0, 1/2, 1/2,
    C: 1-S, 1-S, 0,
    C_2: -1+S, S, 0,
    C_4: S, -1+S, 0,
    D: -1+P, P, 1/2,
    D_2: 1-P, 1-P, 1/2,
    E: -1+Z, 1-Z, 1-H,
    E_2: -Z, Z, H,
    E_4: Z, -Z, 1-H,
];
static KPARAM_mC1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    Z: (2 + a / c * cosbeta) / 4 / sinbeta / sinbeta,
    H: 1/2 - 2 * Z * c * cosbeta / a,
    S: 3/4 - b * b / 4 / a / a / sinbeta / sinbeta,
    P: S - (3/4 - S) * a * cosbeta / c,
);
pub static INFO_mC1: KpathInfo = KpathInfo {
    path: PATH_mC1,
    points: POINTS_mC1,
    kparam: KPARAM_mC1,
};

// ---- mC2 ----
static PATH_mC2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA Y,
    Y M,
    M A,
    A GAMMA,
    L_2 GAMMA,
    GAMMA V_2,
];

static POINTS_mC2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Y: 1/2, 1/2, 0,
    A: 0, 0, 1/2,
    M: 1/2, 1/2, 1/2,
    V_2: 0, 1/2, 0,
    L_2: 0, 1/2, 1/2,
    F: -1+P, 1-P, 1-S,
    F_2: 1-P, P, S,
    F_4: P, 1-P, 1-S,
    H: -Z, Z, X,
    H_2: Z, 1-Z, 1-X,
    H_4: Z, -Z, 1-X,
    G: -M, M, D,
    G_2: M, 1-M, -D,
    G_4: M, -M, -D,
    G_6: 1-M, M, D,
];

static KPARAM_mC2: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    Z: (a * a / b / b + (1 + a / c * cosbeta) / sinbeta / sinbeta) / 4,
    M: (1 + a * a / b / b) / 4,
    D: -a * c * cosbeta / 2 / b / b,
    X: 1/2 - 2 * Z * c * cosbeta / a,
    P: 1 + Z - 2 * M,
    S: X - 2 * D,
);

pub static INFO_mC2: KpathInfo = KpathInfo {
    path: PATH_mC2,
    points: POINTS_mC2,
    kparam: KPARAM_mC2,
};

// ---- mC3 ----

static PATH_mC3: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA A,
    A I_2,
    I M_2,
    M_2 GAMMA,
    GAMMA Y,
    L_2 GAMMA,
    GAMMA V_2,
];
static POINTS_mC3: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Y: 1/2, 1/2, 0,
    A: 0, 0, 1/2,
    M_2: -1/2, 1/2, 1/2,
    V: 1/2, 0, 0,
    V_2: 0, 1/2, 0,
    L_2: 0, 1/2, 1/2,
    I: -1+R, R, 1/2,
    I_2: 1-R, 1-R, 1/2,
    K: -U, U, W,
    K_2: -1+U, 1-U, 1-W,
    K_4: 1-U, U, W,
    H: -Z, Z, E,
    H_2: Z, 1-Z, 1-E,
    H_4: Z, -Z, 1-E,
    N: -F, F, D,
    N_2: F, 1-F, -D,
    N_4: F, -F, -D,
    N_6: 1-F, F, D,
];

static KPARAM_mC3: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    Z: (a * a / b / b + (1 + a / c * cosbeta) / sinbeta / sinbeta) / 4,
    R: 1 - Z * b * b / a / a,
    E: 1/2 - 2 * Z * c * cosbeta / a,
    F: E / 2 + a * a / 4 / b / b + a * c * cosbeta / 2 / b / b,
    U: 2 * F - Z,
    W: c / 2 / a / cosbeta * (1 - 4 * U + a * a * sinbeta * sinbeta / b / b),
    D: -1/4 + W / 2 - Z * c * cosbeta / a,
);
pub static INFO_mC3: KpathInfo = KpathInfo {
    path: PATH_mC3,
    points: POINTS_mC3,
    kparam: KPARAM_mC3,
};

// ---- oP1 ----
static PATH_oP1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X S,
    S Y,
    Y GAMMA,
    GAMMA Z,
    Z U,
    U R,
    R T,
    T Z,
    X U,
    Y T,
    S R,
];
static POINTS_oP1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    X: 1/2, 0, 0,
    Z: 0, 0, 1/2,
    U: 1/2, 0, 1/2,
    Y: 0, 1/2, 0,
    S: 1/2, 1/2, 0,
    T: 0, 1/2, 1/2,
    R: 1/2, 1/2, 1/2,
];
static KPARAM_oP1: &[(FracCoordExpr, FracCoordExpr)] = kparams![];
pub static INFO_oP1: KpathInfo = KpathInfo {
    path: PATH_oP1,
    points: POINTS_oP1,
    kparam: KPARAM_oP1,
};

// ---- oA1 ----

static PATH_oA1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA Y,
    Y C_0,
    SIGMA_0 GAMMA,
    GAMMA Z,
    Z A_0,
    E_0 T,
    T Y,
    GAMMA S,
    S R,
    R Z,
    Z T,
];
static POINTS_oA1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Y: -1/2, 1/2, 0,
    T: -1/2, 1/2, 1/2,
    Z: 0, 0, 1/2,
    S: 0, 1/2, 0,
    R: 0, 1/2, 1/2,
    SIGMA_0: X, X, 0,
    C_0: -X, 1-X, 0,
    A_0: X, X, 1/2,
    E_0: -X, 1-X, 1/2,
];
static KPARAM_oA1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    X: (1 + b * b / c / c) / 4,
);

pub static INFO_oA1: KpathInfo = KpathInfo {
    path: PATH_oA1,
    points: POINTS_oA1,
    kparam: KPARAM_oA1,
};

// ---- oA2 ----
static PATH_oA2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA Y,
    Y F_0,
    DELTA_0 GAMMA,
    GAMMA Z,
    Z B_0,
    G_0 T,
    T Y,
    GAMMA S,
    S R,
    R Z,
    Z T,
];
static POINTS_oA2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Y: 1/2, 1/2, 0,
    T: 1/2, 1/2, 1/2,
    T_2: 1/2, 1/2, -1/2,
    Z: 0, 0, 1/2,
    Z_2: 0, 0, -1/2,
    S: 0, 1/2, 0,
    R: 0, 1/2, 1/2,
    R_2: 0, 1/2, -1/2,
    DELTA_0: -X, X, 0,
    F_0: X, 1-X, 0,
    B_0: -X, X, 1/2,
    B_2: -X, X, -1/2,
    G_0: X, 1-X, 1/2,
    G_2: X, 1-X, -1/2,
];
static KPARAM_oA2: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    X: (1 + c * c / b / b) / 4,
);
pub static INFO_oA2: KpathInfo = KpathInfo {
    path: PATH_oA2,
    points: POINTS_oA2,
    kparam: KPARAM_oA2,
};

// ---- oC1 ----
static PATH_oC1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA Y,
    Y C_0,
    SIGMA_0 GAMMA,
    GAMMA Z,
    Z A_0,
    E_0 T,
    T Y,
    GAMMA S,
    S R,
    R Z,
    Z T,
];
static POINTS_oC1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Y: -1/2, 1/2, 0,
    T: -1/2, 1/2, 1/2,
    Z: 0, 0, 1/2,
    S: 0, 1/2, 0,
    R: 0, 1/2, 1/2,
    SIGMA_0: X, X, 0,
    C_0: -X, 1-X, 0,
    A_0: X, X, 1/2,
    E_0: -X, 1-X, 1/2,
];
static KPARAM_oC1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    X: (1 + a * a / b / b) / 4,
);
pub static INFO_oC1: KpathInfo = KpathInfo {
    path: PATH_oC1,
    points: POINTS_oC1,
    kparam: KPARAM_oC1,
};

// ---- oC2 ----
static PATH_oC2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA Y,
    Y F_0,
    DELTA_0 GAMMA,
    GAMMA Z,
    Z B_0,
    G_0 T,
    T Y,
    GAMMA S,
    S R,
    R Z,
    Z T,
];
static POINTS_oC2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Y: 1/2, 1/2, 0,
    T: 1/2, 1/2, 1/2,
    T_2: 1/2, 1/2, -1/2,
    Z: 0, 0, 1/2,
    Z_2: 0, 0, -1/2,
    S: 0, 1/2, 0,
    R: 0, 1/2, 1/2,
    R_2: 0, 1/2, -1/2,
    DELTA_0: -X, X, 0,
    F_0: X, 1-X, 0,
    B_0: -X, X, 1/2,
    B_2: -X, X, -1/2,
    G_0: X, 1-X, 1/2,
    G_2: X, 1-X, -1/2,
];
static KPARAM_oC2: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    X: (1 + b * b / a / a) / 4,
);
pub static INFO_oC2: KpathInfo = KpathInfo {
    path: PATH_oC2,
    points: POINTS_oC2,
    kparam: KPARAM_oC2,
};

// ---- oF1 ----
static PATH_oF1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA Y,
    Y T,
    T Z,
    Z GAMMA,
    GAMMA SIGMA_0,
    U_0 T,
    Y C_0,
    A_0 Z,
    GAMMA L,
];
static POINTS_oF1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    T: 1, 1/2, 1/2,
    Z: 1/2, 1/2, 0,
    Y: 1/2, 0, 1/2,
    SIGMA_0: 0, H, H,
    U_0: 1, 1-H, 1-H,
    A_0: 1/2, 1/2+J, J,
    C_0: 1/2, 1/2-J, 1-J,
    L: 1/2, 1/2, 1/2,
];
static KPARAM_oF1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    J: (1 + a * a / b / b - a * a / c / c) / 4,
    H: (1 + a * a / b / b + a * a / c / c) / 4,
);
pub static INFO_oF1: KpathInfo = KpathInfo {
    path: PATH_oF1,
    points: POINTS_oF1,
    kparam: KPARAM_oF1,
};

// ---- oF2 ----
static PATH_oF2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA T,
    T Z,
    Z Y,
    Y GAMMA,
    GAMMA LAMBDA_0,
    Q_0 Z,
    T G_0,
    H_0 Y,
    GAMMA L,
];
static POINTS_oF2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    T: 0, 1/2, 1/2,
    Z: 1/2, 1/2, 1,
    Y: 1/2, 0, 1/2,
    LAMBDA_0: K, K, 0,
    Q_0: 1-K, 1-K, 1,
    G_0: 1/2-J, 1-J, 1/2,
    H_0: 1/2+J, J, 1/2,
    L: 1/2, 1/2, 1/2,
];
static KPARAM_oF2: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    J: (1 + c * c / a / a - c * c / b / b) / 4,
    K: (1 + c * c / a / a + c * c / b / b) / 4,
);
pub static INFO_oF2: KpathInfo = KpathInfo {
    path: PATH_oF2,
    points: POINTS_oF2,
    kparam: KPARAM_oF2,
};

// ---- oF3 ----
static PATH_oF3: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA Y,
    Y C_0,
    A_0 Z,
    Z B_0,
    D_0 T,
    T G_0,
    H_0 Y,
    T GAMMA,
    GAMMA Z,
    GAMMA L,
];
static POINTS_oF3: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    T: 0, 1/2, 1/2,
    Z: 1/2, 1/2, 0,
    Y: 1/2, 0, 1/2,
    A_0: 1/2, 1/2+H, H,
    C_0: 1/2, 1/2-H, 1-H,
    B_0: 1/2+K, 1/2, K,
    D_0: 1/2-K, 1/2, 1-K,
    G_0: P, 1/2+P, 1/2,
    H_0: 1-P, 1/2-P, 1/2,
    L: 1/2, 1/2, 1/2,
];
static KPARAM_oF3: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    H: (1 + a * a / b / b - a * a / c / c) / 4,
    K: (1 + b * b / a / a - b * b / c / c) / 4,
    P: (1 + c * c / b / b - c * c / a / a) / 4,
);
pub static INFO_oF3: KpathInfo = KpathInfo {
    path: PATH_oF3,
    points: POINTS_oF3,
    kparam: KPARAM_oF3,
};

// ---- oI1 ----
static PATH_oI1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X F_2,
    SIGMA_0 GAMMA,
    GAMMA Y_0,
    U_0 X,
    GAMMA R,
    R W,
    W S,
    S GAMMA,
    GAMMA T,
    T W,
];
static POINTS_oI1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    X: 1/2, 1/2, -1/2,
    S: 1/2, 0, 0,
    R: 0, 1/2, 0,
    T: 0, 0, 1/2,
    W: 1/4, 1/4, 1/4,
    SIGMA_0: -Z, Z, Z,
    F_2: Z, 1-Z, -Z,
    Y_0: H, -H, H,
    U_0: 1-H, H, -H,
    L_0: -N, N, 1/2-D,
    M_0: N, -N, 1/2+D,
    J_0: 1/2-D, 1/2+D, -N,
];
static KPARAM_oI1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    Z: (1 + a * a / c / c) / 4,
    H: (1 + b * b / c / c) / 4,
    D: (b * b - a * a) / 4 / c / c,
    N: (a * a + b * b) / 4 / c / c,
);
pub static INFO_oI1: KpathInfo = KpathInfo {
    path: PATH_oI1,
    points: POINTS_oI1,
    kparam: KPARAM_oI1,
};

// ---- oI2 ----
static PATH_oI2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X U_2,
    Y_0 GAMMA,
    GAMMA LAMBDA_0,
    G_2 X,
    GAMMA R,
    R W,
    W S,
    S GAMMA,
    GAMMA T,
    T W,
];
static POINTS_oI2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    X: -1/2, 1/2, 1/2,
    S: 1/2, 0, 0,
    R: 0, 1/2, 0,
    T: 0, 0, 1/2,
    W: 1/4, 1/4, 1/4,
    Y_0: Z, -Z, Z,
    U_2: -Z, Z, 1-Z,
    LAMBDA_0: H, H, -H,
    G_2: -H, 1-H, H,
    K: 1/2-D, -N, N,
    K_2: 1/2+D, N, -N,
    K_4: -N, 1/2-D, 1/2+D,
];
static KPARAM_oI2: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    Z: (1 + b * b / a / a) / 4,
    H: (1 + c * c / a / a) / 4,
    D: (c * c - b * b) / 4 / a / a,
    N: (b * b + c * c) / 4 / a / a,
);
pub static INFO_oI2: KpathInfo = KpathInfo {
    path: PATH_oI2,
    points: POINTS_oI2,
    kparam: KPARAM_oI2,
};

// ---- oI3 ----
static PATH_oI3: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X F_0,
    SIGMA_0 GAMMA,
    GAMMA LAMBDA_0,
    G_0 X,
    GAMMA R,
    R W,
    W S,
    S GAMMA,
    GAMMA T,
    T W,
];
static POINTS_oI3: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    X: 1/2, -1/2, 1/2,
    S: 1/2, 0, 0,
    R: 0, 1/2, 0,
    T: 0, 0, 1/2,
    W: 1/4, 1/4, 1/4,
    SIGMA_0: -Y, Y, Y,
    F_0: Y, -Y, 1-Y,
    LAMBDA_0: Z, Z, -Z,
    G_0: 1-Z, -Z, Z,
    V_0: M, 1/2-D, -M,
    H_0: -M, 1/2+D, M,
    H_2: 1/2+D, -M, 1/2-D,
];
static KPARAM_oI3: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    Z: (1 + c * c / b / b) / 4,
    Y: (1 + a * a / b / b) / 4,
    D: (a * a - c * c) / 4 / b / b,
    M: (c * c + a * a) / 4 / b / b,
);
pub static INFO_oI3: KpathInfo = KpathInfo {
    path: PATH_oI3,
    points: POINTS_oI3,
    kparam: KPARAM_oI3,
};

// ---- tP1 ----
static PATH_tP1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X M,
    M GAMMA,
    GAMMA Z,
    Z R,
    R A,
    A Z,
    X R,
    M A,
];
static POINTS_tP1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    Z: 0, 0, 1/2,
    M: 1/2, 1/2, 0,
    A: 1/2, 1/2, 1/2,
    R: 0, 1/2, 1/2,
    X: 0, 1/2, 0,
];
static KPARAM_tP1: &[(FracCoordExpr, FracCoordExpr)] = kparams![];
pub static INFO_tP1: KpathInfo = KpathInfo {
    path: PATH_tP1,
    points: POINTS_tP1,
    kparam: KPARAM_tP1,
};

// ---- tI1 ----
static PATH_tI1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X M,
    M GAMMA,
    GAMMA Z,
    Z_0 M,
    X P,
    P N,
    N GAMMA,
];
static POINTS_tI1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    M: -1/2, 1/2, 1/2,
    X: 0, 0, 1/2,
    P: 1/4, 1/4, 1/4,
    Z: H, H, -H,
    Z_0: -H, 1-H, H,
    N: 0, 1/2, 0,
];
static KPARAM_tI1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    H: (1 + c * c / a / a) / 4,
);
pub static INFO_tI1: KpathInfo = KpathInfo {
    path: PATH_tI1,
    points: POINTS_tI1,
    kparam: KPARAM_tI1,
};

// ---- tI2 ----
static PATH_tI2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X P,
    P N,
    N GAMMA,
    GAMMA M,
    M S,
    S_0 GAMMA,
    X R,
    G M,
];
static POINTS_tI2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    M: 1/2, 1/2, -1/2,
    X: 0, 0, 1/2,
    P: 1/4, 1/4, 1/4,
    N: 0, 1/2, 0,
    S_0: -H, H, H,
    S: H, 1-H, -H,
    R: -Z, Z, 1/2,
    G: 1/2, 1/2, -Z,
];
static KPARAM_tI2: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    H: (1 + a * a / c / c) / 4,
    Z: a * a / 2 / c / c,
);
pub static INFO_tI2: KpathInfo = KpathInfo {
    path: PATH_tI2,
    points: POINTS_tI2,
    kparam: KPARAM_tI2,
};

// ---- hR1 ----
static PATH_hR1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA T,
    T H_2,
    H_0 L,
    L GAMMA,
    GAMMA S_0,
    S_2 F,
    F GAMMA,
];

static POINTS_hR1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    T: 1/2, 1/2, 1/2,
    L: 1/2, 0, 0,
    L_2: 0, -1/2, 0,
    L_4: 0, 0, -1/2,
    F: 1/2, 0, 1/2,
    F_2: 1/2, 1/2, 0,
    S_0: N, -N, 0,
    S_2: 1-N, 0, N,
    S_4: N, 0, -N,
    S_6: 1-N, N, 0,
    H_0: 1/2, -1+Y, 1-Y,
    H_2: Y, 1-Y, 1/2,
    H_4: Y, 1/2, 1-Y,
    H_6: 1/2, 1-Y, -1+Y,
    M_0: N, -1+Y, N,
    M_2: 1-N, 1-Y, 1-N,
    M_4: Y, N, N,
    M_6: 1-N, 1-N, 1-Y,
    M_8: N, N, -1+Y,
];

static KPARAM_hR1: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    D: a * a / 4 / c / c,
    Y: 5 / 6 - 2 *D,
    N: 1 / 3 + D,
);

pub static INFO_hR1: KpathInfo = KpathInfo {
    path: PATH_hR1,
    points: POINTS_hR1,
    kparam: KPARAM_hR1,
};

// ---- hR2 ----
static PATH_hR2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA L,
    L T,
    T P_0,
    P_2 GAMMA,
    GAMMA F,
];
static POINTS_hR2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    T: 1/2, -1/2, 1/2,
    P_0: H, -1+H, H,
    P_2: H, H, H,
    R_0: 1-H, -H, -H,
    M: 1-N, -N, 1-N,
    M_2: N, -1+N, -1+N,
    L: 1/2, 0, 0,
    F: 1/2, -1/2, 0,
];
static KPARAM_hR2: &[(FracCoordExpr, FracCoordExpr)] = kparams!(
    Z: 1/6 - c * c / 9 / a / a,
    H: 1/2 - 2 * Z,
    N: 1/2 + Z,
);

pub static INFO_hR2: KpathInfo = KpathInfo {
    path: PATH_hR2,
    points: POINTS_hR2,
    kparam: KPARAM_hR2,
};

// ---- hP1 ----
static PATH_hP1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA M,
    M K,
    K GAMMA,
    GAMMA A,
    A L,
    L H,
    H A,
    L M,
    H K,
    K H_2,
];
static POINTS_hP1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    A: 0, 0, 1/2,
    K: 1/3, 1/3, 0,
    H: 1/3, 1/3, 1/2,
    H_2: 1/3, 1/3, -1/2,
    M: 1/2, 0, 0,
    L: 1/2, 0, 1/2,
];
static KPARAM_hP1: &[(FracCoordExpr, FracCoordExpr)] = kparams![];
pub static INFO_hP1: KpathInfo = KpathInfo {
    path: PATH_hP1,
    points: POINTS_hP1,
    kparam: KPARAM_hP1,
};

// ---- hP2 ----
static PATH_hP2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA M,
    M K,
    K GAMMA,
    GAMMA A,
    A L,
    L H,
    H A,
    L M,
    H K,
];
static POINTS_hP2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    A: 0, 0, 1/2,
    K: 1/3, 1/3, 0,
    H: 1/3, 1/3, 1/2,
    H_2: 1/3, 1/3, -1/2,
    M: 1/2, 0, 0,
    L: 1/2, 0, 1/2,
];
static KPARAM_hP2: &[(FracCoordExpr, FracCoordExpr)] = kparams![];
pub static INFO_hP2: KpathInfo = KpathInfo {
    path: PATH_hP2,
    points: POINTS_hP2,
    kparam: KPARAM_hP2,
};

// ---- cP1 ----
static PATH_cP1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X M,
    M GAMMA,
    GAMMA R,
    R X,
    R M,
    M X_1,
];
static POINTS_cP1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    R: 1/2, 1/2, 1/2,
    M: 1/2, 1/2, 0,
    X: 0, 1/2, 0,
    X_1: 1/2, 0, 0,
];
static KPARAM_cP1: &[(FracCoordExpr, FracCoordExpr)] = kparams![];
pub static INFO_cP1: KpathInfo = KpathInfo {
    path: PATH_cP1,
    points: POINTS_cP1,
    kparam: KPARAM_cP1,
};

// ---- cP2 ----
static PATH_cP2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X M,
    M GAMMA,
    GAMMA R,
    R X,
    R M,
];
static POINTS_cP2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    R: 1/2, 1/2, 1/2,
    M: 1/2, 1/2, 0,
    X: 0, 1/2, 0,
    X_1: 1/2, 0, 0,
];
static KPARAM_cP2: &[(FracCoordExpr, FracCoordExpr)] = kparams![];
pub static INFO_cP2: KpathInfo = KpathInfo {
    path: PATH_cP2,
    points: POINTS_cP2,
    kparam: KPARAM_cP2,
};

// ---- cF1 ----
static PATH_cF1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X U,
    K GAMMA,
    GAMMA L,
    L W,
    W X,
    X W_2,
];
static POINTS_cF1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    X: 1/2, 0, 1/2,
    L: 1/2, 1/2, 1/2,
    W: 1/2, 1/4, 3/4,
    W_2: 3/4, 1/4, 1/2,
    K: 3/8, 3/8, 3/4,
    U: 5/8, 1/4, 5/8,
];
static KPARAM_cF1: &[(FracCoordExpr, FracCoordExpr)] = kparams![];
pub static INFO_cF1: KpathInfo = KpathInfo {
    path: PATH_cF1,
    points: POINTS_cF1,
    kparam: KPARAM_cF1,
};

// ---- cF2 ----
static PATH_cF2: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA X,
    X U,
    K GAMMA,
    GAMMA L,
    L W,
    W X,
];
static POINTS_cF2: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    X: 1/2, 0, 1/2,
    L: 1/2, 1/2, 1/2,
    W: 1/2, 1/4, 3/4,
    W_2: 3/4, 1/4, 1/2,
    K: 3/8, 3/8, 3/4,
    U: 5/8, 1/4, 5/8,
];
static KPARAM_cF2: &[(FracCoordExpr, FracCoordExpr)] = kparams![];
pub static INFO_cF2: KpathInfo = KpathInfo {
    path: PATH_cF2,
    points: POINTS_cF2,
    kparam: KPARAM_cF2,
};

// ---- cI1 ----
static PATH_cI1: &[(HighSymmetryPoint, HighSymmetryPoint)] = kpath![
    GAMMA H,
    H N,
    N GAMMA,
    GAMMA P,
    P H,
    P N,
];
static POINTS_cI1: &[(
    HighSymmetryPoint,
    (FracCoordExpr, FracCoordExpr, FracCoordExpr),
)] = kpoints![
    GAMMA: 0, 0, 0,
    H: 1/2, -1/2, 1/2,
    P: 1/4, 1/4, 1/4,
    N: 0, 0, 1/2,
];

static KPARAM_cI1: &[(FracCoordExpr, FracCoordExpr)] = kparams![];

pub static INFO_cI1: KpathInfo = KpathInfo {
    path: PATH_cI1,
    points: POINTS_cI1,
    kparam: KPARAM_cI1,
};

// ----------- end -----------

#[allow(
    clippy::unreadable_literal,
    clippy::excessive_precision,
    non_snake_case
)]
#[cfg(test)]
mod tests {
    use ccmat_core::{atomic_number, lattice_angstrom, sites_frac_coord, CrystalBuilder};

    use super::*;

    #[test]
    fn eval_aP2() {
        // inv
        let lattice = lattice_angstrom![
            (-1.2489688500000000, 2.4101964100000002, -1.1084066800000001),
            (-3.8805945199999998, -0.4578071300000000, 0.1425430700000000),
            (-1.4313846100000001, 2.5585931100000003, 3.3956427200000006),
        ];
        let sites = sites_frac_coord![
            (0.5000000000000000, 0.5000000000000000, 0.5000000000000000), atomic_number!(Cu);
            (0.0000000000000000, 0.0000000000000000, 0.0000000000000000), atomic_number!(Ag);
            (0.1691864300000000, 0.3370706700000000, 0.2835515499999998), atomic_number!(O);
            (0.8308135699999999, 0.6629293300000000, 0.7164484500000001), atomic_number!(O);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let kpatheval = eval(lookup(&ExtBravaisClass::aP2), s.lattice().lattice_params()).unwrap();
        let point_R = kpatheval.points[7];

        assert_eq!(point_R.0, &HighSymmetryPoint::R);
        assert_eq!(point_R.1, (FracCoord(0.5), FracCoord(0.5), FracCoord(0.5)));
    }

    #[test]
    fn eval_hR1() {
        // inv
        let lattice = lattice_angstrom![
            (3.7622432546763140, 0.0000000000000000, 0.0000000000000000),
            (-1.8811216273381570, 3.2581982337663353, 0.0000000000000000),
            (0.0000000000000000, 0.0000000000000000, 13.5138670706735571),
        ];
        let sites = sites_frac_coord![
            (0.6666666666666667, 0.3333333333333333, 0.8333333333333333), atomic_number!(In);
            (0.3333333333333333, 0.6666666666666666, 0.1666666666666666), atomic_number!(In);
            (0.0000000000000000, 0.0000000000000000, 0.4999999999999999), atomic_number!(In);
            (0.0000000000000000, 0.0000000000000000, 0.0000000000000000), atomic_number!(Hg);
            (0.6666666666666666, 0.3333333333333333, 0.3333333333333333), atomic_number!(Hg);
            (0.3333333333333333, 0.6666666666666666, 0.6666666666666666), atomic_number!(Hg);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let _ = eval(lookup(&ExtBravaisClass::hR1), s.lattice().lattice_params()).unwrap();

        // noinv
        let lattice = lattice_angstrom![
            (4.3089558515497082, 0.0000000000000000, 0.0000000000000000),
            (-2.1544779257748541, 3.7316652312276553, 0.0000000000000000),
            (0.0000000000000000, 0.0000000000000000, 10.2486516858249015),
        ];
        let sites = sites_frac_coord![
            (0.0000000000000000, 0.0000000000000000, 0.2797794000000000), atomic_number!(Cu);
            (0.6666666666666666, 0.3333333333333333, 0.6131127333333333), atomic_number!(Cu);
            (0.3333333333333333, 0.6666666666666666, 0.9464460666666666), atomic_number!(Cu);
            (0.0000000000000000, 0.0000000000000000, 0.0252206000000002), atomic_number!(I);
            (0.6666666666666666, 0.3333333333333333, 0.3585539333333334), atomic_number!(I);
            (0.3333333333333333, 0.6666666666666666, 0.6918872666666668), atomic_number!(I);
        ];

        let s = CrystalBuilder::new()
            .with_lattice(&lattice)
            .with_frac_sites(sites)
            .build()
            .unwrap();

        let _ = eval(lookup(&ExtBravaisClass::hR1), s.lattice().lattice_params()).unwrap();
    }
}
