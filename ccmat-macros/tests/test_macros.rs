use ccmat_macros::matrix_3x3;

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
        1 2 3;
        4 5 6.1;
        7 8 9;
    ];

    assert_eq_approx!(mat[0][0], 1.0);
    assert_eq_approx!(mat[1][2], 6.1);

    let mat = matrix_3x3![
        1 + 2 + 1 2 3;  // althrough ugly and not recommemded
        4 5 6.1;
        7 8 9;
    ];

    assert_eq_approx!(mat[0][0], 4.0);
    assert_eq_approx!(mat[1][2], 6.1);

    let mat = matrix_3x3![
        1 + 2 + 1, 2, 3; // better
        4, 5, 6.1;
        7, 8, 9;
    ];

    assert_eq_approx!(mat[0][0], 4.0);
    assert_eq_approx!(mat[1][2], 6.1);

    let mat = matrix_3x3![
        cos(PI/3.), PI, f64::cos(PI/6.);
        4, 5, 6.1;
        7, 8, sqrt(9.0);
    ];

    assert_eq_approx!(mat[0][0], 0.5);
    assert_eq_approx!(mat[1][2], 6.1);
    assert_eq_approx!(mat[2][2], 3.0);
}

#[test]
fn ui() {
    let t = trybuild::TestCases::new();
    t.compile_fail("tests/ui/*.rs");
}
