/*
*
* The **math related** macros defined in the crate are major for internal usage.
* The macros that for users will be re-exported by `ccmat`.
*/

use proc_macro::TokenStream;
use quote::quote;
use syn::parse::{Parse, ParseStream};
use syn::punctuated::Punctuated;
use syn::spanned::Spanned;
use syn::visit_mut::VisitMut;
use syn::{BinOp, ExprBinary, ExprCall, ExprLit, ExprPath, Lit, Path, PathSegment};
use syn::{Expr, Result, Token};

#[rustfmt::skip]
struct MatrixInput {
    mat: [[Expr; 3]; 3],
}

fn detect_integer_division(expr: &Expr) -> Result<()> {
    match expr {
        Expr::Binary(ExprBinary {
            left, op, right, ..
        }) => {
            // If op is `/`, check if both sides look integer-ish
            if matches!(op, BinOp::Div(_)) {
                is_integer_expr(left)?;
                is_integer_expr(right)?;
            }

            // Recurse down left/right anyway
            detect_integer_division(left)?;
            detect_integer_division(right)?;
            Ok(())
        }

        // Parentheses
        Expr::Paren(p) => detect_integer_division(&p.expr),

        // Unary ops: `-1/2`
        Expr::Unary(u) => detect_integer_division(&u.expr),

        // Function call arguments, include Path Call e.g f64::cos(1/2)
        Expr::Call(call) => {
            for arg in &call.args {
                detect_integer_division(arg)?;
            }
            Ok(())
        }

        // Method call args: foo.bar(1/2)
        Expr::MethodCall(call) => {
            for arg in &call.args {
                detect_integer_division(arg)?;
            }
            detect_integer_division(&call.receiver)?;

            Ok(())
        }

        // Other expressions — no division here
        _ => Ok(()),
    }
}

fn is_integer_expr(expr: &Expr) -> Result<()> {
    match expr {
        Expr::Path(ExprPath { path, .. }) => {
            if path.segments.len() > 2
                && path.segments[0].ident == "std"
                && path.segments[1].ident == "f64"
            {
                Ok(())
            } else {
                Err(syn::Error::new(
                    expr.span(),
                    "cannot determine float type at compile time; convert explicitly using `f64::from(...)`"
                ))
            }
        }
        Expr::Lit(ExprLit {
            lit: Lit::Int(_), ..
        }) => {
            Err(syn::Error::new(
                expr.span(),
                "integer division detected; Rust performs integer division here, but this macro expects floats. 
Write `1.0 / 2.0` instead of `1 / 2`."
            ))
        },
        Expr::Paren(p) => is_integer_expr(&p.expr),
        Expr::Unary(u) => is_integer_expr(&u.expr), // handle -1
        _ => Ok(()),
    }
}

struct PathPrefixF64Attach;

impl VisitMut for PathPrefixF64Attach {
    fn visit_expr_mut(&mut self, i: &mut syn::Expr) {
        // PI -> f64::const::PI
        if let Expr::Path(ExprPath { path, .. }) = i {
            if path.segments.len() == 1 {
                let ident = &path.segments[0].ident;
                if ident == "PI" || ident == "Pi" || ident == "pi" {
                    *i = syn::parse_quote!(std::f64::consts::PI);
                    return;
                }
            }
        }

        if let Expr::Call(ExprCall { func, .. }) = i {
            // cos(...) -> f64::cos(...)
            if let Expr::Path(ExprPath { path, .. }) = &mut **func {
                if path.segments.len() == 1 {
                    let ident = &path.segments[0].ident;
                    if ident == "abs"
                        || ident == "cos"
                        || ident == "sin"
                        || ident == "tan"
                        || ident == "acos"
                        || ident == "asin"
                        || ident == "atan"
                        || ident == "acosh"
                        || ident == "asinh"
                        || ident == "atanh"
                        || ident == "sqrt"
                        || ident == "log"
                    {
                        let mut segs = Punctuated::new();
                        segs.push(PathSegment::from(syn::Ident::new("f64", path.span())));
                        segs.push(PathSegment::from(ident.clone()));

                        *path = Path {
                            leading_colon: None,
                            segments: segs,
                        };
                    }
                }
            }
        }

        syn::visit_mut::visit_expr_mut(self, i);
    }
}

impl Parse for MatrixInput {
    fn parse(input: ParseStream) -> Result<Self> {
        let mut rows = Vec::new();
        let mut current_row = Vec::new();
        let span = input.span();

        while !input.is_empty() {
            if let Ok(mut expr) = input.parse::<Expr>() {
                // f64::op
                PathPrefixF64Attach.visit_expr_mut(&mut expr);

                // detect integer devide as an error e.g raise on 1/2
                detect_integer_division(&expr)?;

                if current_row.len() > 2 {
                    return Err(syn::Error::new(
                        expr.span(),
                        "cannot have more than 3 items per row",
                    ));
                }
                current_row.push(expr);
            }

            if input.peek(Token!(;)) {
                input.parse::<Token!(;)>()?;
                if current_row.len() != 3 {
                    return Err(syn::Error::new(
                        current_row[0].span(),
                        "each row must contain exactly 3 elements — this row has fewer; the spacing may be ambiguous. Consider using `,` to separate elements explicitly."
                    ));
                }

                rows.push(std::mem::take(&mut current_row));
            } else if input.peek(Token!(,)) {
                // allow ',' commas optionally
                input.parse::<Token!(,)>()?;
            }
        }

        // after the last
        if !current_row.is_empty() {
            rows.push(current_row);
        }

        if rows.len() != 3 {
            return Err(syn::Error::new(span, "expect 3 rows for a 3x3 matrix"));
        }

        let mat = [
            [rows[0][0].clone(), rows[0][1].clone(), rows[0][2].clone()],
            [rows[1][0].clone(), rows[1][1].clone(), rows[1][2].clone()],
            [rows[2][0].clone(), rows[2][1].clone(), rows[2][2].clone()],
        ];
        Ok(MatrixInput { mat })
    }
}

#[proc_macro]
pub fn matrix_3x3(tokens: TokenStream) -> TokenStream {
    let tokens: proc_macro2::TokenStream = tokens.into();
    let mat_input = syn::parse2::<MatrixInput>(tokens);

    let mat_input = match mat_input {
        Ok(input) => input,
        Err(err) => return err.to_compile_error().into(),
    };

    let MatrixInput { mat } = mat_input;

    let row0 = mat[0].iter().map(|x| quote!(f64::from(#x)));
    let row1 = mat[1].iter().map(|x| quote!(f64::from(#x)));
    let row2 = mat[2].iter().map(|x| quote!(f64::from(#x)));

    let expand = quote! {{
        let mat: [[f64; 3]; 3] = [
            [#(#row0,)*],
            [#(#row1,)*],
            [#(#row2,)*],
        ];
        mat
    }};

    expand.into()
}
