use rand::Rng;
use zkp_curve::Curve;
use std::sync::Arc;

use crate::sumfold::multilinear::MultilinearPolynomial;
use crate::sumfold::q_poly::build_Q_polynomial;

use super::fj_poly::build_fj_polynomial;

#[allow(non_snake_case)]
#[derive(Clone)]
pub struct SumfoldInstance<G: Curve> {
    /// A function that takes a slice of Scalars (e.g. g0(x), g1(x), ...) and returns a single Scalar.
    pub F_func: Arc<dyn Fn(&[<G as Curve>::Fr]) -> <G as Curve>::Fr + Send + Sync>,

    /// A vector of multilinear polynomials: e.g. [g0, g1, ...].
    /// sumcheckのx軸方向に取るtこのgベクトル
    pub g_vec: Vec<MultilinearPolynomial<<G as Curve>::Fr>>,
}

/// Implements sumfold() following the requested steps:
/// 1. F = instances[0].f_func
/// 2. Ensure all g_vec have the same length.
/// 3. Define n, t, x, etc.
/// 4. Prepare f_js from g_bj.
/// 5. Pick random rho in [0..n).
/// 6. Call build_Q_polynomial.
/// 7. Return.
///
/// Output type: (SumfoldInstance<G>, G::Fr, MultilinearPolynomial<G::Fr>)
#[allow(non_snake_case)]
pub fn sumfold<G: Curve>(
    instances: Vec<SumfoldInstance<G>>,
) -> (
    SumfoldInstance<G>,
    <G as Curve>::Fr,
    MultilinearPolynomial<<G as Curve>::Fr>,
) {
    // Step 1: F = instances[0].f_func
    let f_cloned = instances[0].F_func.clone();

    // Step 2: Ensure all g_vec have the same length
    let first_len = instances[0].g_vec.len();
    for inst in &instances {
        assert_eq!(
            inst.g_vec.len(),
            first_len,
            "All instances must have the same number of polynomials in g_vec"
        );
    }

    // Step 3: define n, t, x, etc.
    let n = instances.len();
    let t = instances[0].g_vec.len();
    let x = instances[0].g_vec[0].z.len();

    // Make sure n and x are powers of 2
    let n_pow2 = n.next_power_of_two();
    let x_pow2 = x.next_power_of_two();

    // nu, l for interpolation
    let nu = (n_pow2 as f64).log2() as usize;
    let l  = (x_pow2 as f64).log2() as usize;

    // Step 4: Prepare f_js from g_bj
    let gs: Vec<Vec<_>> = (0..t)
        .map(|j| instances.iter().map(|inst| inst.g_vec[j].clone()).collect())
        .collect();

    assert_eq!(gs.len(), t, "gs must have t elements");

    let f_js: Vec<_> = gs.iter().map(|gs_for_j| build_fj_polynomial(gs_for_j)).collect();


    // Step 5: pick random rho in [0..n)
    let mut rng = rand::thread_rng();
    let rho_int = rng.gen_range(0, n); // for rand >= 0.7, use two-arg version e.g. rng.gen_range(0..n)
    let rho_field = <G as Curve>::Fr::from(rho_int as u64);

    // Step 6: call build_q_polynomial
    let Q_b = build_Q_polynomial(&f_js, &*f_cloned, rho_int, nu, l);

    // Step 7: return
    (
        instances[rho_int].clone(),
        rho_field,
        Q_b
    )
}
