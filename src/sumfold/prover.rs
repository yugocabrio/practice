use ark_ff::Zero;
use rand::Rng;
use zkp_curve::Curve;
use std::sync::Arc;

use crate::sumfold::multilinear::MultilinearPolynomial;
use crate::sumfold::q_poly::build_Q_polynomial;

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
/// 4. Prepare g_bj from instances.
/// 5. Prepare f_js from g_bj.
/// 6. Pick a random rho in [0..n).
/// 7. Call build_q_polynomial.
/// 8. Return (chosen_instance, rho_field, q_b).
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

    // Step 4: Prepare g_bj from instances
    let mut g_bj = Vec::with_capacity(n);
    for inst in &instances {
        g_bj.push(inst.g_vec.clone());
    }

    // Step 5: Prepare f_js from g_bj
    let size = 1 << (nu + l);
    let mut f_js = Vec::with_capacity(t);
    // index = (b_val << l) + x_val
    for j in 0..t {
        let mut f_eval = vec![<G as Curve>::Fr::zero(); size];
        for b_val in 0..n {
            for x_val in 0..x {
                let idx = (b_val << l) + x_val;
                if idx < size {
                    f_eval[idx] = g_bj[b_val][j].z[x_val];
                }
            }
        }
        f_js.push(MultilinearPolynomial::new(f_eval));
    }

    // Step 6: pick random rho in [0..n)
    let mut rng = rand::thread_rng();
    let rho_int = rng.gen_range(0, n); // for rand >= 0.7, use two-arg version e.g. rng.gen_range(0..n)
    let rho_field = <G as Curve>::Fr::from(rho_int as u64);

    // Step 7: call build_q_polynomial
    let Q_b = build_Q_polynomial(&f_js, &*f_cloned, rho_int, nu, l);

    // Step 8: return
    (
        instances[rho_int].clone(),
        rho_field,
        Q_b
    )
}
