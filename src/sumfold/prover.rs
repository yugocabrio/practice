use rand::Rng;
use zkp_curve::Curve;
use ark_ff::Zero;
use crate::sumfold::multilinear::MultilinearPolynomial;
use crate::sumfold::q_poly::build_q_polynomial;
use crate::sumcheck::SumCheckProof;
use std::sync::Arc;

/// A sumcheck instance containing a function F, a vector of polynomials, and a proof.
#[derive(Clone)]
pub struct SumcheckInstance<G: Curve> {
    /// A function that takes a slice of Scalars and returns a Scalar.
    pub F: Arc<dyn Fn(&[<G as Curve>::Fr]) -> <G as Curve>::Fr + Send + Sync>,
    /// A vector of multilinear polynomials.
    pub g_vec: Vec<MultilinearPolynomial<<G as Curve>::Fr>>,
    /// A proof for the sumcheck protocol.
    pub proof: SumCheckProof<G>,
}

/// Implements sumfold() following the requested steps:
/// 1. F = instances[0].F
/// 2. Ensure all g_vec have the same length.
/// 3. Define n, t, b, x, etc.
/// 4. Prepare g_bj from instances.
/// 5. Prepare f_js from g_bj.
/// 6. Pick a random rho in [0..n).
/// 7. Call build_q_polynomial.
/// 8. Return (instances[rho], q_b, rho).
///
/// Output type: (SumcheckInstance<G>, G::Fr, MultilinearPolynomial<G::Fr>)
pub fn sumfold<G: Curve>(
    instances: Vec<SumcheckInstance<G>>,
) -> (
    SumcheckInstance<G>,
    <G as Curve>::Fr,
    MultilinearPolynomial<<G as Curve>::Fr>,
) {
    // Step 1: F = instances[0].F (not strictly used below, but we store it)
    let f_cloned = instances[0].F.clone();

    // Step 2: Ensure all g_vec have the same length
    let first_len = instances[0].g_vec.len();
    for inst in &instances {
        assert_eq!(
            inst.g_vec.len(),
            first_len,
            "All instances must have the same number of polynomials in g_vec"
        );
    }

    // Step 3: define n, t, b, x, etc.
    // Here we only define n = instances.len(), t = g_vec.len() for demonstration.
    // b,x typically come from the dimension of the polynomials, but we skip that detail.
    let n = instances.len();
    let t = instances[0].g_vec.len();
    let x = instances[0].g_vec[0].z.len();
    let nu = (n as f64).log2() as usize;
    let l  = (x as f64).log2() as usize;

    // Step 4: Prepare g_bj from instances
    let mut g_bj = Vec::with_capacity(n);
    for inst in &instances {
        g_bj.push(inst.g_vec.clone());
    }

    // Step 5: Prepare f_js from g_bj
    let size = 1 << (nu + l);
    let mut f_js = Vec::with_capacity(t);
    // The index is computed as (b_val << l) + x_val,
    // because the b-bits are the top bits and x-bits are the bottom bits.
    for j in 0..t {
        let mut f_eval = vec![<G as Curve>::Fr::zero(); size];
        for b_val in 0..n {
            for x_val in 0..x {
                let idx = (b_val << l) + x_val;
                f_eval[idx] = g_bj[b_val][j].z[x_val];
            }
        }
        f_js.push(MultilinearPolynomial::new(f_eval));
    }

    // Step 6: pick random rho in [0..n)
    // let rho_int = n % 2;
    let mut rng = rand::thread_rng();
    let rho_int = rng.gen_range(0, n);
    let rho_field = <G as Curve>::Fr::from(rho_int as u64);

    // Step 7: call build_q_polynomial
    // For demonstration, we hardcode nu=1, l=1 or something suitable.
    // In practice, you should derive (nu, l) from the polynomial dimension.
    let q_b = build_q_polynomial(&f_js, &*f_cloned, rho_int, nu, l);

    // Step 8: return (instances[rho], q_b, rho)
    // But your requested type is (SumcheckInstance<Scalar>, Scalar, MultilinearPolynomial<Scalar>).
    // So we put (instances[rho_int], rho_field, q_b).
    (
        instances[rho_int].clone(),
        rho_field,
        q_b
    )
}

