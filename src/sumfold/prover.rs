use ark_ff::Zero;
use rand::Rng;
use zkp_curve::Curve;
use std::sync::Arc;

use crate::sumfold::multilinear::MultilinearPolynomial;
use crate::sumfold::q_poly::build_q_polynomial;

#[derive(Clone)]
pub struct SumfoldInstance<G: Curve> {
    /// A function that takes a slice of Scalars (e.g. g0(x), g1(x), ...) and returns a single Scalar.
    pub f_func: Arc<dyn Fn(&[<G as Curve>::Fr]) -> <G as Curve>::Fr + Send + Sync>,

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
pub fn sumfold<G: Curve>(
    instances: Vec<SumfoldInstance<G>>,
) -> (
    SumfoldInstance<G>,
    <G as Curve>::Fr,
    MultilinearPolynomial<<G as Curve>::Fr>,
) {
    // Step 1: F = instances[0].f_func
    let f_cloned = instances[0].f_func.clone();

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
    // nu, l for interpolation
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
    // index = (b_val << l) + x_val
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
    let mut rng = rand::thread_rng();
    let rho_int = rng.gen_range(0, n); // for rand >= 0.7, use two-arg version e.g. rng.gen_range(0..n)
    let rho_field = <G as Curve>::Fr::from(rho_int as u64);

    // Step 7: call build_q_polynomial
    let q_b = build_q_polynomial(&f_js, &*f_cloned, rho_int, nu, l);

    // Step 8: return
    (
        instances[rho_int].clone(),
        rho_field,
        q_b
    )
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::Arc;

    use ark_bls12_381::{Bls12_381 as E, Fr as BLSFr};
    use ark_ff::{PrimeField, Zero};
    use merlin::Transcript;
    use crate::sumfold::multilinear::MultilinearPolynomial;
    use crate::sumcheck::SumCheckProof;

    /// small helper: build a multilinear poly of dimension l => size=2^l
    fn build_small_poly(l: usize) -> MultilinearPolynomial<BLSFr> {
        let size = 1 << l;
        let mut z = Vec::with_capacity(size);
        let mut val = BLSFr::zero();
        for _ in 0..size {
            val += BLSFr::from(1u64);
            z.push(val);
        }
        MultilinearPolynomial::new(z)
    }

    /// for demonstration: compute product poly
    // fn build_product_poly(
    //     g0: &MultilinearPolynomial<BLSFr>,
    //     g1: &MultilinearPolynomial<BLSFr>,
    // ) -> MultilinearPolynomial<BLSFr>
    // {
    //     assert_eq!(g0.len(), g1.len());
    //     let size = g0.len();
    //     let mut z = g0.z.clone();
    //     for i in 0..size {
    //         z[i] *= g1.z[i];
    //     }
    //     MultilinearPolynomial::new(z)
    // }

    /// A basic test that uses sumfold to pick one instance's polynomial, build Q(b),
    /// and check Q(rho)= sum_x [F(g0(x), g1(x))].
    #[test]
    fn test_sumfold() {
        let n = 2;
        let l = 2;

        // define F as product
        let f_arc: Arc<dyn Fn(&[BLSFr])->BLSFr + Send + Sync> =
            Arc::new(|vals: &[BLSFr]| vals.iter().product());

        // build 2 instances
        let mut instances = Vec::with_capacity(n);
        for _ in 0..n {
            // build g0,g1
            let g0= build_small_poly(l);
            let g1= build_small_poly(l);
            // store
            let inst = SumfoldInstance {
                f_func: f_arc.clone(),
                g_vec: vec![g0, g1],
            };
            instances.push(inst);
        }

        // call sumfold
        let (chosen_inst, rho_field, q_b) = sumfold::<E>(instances);

        // 1) compute the actual sum_x F(g0(x), g1(x)) for chosen_inst
        let size = chosen_inst.g_vec[0].len();
        let mut t_val = BLSFr::zero();
        for i in 0..size {
            let val = (chosen_inst.f_func)(&[
                chosen_inst.g_vec[0].z[i],
                chosen_inst.g_vec[1].z[i],
            ]);
            t_val += val;
        }

        // 2) check Q(rho) == that sum
        let rho_usize = (rho_field.into_repr().as_ref()[0]) as usize;
        let qb_rho = q_b.z[rho_usize];
        assert_eq!(
            qb_rho, t_val,
            "q_b(rho) must match sum_x of F(g_vec)"
        );

        // 3) naive sumcheck of Q(b):
        let total_sum: BLSFr = q_b.z.iter().copied().sum();

        let mut q_b_clone = q_b.clone();

        let mut transcript = Transcript::new(b"test_sumcheck");
        let (proof, _challenges)= SumCheckProof::<E>::simple_prover(&mut q_b_clone, total_sum, &mut transcript);

        let mut verify_transcript = Transcript::new(b"test_sumcheck");
        let ok = proof.simple_verify(total_sum, &mut verify_transcript);
        assert!(ok, "Naive sumcheck on q_b should pass as well");
    }

}
