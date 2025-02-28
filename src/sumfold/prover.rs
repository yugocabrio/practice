//! sumfold の "prover" 的な例
//! G: zkp_curve::Curve

use rand::Rng;
use zkp_curve::Curve; // あなたの定義
use ark_ff::{Zero};
use crate::sumfold::multilinear::MultilinearPolynomial;
use crate::sumfold::q_poly::build_q_polynomial;
use crate::sumcheck::SumCheckProof;
use std::sync::Arc;

/// sumcheck instance
#[derive(Clone)]
pub struct SumcheckInstance<G: Curve> {
    /// F: takes slice of G::Fr -> G::Fr
    pub F: Arc<dyn Fn(&[<G as Curve>::Fr]) -> <G as Curve>::Fr + Send + Sync>,
    /// マルチリニア多項式群
    pub g_vec: Vec<MultilinearPolynomial<<G as Curve>::Fr>>,
    /// sumcheckの証明
    pub proof: SumCheckProof<G>,
}

/// sumfold: ランダムにrho選んで Q(b) 構築
pub fn sumfold<G: Curve>(
    instances: Vec<SumcheckInstance<G>>,
) -> (
    SumcheckInstance<G>,
    <G as Curve>::Fr,
    MultilinearPolynomial<<G as Curve>::Fr>,
) {
    let n = instances.len();
    let mut rng = rand::thread_rng();
    // rand=0.7 => rng.gen_range(low, high)
    let rho_int = rng.gen_range(0, n);

    // step1
    let f_cloned = instances[0].F.clone();

    // step2
    let first_len = instances[0].g_vec.len();
    for inst in &instances {
        assert_eq!(inst.g_vec.len(), first_len);
    }

    // step3
    let t = instances[0].g_vec.len();
    let x = instances[0].g_vec[0].z.len();
    let nu = (n as f64).log2() as usize;
    let l  = (x as f64).log2() as usize;

    // step4
    let mut g_bj = Vec::with_capacity(n);
    for inst in &instances {
        g_bj.push(inst.g_vec.clone());
    }

    // step5
    let size = 1 << (nu + l);
    let mut f_js = Vec::with_capacity(t);
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

    // step7
    let q_b = build_q_polynomial(&f_js, &*f_cloned, rho_int, nu, l);

    // step8
    let rho_field = <G as Curve>::Fr::from(rho_int as u64);
    (instances[rho_int].clone(), rho_field, q_b)
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, PrimeField, UniformRand, Zero, One};
    use rand::rngs::StdRng;
    use rand::SeedableRng;
    use std::sync::Arc;

    // PairingEngine (ark) + zkp_curve::Curve
    use ark_bls12_381::Bls12_381 as E;
    // sumfold
    use crate::sumfold::multilinear::MultilinearPolynomial;
    use crate::sumcheck::SumCheckProof;

    #[test]
    fn test_sumfold() {
        type FE = <E as Curve>::Fr; // 型エイリアス: <Bls12_381 as zkp_curve::Curve>::Fr
        // (Bls12_381 は PairingEngine, かつ zkp_curve::Curve のimplあり)

        let n = 2;
        let t = 2;
        let x = 4;

        // F: product
        let f_arc: Arc<dyn Fn(&[FE]) -> FE + Send + Sync> =
            Arc::new(|vals| vals.iter().product::<FE>());

        let mut rng = StdRng::seed_from_u64(12345);
        let mut instances: Vec<SumcheckInstance<E>> = Vec::with_capacity(n);

        for _ in 0..n {
            // g_vec
            let mut g_vec = Vec::with_capacity(t);
            for _ in 0..t {
                let evals = (0..x).map(|_| FE::rand(&mut rng)).collect();
                g_vec.push(MultilinearPolynomial::new(evals));
            }

            // sumcheck proof dummy
            let proof = SumCheckProof::<E> {
                polys: Vec::new(),
                poly_value_at_r: Vec::new(),
            };

            let inst = SumcheckInstance {
                F: f_arc.clone(),
                g_vec,
                proof,
            };
            instances.push(inst);
        }

        let (chosen, rho_f, q_b) = sumfold::<E>(instances);

        // sum-of-products
        let mut T = FE::zero();
        for i in 0..x {
            T += (chosen.F)(&[
                chosen.g_vec[0].z[i],
                chosen.g_vec[1].z[i],
            ]);
        }

        // into_repr() => import ark_ff::PrimeField
        let rho_u = rho_f.into_repr().as_ref()[0] as usize;  
        //  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ `as_ref()` で BigInt sliceを取り、[0]番目

        assert_eq!(q_b.z[rho_u], T, "Q(rho) must match sum-of-products");
    }
}
