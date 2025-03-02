use ark_ff::Field;
use crate::sumfold::multilinear::MultilinearPolynomial;

#[derive(Clone)]
pub struct SimpleSumcheckProof<F: Field> {
    /// ラウンドごとにやりとりする1変数多項式
    pub univariate_polys: Vec<Vec<F>>,
    /// ラウンドでの乱数
    pub r_challenges: Vec<F>,
    /// sumcheck終了後の定数
    pub final_value: F,
}

pub struct SimpleSumcheck;

impl SimpleSumcheck {
    /// ProveではmultilinearPolynomialとlaimed_sum = Σ_x P(x) (x in {0,1}^l)を受け取る
    pub fn prove<F: Field>(
        poly: &mut MultilinearPolynomial<F>,
        claimed_sum: F,
    ) -> SimpleSumcheckProof<F> {
        let l = poly.get_num_vars();
        let mut current_claim = claimed_sum;

        let mut univariate_polys = Vec::new();
        let mut r_challenges = Vec::new();

        // 変数の数分のラウンド
        for round in 0..l {
            // P(0) P(1) に分ける
            let size = poly.len() / 2;
            let (poly_lo, poly_hi) = poly.z.split_at_mut(size);

            // P(0)とP(1)それぞれのsum
            let eval0: F = poly_lo.iter().copied().sum();
            let eval1: F = poly_hi.iter().copied().sum();

            // P0, P1 を満たす1次多項式 u(t) = eval0 + (eval1 - eval0)*t
            let c = eval0;
            let b = eval1 - eval0;
            univariate_polys.push(vec![c, b]);

            // 乱数
            // TODO: Fiat shamir
            let r_i = F::from(round as u64);
            r_challenges.push(r_i);

            // 変数xをrで固定して次元を減らす
            for i in 0..size {
                // new_poly[i] = poly_lo[i] + r_i * (poly_hi[i] - poly_lo[i])
                let tmp = poly_hi[i] - poly_lo[i];
                poly_lo[i] += r_i * tmp;
            }
            poly.z.truncate(size);

            // claimed_sum を u(r_i) に更新
            let next_claim = c + b * r_i;
            current_claim = next_claim;
        }

        SimpleSumcheckProof {
            univariate_polys,
            r_challenges,
            final_value: current_claim,
        }
    }

    pub fn verify<F: Field>(
        proof: &SimpleSumcheckProof<F>,
        claimed_sum: F,
        l: usize,
    ) -> bool {
        if proof.univariate_polys.len() != l {
            return false;
        }

        let mut round_claim = claimed_sum;
        let mut next_claim = F::zero();

        for (round, poly_coefs) in proof.univariate_polys.iter().enumerate() {
            // 1次多項式 (c + b*x) であるか
            if poly_coefs.len() != 2 {
                return false;
            }
            let c = poly_coefs[0];
            let b = poly_coefs[1];

            // poly(0)+poly(1) = round_claimを満たすか 
            if c + (c + b) != round_claim {
                return false;
            }

            // 次ラウンドを c + b*r_i に更新
            let r_i = proof.r_challenges[round];
            next_claim = c + b * r_i;
            round_claim = next_claim;
        }

        // Proverのfinal_valueと等しいか
        proof.final_value == next_claim
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, One, Zero};
    use ark_bls12_381::Fr;
    use crate::sumfold::multilinear::MultilinearPolynomial;

    // 評価点の計算
    fn build_small_poly<F: Field>(l: usize) -> MultilinearPolynomial<F> {
        let size = 1 << l;
        let mut z = Vec::with_capacity(size);

        let mut val = F::zero();
        for _ in 0..size {
            val += F::one();
            z.push(val);
        }
        MultilinearPolynomial::new(z)
    }

    #[test]
    fn test_simple_sumcheck_correct() {
        let l = 10;
        let poly = build_small_poly::<Fr>(l);

        let claimed_sum: Fr = poly.z.iter().copied().sum();

        let proof = SimpleSumcheck::prove(&mut poly.clone(), claimed_sum);

        let ok = SimpleSumcheck::verify(&proof, claimed_sum, l);
        assert!(ok, "sumcheck pass");
    }

    #[test]
    fn test_simple_sumcheck_incorrect_claim() {
        let l = 10;
        let poly = build_small_poly::<Fr>(l);

        let real_sum: Fr = poly.z.iter().copied().sum();

        // sumを水増し
        let claimed_sum = real_sum + Fr::one();

        let proof = SimpleSumcheck::prove(&mut poly.clone(), claimed_sum);

        let ok = SimpleSumcheck::verify(&proof, claimed_sum, l);
        assert!(!ok, "sumcheck fail");
    }
}
