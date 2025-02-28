use ark_ff::{to_bytes, Field};
use ark_poly::{polynomial::univariate::DensePolynomial, Polynomial, UVPolynomial};
use ark_serialize::*;
use merlin::Transcript;
use zkp_curve::Curve;

use crate::evaluate::{combine_with_r, random_bytes_to_fr};
use crate::polynomial_to_bytes;
use crate::Vec;

#[derive(CanonicalSerialize, CanonicalDeserialize, Clone)]
pub struct SumCheckProof<G: Curve> {
    pub polys: Vec<DensePolynomial<G::Fr>>,
    pub poly_value_at_r: Vec<G::Fr>,
}

// 乗算についてのメモ
// f_vcは レイヤーi+1のV_i+1(x)を並べたベクトル
// mul_hg_vec はg_vec.0と等しく、これはmulの部分和 h_i+1(x)=Σ_y mul(z, x, y)V_i+1(y)) 

impl<G: Curve> SumCheckProof<G> {
    pub fn phase_one_prover(
        f_vec: &Vec<G::Fr>,
        g_vec: &(Vec<G::Fr>, Vec<G::Fr>, Vec<G::Fr>),
        bit_size: usize,
        claim: G::Fr,
        transcript: &mut Transcript,
    ) -> (Self, Vec<G::Fr>) {
        let mut size = f_vec.len();
        let (mut mul_hg_vec, mut add_hg_vec1, mut add_hg_vec2) =
            (g_vec.0.clone(), g_vec.1.clone(), g_vec.2.clone());
        let mut f_vec = f_vec.clone();
        assert_eq!(size, g_vec.0.len());
        assert_eq!(size, g_vec.1.len());
        assert_eq!(size, g_vec.2.len());
        assert_eq!(size, (2usize).pow(bit_size as u32));
        let fr_two: G::Fr = 2u32.into();

        let mut claim = claim;
        let mut ru = Vec::new();
        let mut polys = Vec::new();

        for _ in 0..bit_size {
            size /= 2;
            let eval_0: G::Fr = (0..size)
                .map(|j| {
                    // 1番目の項がh_i+1(x)に対応？
                    f_vec[j] * &mul_hg_vec[j] + &(f_vec[j] * &add_hg_vec1[j]) + &add_hg_vec2[j]
                })
                .sum();
            let eval_1 = claim - &eval_0;

            let f_vec_tmp = combine_with_r::<G>(&f_vec, fr_two);
            let mul_hg_vec_tmp = combine_with_r::<G>(&mul_hg_vec, fr_two);
            let add_hg_vec1_tmp = combine_with_r::<G>(&add_hg_vec1, fr_two);
            let add_hg_vec2_tmp = combine_with_r::<G>(&add_hg_vec2, fr_two);

            let eval_2: G::Fr = (0..size)
                .map(|j| {
                    f_vec_tmp[j] * &mul_hg_vec_tmp[j]
                        + &(f_vec_tmp[j] * &add_hg_vec1_tmp[j])
                        + &add_hg_vec2_tmp[j]
                })
                .sum();

            // degree = 2
            // f(x) = ax^2 + bx + c
            // a = (eval_0 - 2eval_1 + eval_2)/2
            let a_coeff = (eval_0 - &eval_1.double() + &eval_2) * &fr_two.inverse().unwrap();
            // c = eval_0
            let c_coeff = eval_0;
            // b = eval_1 - a - c
            let b_coeff = eval_1 - &a_coeff - &c_coeff;
            let poly = DensePolynomial::from_coefficients_vec(vec![c_coeff, b_coeff, a_coeff]);

            transcript.append_message(b"poly", &polynomial_to_bytes::<G>(&poly));
            let mut buf = [0u8; 32];
            transcript.challenge_bytes(b"challenge_nextround", &mut buf);
            let r_i = random_bytes_to_fr::<G>(&buf);

            mul_hg_vec = combine_with_r::<G>(&mul_hg_vec, r_i);
            add_hg_vec1 = combine_with_r::<G>(&add_hg_vec1, r_i);
            add_hg_vec2 = combine_with_r::<G>(&add_hg_vec2, r_i);
            f_vec = combine_with_r::<G>(&f_vec, r_i);

            claim = poly.evaluate(&r_i);
            ru.push(r_i);
            polys.push(poly);
        }

        let poly_value_at_r = vec![f_vec[0], mul_hg_vec[0], add_hg_vec1[0], add_hg_vec2[0]];
        transcript.append_message(b"claim_final", &to_bytes!(poly_value_at_r).unwrap());

        let proof = Self {
            polys,
            poly_value_at_r,
        };
        (proof, ru)
    }

    pub fn phase_two_prover(
        f_vec: &Vec<G::Fr>,
        g_vec: &(Vec<G::Fr>, Vec<G::Fr>, G::Fr),
        bit_size: usize,
        claim: G::Fr,
        transcript: &mut Transcript,
    ) -> (Self, Vec<G::Fr>) {
        let mut size = f_vec.len();
        let (mut mul_hg_vec, mut add_hg_vec, fu) = (g_vec.0.clone(), g_vec.1.clone(), g_vec.2);
        let mut f_vec = f_vec.clone();
        assert_eq!(size, g_vec.0.len());
        assert_eq!(size, g_vec.1.len());
        assert_eq!(size, (2usize).pow(bit_size as u32));
        let fr_two: G::Fr = 2u32.into();

        let mut claim = claim;
        let mut rv = Vec::new();
        let mut polys = Vec::new();

        for _ in 0..bit_size {
            size /= 2;
            let eval_0: G::Fr = (0..size)
                .map(|j| {
                    mul_hg_vec[j] * &f_vec[j] * &fu
                        + &(add_hg_vec[j] * &fu)
                        + &(add_hg_vec[j] * &f_vec[j])
                })
                .sum();
            let eval_1 = claim - &eval_0;

            let f_vec_tmp = combine_with_r::<G>(&f_vec, fr_two);
            let mul_hg_vec_tmp = combine_with_r::<G>(&mul_hg_vec, fr_two);
            let add_hg_vec_tmp = combine_with_r::<G>(&add_hg_vec, fr_two);
            // let add_hg_vec2_tmp = combine_with_r::<G>(&add_hg_vec2, fr_two);
            let eval_2: G::Fr = (0..size)
                .map(|j| {
                    mul_hg_vec_tmp[j] * &f_vec_tmp[j] * &fu
                        + &(add_hg_vec_tmp[j] * &fu)
                        + &(add_hg_vec_tmp[j] * &f_vec_tmp[j])
                })
                .sum();

            // degree = 2
            // f(x) = ax^2 + bx + c
            // a = (eval_0 - 2eval_1 + eval_2)/2
            let a_coeff = (eval_0 - &eval_1.double() + &eval_2) * &fr_two.inverse().unwrap();
            // c = eval_0
            let c_coeff = eval_0;
            // b = eval_1 - a - c
            let b_coeff = eval_1 - &a_coeff - &c_coeff;
            let poly = DensePolynomial::from_coefficients_vec(vec![c_coeff, b_coeff, a_coeff]);

            transcript.append_message(b"poly", &polynomial_to_bytes::<G>(&poly));
            let mut buf = [0u8; 32];
            transcript.challenge_bytes(b"challenge_nextround", &mut buf);
            let r_i = random_bytes_to_fr::<G>(&buf);
            mul_hg_vec = combine_with_r::<G>(&mul_hg_vec, r_i);
            add_hg_vec = combine_with_r::<G>(&add_hg_vec, r_i);
            f_vec = combine_with_r::<G>(&f_vec, r_i);

            claim = poly.evaluate(&r_i);
            rv.push(r_i);
            polys.push(poly);
        }

        let poly_value_at_r = vec![f_vec[0], mul_hg_vec[0], add_hg_vec[0]];
        transcript.append_message(b"claim_final", &to_bytes!(poly_value_at_r).unwrap());

        let proof = Self {
            polys,
            poly_value_at_r,
        };

        (proof, rv)
    }
}

