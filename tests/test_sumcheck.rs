use zkp_libra::sumcheck::SumCheckProof;
use std::time::Instant;
use ark_ff::{Field, One};
use ark_bls12_381::Fr;
use merlin::Transcript;
use zkp_libra::sumfold::multilinear::MultilinearPolynomial;
use ark_bls12_381::Bls12_381;
use rand::{rngs::ThreadRng, thread_rng, Rng};

fn build_random_poly<F: Field, R: Rng>(n: usize, rng: &mut R) -> MultilinearPolynomial<F> {
    let size = 1 << n;

    MultilinearPolynomial::new(
        (0..size)
          .into_iter()
          .map(|_| {
              F::from(rng.gen_range(1, 50) as u64)
          })
          .collect())
}

#[test]
fn test_sumcheck_correctness() {
    let ns = [2, 4, 8, 16];
    let num_tries = 100;
    for &n in &ns {
        for _ in 0..num_tries {
            let mut rng = thread_rng();
            let poly = build_random_poly::<Fr, ThreadRng>(n, &mut rng);
            let claimed_sum: Fr = poly.z.iter().copied().sum();

            let mut transcript = Transcript::new(b"test_sumcheck");

            let start = Instant::now();
            let (proof, _challenges) = SumCheckProof::<Bls12_381>::prove(
                &mut poly.clone(),
                claimed_sum,
                &mut transcript,
            );
            let prove_time = start.elapsed().as_secs_f64() * 1000.0;

            let mut verify_transcript = Transcript::new(b"test_sumcheck");
            let start = Instant::now();
            let ok = proof.verify(claimed_sum, &mut verify_transcript);
            let verify_time = start.elapsed().as_secs_f64() * 1000.0;
            assert!(ok, "sumcheck pass");
            println!("n = {}, prove time = {:.2} ms, verify time = {:.2} ms", n, prove_time, verify_time);
        }
    }
}

#[test]
fn test_sumcheck_soundness() {
    let ns = [2, 4, 8, 16];
    let num_tries = 100;
    for &n in &ns {
        for _ in 0..num_tries {
            let mut rng = thread_rng();
            let poly = build_random_poly::<Fr, ThreadRng>(n, &mut rng);
            let claimed_sum: Fr = poly.z.iter().copied().sum();

            let mut transcript = Transcript::new(b"test_sumcheck");

            let (proof, _challenges) = SumCheckProof::<Bls12_381>::prove(
                &mut poly.clone(),
                claimed_sum,
                &mut transcript,
            );

            // Verify invalid claim
            let mut verify_transcript = Transcript::new(b"test_sumcheck");
            let invalid_claim = claimed_sum + Fr::one();
            let ok = proof.verify(invalid_claim, &mut verify_transcript);
            assert!(!ok, "sumcheck fail");

            // Verify invalid proof
            let mut verify_transcript = Transcript::new(b"test_sumcheck");
            let invalid_proof = SumCheckProof::<Bls12_381> {
                polys: proof.polys.clone(),
                poly_value_at_r: proof.poly_value_at_r.clone().into_iter().map(|x| x + Fr::one()).collect(),
            };
            let ok = invalid_proof.verify(claimed_sum, &mut verify_transcript);
            assert!(!ok, "sumcheck fail");

            // Verify invalid transcript
            let mut verify_transcript = Transcript::new(b"test_sumcheck_invalid");
            let ok = proof.verify(claimed_sum, &mut verify_transcript);
            assert!(!ok, "sumcheck fail");
        }
    }
}