use ark_ff::{to_bytes, One, Zero};
use merlin::Transcript;
use std::sync::Arc;
use zkp_curve::Curve;

use crate::circuit::Circuit;
use crate::evaluate::{eval_output, eval_value, random_bytes_to_fr};
use crate::libra_linear_gkr::{initialize_phase_one, initialize_phase_two};
use crate::sumcheck::SumCheckProof;
use crate::sumfold::multilinear::MultilinearPolynomial;
use crate::sumfold::prover::{sumfold, SumfoldInstance};
use crate::Vec;

/// GKRFold implementation that compresses multiple GKR proofs (only multiplication gates) into a single SumCheckProof
///
///  - We ignore add_hg_vec, focusing only on `mul_hg_vec` (乗算ゲート) for demonstration.
///  - We forcibly unify polynomial length by zero-padding to the maximum bit_size among all layers.
///  - We also pad the number of polynomials to 2^ν so that sumfold doesn't panic.
///
/// # Arguments
///
/// * `circuits` - A vector of circuits
/// * `inputs` - A vector of inputs for each circuit
/// * `witnesses` - A vector of witnesses for each circuit
/// * `circuit_to_hash` - A vector of circuit hashes
///
/// # Returns
///
/// * `SumCheckProof` - A single SumCheckProof that represents all the GKR proofs
pub fn gkrfold_mul_only<G: Curve>(
    circuits: Vec<&Circuit>,
    inputs: Vec<&Vec<G::Fr>>,
    witnesses: Vec<Vec<G::Fr>>,
    circuit_to_hash: Vec<G::Fr>,
) -> SumCheckProof<G> {
    let mut transcript = Transcript::new(b"libra - gkrfold-mul-only");

    // 1) write circuit_to_hash
    for hash in &circuit_to_hash {
        transcript.append_message(b"circuit_to_hash", &to_bytes!(hash).unwrap());
    }

    // 2) evaluate circuits
    let mut circuit_evals = Vec::new();
    let mut max_bit_size = 0; // track the max bit_size across all layers
    for (i, circuit) in circuits.iter().enumerate() {
        let evals = circuit.evaluate::<G>(inputs[i], &witnesses[i]).unwrap();

        transcript.append_message(b"input", &to_bytes!(evals[0]).unwrap());
        transcript.append_message(b"output", &to_bytes!(evals[evals.len() - 1]).unwrap());

        // find max bit_size
        for d in 0..circuit.depth {
            let bsize = circuit.layers[d].bit_size;
            if bsize > max_bit_size {
                max_bit_size = bsize;
            }
        }

        circuit_evals.push(evals);
    }

    // We'll unify all polynomials to length = 2^max_bit_size
    let unified_len = 1 << max_bit_size;

    // 3) build SumFoldInstances
    let mut sumfold_instances = Vec::new();

    for (cir_idx, circuit) in circuits.iter().enumerate() {
        let evals = &circuit_evals[cir_idx];

        let mut alpha = G::Fr::one();
        let mut beta = G::Fr::zero();

        // Evaluate output
        let (mut result_u, mut gu) = eval_output::<G>(
            &evals[circuit.depth - 1],
            circuit.layers[circuit.depth - 1].bit_size,
            &mut transcript,
        );
        let mut gv = vec![G::Fr::zero(); gu.len()];
        let mut result_v = G::Fr::zero();

        // from top layer down
        for d in (1..circuit.depth).rev() {
            // skip the claim stuff
            let uv_size = circuit.layers[d - 1].bit_size;

            // Phase1 => only take mul_hg_vec
            let (mul_hg_vec, _add1, _add2) = initialize_phase_one::<G>(
                &gu, &gv,
                &circuit.layers[d].gates,
                &evals[d - 1],
                uv_size,
                alpha, beta,
            );
            // => mul_hg_vec is dimension 2^uv_size
            // unify dimension
            let mul_hg_vec_padded = pad_with_zero::<G>(mul_hg_vec, unified_len);

            // Make instance for phase1
            sumfold_instances.push(
                create_sumfold_instance_phase::<G>(&mul_hg_vec_padded)
            );

            // random ru
            let mut ru = Vec::new();
            for _ in 0..uv_size {
                let mut buf = [0u8; 32];
                transcript.challenge_bytes(b"challenge_nextround", &mut buf);
                let r_u = random_bytes_to_fr::<G>(&buf);
                ru.push(r_u);
            }

            // Phase2 => we only care about mul_hg_vec from phase2
            let (mul_hg2, _add_hg2, fu) = initialize_phase_two::<G>(
                &gu, &gv, &ru,
                &circuit.layers[d].gates,
                &evals[d - 1],
                uv_size,
                alpha, beta,
            );
            let mul_hg2_padded = pad_with_zero::<G>(mul_hg2, unified_len);

            // sumfold instance for phase2
            sumfold_instances.push(
                create_sumfold_instance_phase::<G>(&mul_hg2_padded)
            );

            // random rv
            let mut rv = Vec::new();
            for _ in 0..uv_size {
                let mut buf = [0u8; 32];
                transcript.challenge_bytes(b"challenge_nextround", &mut buf);
                let r_v = random_bytes_to_fr::<G>(&buf);
                rv.push(r_v);
            }

            // update
            if d > 1 {
                gu = ru;
                gv = rv;
                result_u = fu;
                result_v = eval_value::<G>(&evals[d-1], &gv);

                let mut buf = [0u8; 32];
                transcript.challenge_bytes(b"challenge_alpha", &mut buf);
                alpha = random_bytes_to_fr::<G>(&buf);

                let mut buf2 = [0u8; 32];
                transcript.challenge_bytes(b"challenge_beta", &mut buf2);
                beta = random_bytes_to_fr::<G>(&buf2);
            }
        }
    }

    if sumfold_instances.is_empty() {
        panic!("No SumFold instances");
    }

    // 4) SumFold には "polynomial 数が 2^ν" 必須なので、必要なら dummy 追加
    let needed_power2 = next_power_of_two(sumfold_instances.len());
    while sumfold_instances.len() < needed_power2 {
        sumfold_instances.push(dummy_sumfold_instance::<G>(unified_len));
    }

    // 5) run sumfold
    let (_chosen, _rho_field, q_b) = sumfold::<G>(sumfold_instances);

    // 6) final sumcheck => total sum of Q(b)
    let total_sum: G::Fr = q_b.z.iter().copied().sum();
    let mut poly_clone = q_b.clone();
    let mut final_transcript = Transcript::new(b"gkrfold_sumcheck - mulonly");
    let (proof, _) = SumCheckProof::<G>::prove(&mut poly_clone, total_sum, &mut final_transcript);

    proof
}

// ============ Helpers ==============

/// zero-pad a 1D array up to length = new_len
fn pad_with_zero<G: Curve>(orig: Vec<G::Fr>, new_len: usize) -> Vec<G::Fr> {
    let mut v = orig;
    if v.len() < new_len {
        v.resize(new_len, G::Fr::zero());
    }
    v
}

/// find next power-of-two >= x
fn next_power_of_two(mut x: usize) -> usize {
    if x == 0 {
        return 1;
    }
    x -= 1;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    // up to 64 bits if needed
    x += 1;
    x
}

/// create dummy polynomial => dimension= new_len => all zero
/// for sumfold: g0*g1 => 0 * 0 => 0
fn dummy_sumfold_instance<G: Curve>(new_len: usize) -> SumfoldInstance<G> {
    let zeros = vec![G::Fr::zero(); new_len];
    let g0 = MultilinearPolynomial::new(zeros.clone());
    let g1 = MultilinearPolynomial::new(zeros);

    let f_prod: Arc<dyn Fn(&[G::Fr]) -> G::Fr + Send + Sync> =
        Arc::new(|vals: &[G::Fr]| vals[0]*vals[1]);

    SumfoldInstance {
        F_func: f_prod,
        g_vec: vec![g0, g1],
    }
}

/// Phase1 sumfold => g0= mul_hg_vec, g1= 1 array => but we only keep mul_hg_vec
fn create_sumfold_instance_phase<G: Curve>(poly: &Vec<G::Fr>) -> SumfoldInstance<G> {
    let g0 = MultilinearPolynomial::new(poly.clone());
    let ones = vec![G::Fr::one(); poly.len()];
    let g1 = MultilinearPolynomial::new(ones);

    let f_prod: Arc<dyn Fn(&[G::Fr]) -> G::Fr + Send + Sync> =
        Arc::new(|vals: &[G::Fr]| vals[0]*vals[1]);

    SumfoldInstance {
        F_func: f_prod,
        g_vec: vec![g0, g1],
    }
}
