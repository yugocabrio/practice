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

/// GKRFold implementation that compresses multiple GKR proofs into a single SumCheckProof
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
pub fn gkrfold<G: Curve>(
    circuits: Vec<&Circuit>,
    inputs: Vec<&Vec<G::Fr>>,
    witnesses: Vec<Vec<G::Fr>>,
    circuit_to_hash: Vec<G::Fr>,
) -> SumCheckProof<G> {
    // 1. Set up the transcript (similar to libra_linear_gkr.rs:L28~36)
    let mut transcript = Transcript::new(b"libra - gkrfold");

    // Append all circuit_to_hash values to the transcript
    for hash in &circuit_to_hash {
        transcript.append_message(b"circuit_to_hash", &to_bytes!(hash).unwrap());
    }

    // Evaluate all circuits
    let mut circuit_evals = Vec::new();
    for (i, circuit) in circuits.iter().enumerate() {
        let evals = circuit.evaluate::<G>(inputs[i], &witnesses[i]).unwrap();

        // Append input and output to transcript
        transcript.append_message(b"input", &to_bytes!(evals[0]).unwrap());
        transcript.append_message(
            b"output",
            &to_bytes!(evals[evals.len() - 1]).unwrap(),
        );

        circuit_evals.push(evals);
    }

    // 2. For each circuit, collect SumFoldInstances
    let mut sumfold_instances = Vec::new();

    for (circuit_idx, circuit) in circuits.iter().enumerate() {
        let evals = &circuit_evals[circuit_idx];

        let mut alpha = G::Fr::one();
        let mut beta = G::Fr::zero();

        // V_0(g^(0)), g^(0)
        let (mut result_u, mut gu) = eval_output::<G>(
            &evals[evals.len() - 1],
            circuit.layers[circuit.depth - 1].bit_size,
            &mut transcript,
        );

        let mut gv = vec![G::Fr::zero(); gu.len()];
        let mut result_v = G::Fr::zero();

        // For each depth in the circuit
        for d in (1..circuit.depth).rev() {
            let _claim = alpha * &result_u + &(beta * &result_v);
            let uv_size = circuit.layers[d - 1].bit_size;

            // Phase 1: Create SumFoldInstances for the first phase
            let (mul_hg_vec, add_hg_vec1, add_hg_vec2) = initialize_phase_one::<G>(
                &gu,
                &gv,
                &circuit.layers[d].gates,
                &evals[d - 1],
                uv_size,
                alpha,
                beta,
            );

            // Create 3 SumFoldInstances for phase one
            // 1. For mul_hg_vec
            let mul_hg_instance = create_sumfold_instance::<G>(&evals[d-1], &mul_hg_vec);
            sumfold_instances.push(mul_hg_instance);

            // 2. For add_hg_vec1
            let add_hg_instance1 = create_sumfold_instance::<G>(&evals[d-1], &add_hg_vec1);
            sumfold_instances.push(add_hg_instance1);

            // 3. For add_hg_vec2
            let add_hg_instance2 = create_sumfold_instance::<G>(&evals[d-1], &add_hg_vec2);
            sumfold_instances.push(add_hg_instance2);

            // Generate a random point for phase one
            let mut ru = Vec::new();
            for _ in 0..uv_size {
                let mut buf = [0u8; 32];
                transcript.challenge_bytes(b"challenge_nextround", &mut buf);
                let r_u = random_bytes_to_fr::<G>(&buf);
                ru.push(r_u);
            }

            // Phase 2: Create SumFoldInstances for the second phase
            let (mul_hg_vec, add_hg_vec, fu) = initialize_phase_two::<G>(
                &gu,
                &gv,
                &ru,
                &circuit.layers[d].gates,
                &evals[d - 1],
                uv_size,
                alpha,
                beta,
            );

            // Create 3 SumFoldInstances for phase two
            // 4. For mul_hg_vec in phase two
            let mul_hg_instance2 = create_sumfold_instance::<G>(&evals[d-1], &mul_hg_vec);
            sumfold_instances.push(mul_hg_instance2);

            // 5. For add_hg_vec in phase two
            let add_hg_instance = create_sumfold_instance::<G>(&evals[d-1], &add_hg_vec);
            sumfold_instances.push(add_hg_instance);

            // 6. For fu (constant value)
            let fu_vec = vec![fu; evals[d-1].len()];
            let fu_instance = create_sumfold_instance::<G>(&evals[d-1], &fu_vec);
            sumfold_instances.push(fu_instance);

            // Generate a random point for phase two
            let mut rv = Vec::new();
            for _ in 0..uv_size {
                let mut buf = [0u8; 32];
                transcript.challenge_bytes(b"challenge_nextround", &mut buf);
                let r_v = random_bytes_to_fr::<G>(&buf);
                rv.push(r_v);
            }

            // Update for next iteration
            if d > 1 {
                gu = ru.clone();
                gv = rv.clone();
                result_u = fu;
                result_v = eval_value::<G>(&evals[d-1], &rv);

                // Generate new alpha and beta
                let mut buf = [0u8; 32];
                transcript.challenge_bytes(b"challenge_alpha", &mut buf);
                alpha = random_bytes_to_fr::<G>(&buf);

                let mut buf = [0u8; 32];
                transcript.challenge_bytes(b"challenge_beta", &mut buf);
                beta = random_bytes_to_fr::<G>(&buf);
            }
        }
    }

    // Make sure we have at least one instance
    if sumfold_instances.is_empty() {
        panic!("No SumFold instances were created");
    }

    // 3. Call sumfold on all instances
    let (_chosen_instance, _rho_field, q_b) = sumfold::<G>(sumfold_instances);

    // 4. Return the result of SumCheckProof::prove()
    let mut q_b_clone = q_b.clone();
    let total_sum: G::Fr = q_b.z.iter().copied().sum();

    let mut final_transcript = Transcript::new(b"gkrfold_sumcheck");
    let (proof, _) = SumCheckProof::<G>::prove(&mut q_b_clone, total_sum, &mut final_transcript);

    proof
}

/// Helper function to create a SumfoldInstance from a vector and a polynomial
fn create_sumfold_instance<G: Curve>(
    values: &Vec<G::Fr>,
    poly: &Vec<G::Fr>,
) -> SumfoldInstance<G> {
    // Create multilinear polynomials
    let g0 = MultilinearPolynomial::new(poly.clone());
    let g1 = MultilinearPolynomial::new(values.clone());

    // Create the function that computes the product
    let f_func: Arc<dyn Fn(&[G::Fr]) -> G::Fr + Send + Sync> =
        Arc::new(|vals: &[G::Fr]| vals[0] * vals[1]);

    SumfoldInstance {
        F_func: f_func,
        g_vec: vec![g0, g1],
    }
}
