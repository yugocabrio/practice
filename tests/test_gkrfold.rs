use ark_bls12_381::{Bls12_381 as E, Fr};
use ark_ff::One;
use std::time::Instant;

use zkp_libra::{
    circuit::Circuit, gkrfold::gkrfold, libra_linear_gkr::LinearGKRProof,
};

/// Create a simple circuit for testing
///
/// input0 * witness0 + input1 * witness1 = output
fn prepare_simple_circuit() -> (Vec<Fr>, Vec<Fr>, Vec<Vec<(u8, usize, usize)>>) {
    let inputs = vec![
        Fr::one() + &Fr::one(),              // 2
        Fr::one() + &Fr::one() + &Fr::one(), // 3
    ];

    let witnesses = vec![
        Fr::one() + &Fr::one() + &Fr::one(),              // 3
        Fr::one() + &Fr::one() + &Fr::one() + &Fr::one(), // 4
    ];

    let mut layers = Vec::new();

    // Layer 1: Multiplication gates
    let mut layer = Vec::new();
    layer.push((1, 0, 2)); // witness0 * input0 = 3 * 2 = 6
    layer.push((1, 1, 3)); // witness1 * input1 = 4 * 3 = 12
    layers.push(layer);

    // Layer 2: Addition gate
    let mut layer = Vec::new();
    layer.push((0, 0, 1)); // 6 + 12 = 18
    layers.push(layer);

    (inputs, witnesses, layers)
}

/// Create another simple circuit for testing
///
/// input0 * witness0 - input1 * witness1 = output
/// (implemented as addition with a negative value)
fn prepare_another_circuit() -> (Vec<Fr>, Vec<Fr>, Vec<Vec<(u8, usize, usize)>>) {
    let inputs = vec![
        Fr::one() + &Fr::one() + &Fr::one() + &Fr::one(), // 4
        Fr::one() + &Fr::one(),                           // 2
    ];

    let witnesses = vec![
        Fr::one() + &Fr::one(),                           // 2
        Fr::one() + &Fr::one() + &Fr::one() + &Fr::one(), // 4
    ];

    // For subtraction, we'll multiply the second value by -1 and then add
    // We can't directly represent -1, so we'll handle this in the circuit evaluation

    let mut layers = Vec::new();

    // Layer 1: Multiplication gates
    let mut layer = Vec::new();
    layer.push((1, 0, 2)); // witness0 * input0 = 2 * 4 = 8
    layer.push((1, 1, 3)); // witness1 * input1 = 4 * 2 = 8
    layers.push(layer);

    // Layer 2: Addition gate (we'll treat the second value as negative during evaluation)
    let mut layer = Vec::new();
    layer.push((0, 0, 1)); // 8 - 8 = 0 (but represented as 8 + (-8))
    layers.push(layer);

    (inputs, witnesses, layers)
}

/// Create a simple Poseidon hash circuit
///
/// This is a simplified version that just does some additions and multiplications
/// to simulate the Poseidon hash function
fn prepare_poseidon_circuit() -> (Vec<Fr>, Vec<Fr>, Vec<Vec<(u8, usize, usize)>>) {
    // For simplicity, we'll use 2 inputs (representing the state of the hash)
    let inputs = vec![
        Fr::one(),                           // 1
        Fr::one() + &Fr::one(),              // 2
    ];

    // And 2 witnesses (representing round constants)
    let witnesses = vec![
        Fr::one() + &Fr::one() + &Fr::one(),              // 3
        Fr::one() + &Fr::one() + &Fr::one() + &Fr::one(), // 4
    ];

    let mut layers = Vec::new();

    // Layer 1: First round - Add round constants
    let mut layer = Vec::new();
    layer.push((0, 0, 2)); // input0 + witness0 = 1 + 3 = 4
    layer.push((0, 1, 3)); // input1 + witness1 = 2 + 4 = 6
    layers.push(layer);

    // Layer 2: Apply S-box (multiplication)
    let mut layer = Vec::new();
    layer.push((1, 0, 0)); // result0 * result0 = 4 * 4 = 16
    layer.push((1, 1, 1)); // result1 * result1 = 6 * 6 = 36
    layers.push(layer);

    // Layer 3: Mix layer (addition)
    let mut layer = Vec::new();
    layer.push((0, 0, 1)); // result0 + result1 = 16 + 36 = 52
    layers.push(layer);

    (inputs, witnesses, layers)
}

#[test]
fn test_gkrfold() {
    println!("Starting GKRFold test...");

    // Prepare the first circuit
    let (inputs1, witnesses1, layers1) = prepare_simple_circuit();
    let circuit1 = Circuit::new(inputs1.len(), witnesses1.len(), &layers1);
    let circuit_to_hash1 = circuit1.circuit_to_hash::<E>();

    // Prepare the second circuit
    let (inputs2, witnesses2, layers2) = prepare_another_circuit();
    let circuit2 = Circuit::new(inputs2.len(), witnesses2.len(), &layers2);
    let circuit_to_hash2 = circuit2.circuit_to_hash::<E>();

    println!("Circuits prepared successfully");

    // Generate individual proofs (not necessary for GKRFold, but useful for verification)
    let (_proof1, output1) = LinearGKRProof::<E>::prover(&circuit1, &inputs1, &witnesses1, circuit_to_hash1);
    let (_proof2, output2) = LinearGKRProof::<E>::prover(&circuit2, &inputs2, &witnesses2, circuit_to_hash2);

    println!("Individual proofs generated successfully");
    println!("Circuit 1 output: {:?}", output1);
    println!("Circuit 2 output: {:?}", output2);

    // Now use GKRFold to compress the proofs
    let folded_proof = gkrfold::<E>(
        vec![&circuit1, &circuit2],
        vec![&inputs1, &inputs2],
        vec![witnesses1.clone(), witnesses2.clone()],
        vec![circuit_to_hash1, circuit_to_hash2],
    );

    println!("GKRFold proof generated successfully");

    // Verify that the folded proof was created successfully
    assert!(folded_proof.polys.len() > 0, "Folded proof should have polynomials");

    println!("GKRFold test completed successfully");
}

#[test]
fn test_gkrfold_multiple_circuits() {
    println!("Starting GKRFold test with multiple circuits...");

    // Test with 4 circuits
    test_gkrfold_with_n_circuits(4);

    // Test with 8 circuits
    test_gkrfold_with_n_circuits(8);

    // Test with 16 circuits
    test_gkrfold_with_n_circuits(16);

    println!("All multiple circuit tests completed successfully");
}

/// Helper function to test GKRFold with N circuits
/// This uses the simple circuit creation functions we already have
fn test_gkrfold_with_n_circuits(n: usize) {
    println!("Testing GKRFold with {} circuits...", n);

    let start = Instant::now();

    // Create vectors to hold all circuits and their data
    let mut circuits = Vec::with_capacity(n);
    let mut all_inputs = Vec::with_capacity(n);
    let mut all_witnesses = Vec::with_capacity(n);
    let mut all_circuit_to_hash = Vec::with_capacity(n);

    // Create n circuits by alternating between the two circuit types we already have
    for i in 0..n {
        if i % 2 == 0 {
            // Even index: use simple circuit
            let (inputs, witnesses, layers) = prepare_simple_circuit();
            let circuit = Circuit::new(inputs.len(), witnesses.len(), &layers);
            let circuit_to_hash = circuit.circuit_to_hash::<E>();

            circuits.push(circuit);
            all_inputs.push(inputs);
            all_witnesses.push(witnesses);
            all_circuit_to_hash.push(circuit_to_hash);
        } else {
            // Odd index: use another circuit
            let (inputs, witnesses, layers) = prepare_another_circuit();
            let circuit = Circuit::new(inputs.len(), witnesses.len(), &layers);
            let circuit_to_hash = circuit.circuit_to_hash::<E>();

            circuits.push(circuit);
            all_inputs.push(inputs);
            all_witnesses.push(witnesses);
            all_circuit_to_hash.push(circuit_to_hash);
        }
    }

    println!("Generated {} circuits in {:?}", n, start.elapsed());

    // Prepare references for gkrfold
    let circuit_refs: Vec<&Circuit> = circuits.iter().collect();
    let input_refs: Vec<&Vec<Fr>> = all_inputs.iter().collect();

    // Use GKRFold to compress the proofs
    let folded_proof = gkrfold::<E>(
        circuit_refs,
        input_refs,
        all_witnesses.clone(),
        all_circuit_to_hash,
    );

    println!("GKRFold proof for {} circuits generated in {:?}", n, start.elapsed());

    // Verify that the folded proof was created successfully
    assert!(folded_proof.polys.len() > 0, "Folded proof should have polynomials");

    println!("GKRFold test with {} circuits completed successfully", n);
}

#[test]
fn test_gkrfold_poseidon() {
    println!("Starting GKRFold test with Poseidon hash circuits...");

    // Prepare two Poseidon hash circuits
    let (inputs1, witnesses1, layers1) = prepare_poseidon_circuit();
    let circuit1 = Circuit::new(inputs1.len(), witnesses1.len(), &layers1);
    let circuit_to_hash1 = circuit1.circuit_to_hash::<E>();

    // Create a second Poseidon circuit with slightly different inputs
    let mut inputs2 = inputs1.clone();
    inputs2[0] = Fr::one() + &Fr::one(); // Change first input to 2
    let witnesses2 = witnesses1.clone();
    let circuit2 = Circuit::new(inputs2.len(), witnesses2.len(), &layers1);
    let circuit_to_hash2 = circuit2.circuit_to_hash::<E>();

    println!("Poseidon circuits prepared successfully");

    // Generate individual proofs
    let (_proof1, output1) = LinearGKRProof::<E>::prover(&circuit1, &inputs1, &witnesses1, circuit_to_hash1);
    let (_proof2, output2) = LinearGKRProof::<E>::prover(&circuit2, &inputs2, &witnesses2, circuit_to_hash2);

    println!("Individual Poseidon proofs generated successfully");
    println!("Poseidon Circuit 1 output: {:?}", output1);
    println!("Poseidon Circuit 2 output: {:?}", output2);

    // Now use GKRFold to compress the proofs
    let folded_proof = gkrfold::<E>(
        vec![&circuit1, &circuit2],
        vec![&inputs1, &inputs2],
        vec![witnesses1.clone(), witnesses2.clone()],
        vec![circuit_to_hash1, circuit_to_hash2],
    );

    println!("GKRFold proof for Poseidon circuits generated successfully");

    // Verify that the folded proof was created successfully
    assert!(folded_proof.polys.len() > 0, "Folded proof should have polynomials");

    println!("GKRFold Poseidon test completed successfully");
}
