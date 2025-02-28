use ark_ff::{Field, One, Zero};
use crate::sumfold::multilinear::MultilinearPolynomial;

/// Given a set of polynomials for the j-th index:
/// g_{0,j}(x), g_{1,j}(x), ..., g_{(2^ν - 1), j}(x),
/// this function constructs and returns the multilinear extension (MLE)
/// defined by:
///   f_j(b,x) = Σ_{i in {0,1}^ν} eq(b, i) * g_{i,j}(x)
/// as an MLE in (ν + m) variables.
///
/// - gs_for_j: A slice of length 2^ν, where each element is a MultilinearPolynomial in m variables.
/// - Returns: A MultilinearPolynomial in (ν + m) variables (dense representation).
pub fn build_fj_polynomial<F: Field>(
    gs_for_j: &[MultilinearPolynomial<F>],
) -> MultilinearPolynomial<F> {
    let num_b = gs_for_j.len();
    let nu = (num_b as f64).log2() as usize;
    assert_eq!(1 << nu, num_b, "gs_for_j.len() must be 2^ν");

    let l = gs_for_j[0].get_num_vars();
    for b in 1..num_b {
        assert_eq!(
            gs_for_j[b].get_num_vars(),
            l,
            "all g_{{b,j}}(x) must have the same number of variables"
        );
    }

    let new_num_vars = nu + l;
    let new_len = 1 << new_num_vars;
    let mut f_j_evals = vec![F::zero(); new_len];

    let block_size = 1 << l;
    for (b, chunk) in f_j_evals.chunks_mut(block_size).enumerate() {
        // Copy gs_for_j[b].Z into the chunk
        chunk.copy_from_slice(&gs_for_j[b].z);
    }

    MultilinearPolynomial::new(f_j_evals)
}

/// Given decimal representations of b and x,
/// this function converts them internally into their bit representation (B1,...,Bν, X1,...,Xm)
/// and evaluates f_j(b,x).
///
/// - f: A MultilinearPolynomial in (ν + m) variables (constructed via build_fj_polynomial)
/// - b: Decimal representation of b. An integer with ν bits (0 <= b < 2^ν)
/// - x: Decimal representation of x. An integer with m bits (0 <= x < 2^m)
/// - nu: The number of bits required to represent b.
/// - m: The number of bits required to represent x.
pub fn evaluate_fj_at_decimals<F: Field>(
    f: &MultilinearPolynomial<F>,
    b: usize,
    x: usize,
    nu: usize,
    l: usize
)->F {
    // Concatenate b_bits and x_bits to form the evaluation input vector.
    // Since the vectors have already been produced in parallel, a sequential extend is sufficient here.
    let point = build_bx_point(b,x,nu,l);
    f.evaluate(&point)
}

/// Builds the assignment vector (b,x) in (nu + l) bits, also from MSB to LSB.
/// The first `nu` bits correspond to `b`, and the next `l` bits to `x`.
pub fn build_bx_point<F: Field>(b_val: usize,x_val: usize,nu: usize,l: usize)->Vec<F>{
    let mut point=Vec::with_capacity(nu+l);

    // b in MSB-first order
    for i in (0..nu).rev() {
        let bit = if ((b_val>>i)&1)==1 {
            F::one()
        } else {
            F::zero()
        };
        point.push(bit);
    }

    // x in MSB-first order
    for i in (0..l).rev() {
        let bit = if ((x_val>>i)&1)==1 {
            F::one()
        } else {
            F::zero()
        };
        point.push(bit);
    }
    point
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::UniformRand;
    use ark_bls12_381::Fr as FF;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_build_fj() {
        test_build_fj_polynomial(2,2);
        test_build_fj_polynomial(4,2);
    }

    fn test_build_fj_polynomial(n:usize, x:usize){
        let nu=(n as f64).log2() as usize;
        let l=(x as f64).log2() as usize;
        let mut rng=StdRng::seed_from_u64(99);

        let gs_for_j: Vec<MultilinearPolynomial<FF>> = (0..n).map(|_| {
            let evals=(0..x).map(|_| FF::rand(&mut rng)).collect();
            MultilinearPolynomial::new(evals)
        }).collect();

        let f_j=build_fj_polynomial(&gs_for_j);
        assert_eq!(f_j.get_num_vars(), nu+l);
        for bv in 0..n {
            for xv in 0..x {
                let expect=gs_for_j[bv].z[xv];
                let actual=evaluate_fj_at_decimals(&f_j,bv,xv,nu,l);
                assert_eq!(actual,expect);
            }
        }
    }
}
