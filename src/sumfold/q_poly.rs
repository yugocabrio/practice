use ark_ff::Field;
use crate::sumfold::eq::EqPolynomial;
use crate::sumfold::fj_poly::build_bx_point;
use crate::sumfold::multilinear::MultilinearPolynomial;

/// Converts a decimal integer `val` into a bit vector of length `num_bits`,
/// from the most significant bit to the least significant bit.
/// E.g., if val=6 (0b110) and num_bits=3, this returns [1, 1, 0].
fn decimal_to_bits_msb_first<F: Field>(val: usize,num_bits: usize)->Vec<F>{
    let mut bits=Vec::with_capacity(num_bits);
    for i in (0..num_bits).rev(){
        let b= if ((val>>i)&1)==1 {
            F::one()
        }else{
            F::zero()
        };
        bits.push(b);
    }
    bits
}

/// Constructs Q(b) = eq(rho, b) * ( sum_{x in {0,1}^l} F( f_1(b,x), ..., f_t(b,x) ) ).
///
/// - `f_js`: [f_1, ..., f_t], each an MLE in (nu + l) variables (bits of b + bits of x).
/// - `F`: A function that takes t field elements and returns a single field element (e.g., product).
/// - `rho`: decimal integer in [0..2^nu).
/// - `nu`: number of bits for b.
/// - `l`: number of bits for x.
///
/// Returns a MultilinearPolynomial in `nu` variables (the variable b).
#[allow(non_snake_case)]
pub fn build_Q_polynomial<F: Field>(
    f_js:&[MultilinearPolynomial<F>],
    F: &(dyn Fn(&[F])->F + Sync),
    rho: usize,
    nu: usize,
    l: usize
)-> MultilinearPolynomial<F> {
    // 1) Build eq(rho, ·) as an EqPolynomial using MSB-first bits of rho.
    let rho_bits=decimal_to_bits_msb_first::<F>(rho,nu);
    let eq_rho=EqPolynomial::new(rho_bits);

    // 2) We'll build Q in dense form for b in [0..2^nu].
    let len_n=1<<nu;
    let len_x=1<<l;
    let t=f_js.len();
    let mut q_evals= vec![F::zero();len_n];

    // 3) For each b in [0..2^nu), compute eq(rho,b) and sum over x.
    //    Q(b) = eq(rho,b) * Σ_x F( f_1(b,x), ..., f_t(b,x) ).
    for (b_val,q_out) in q_evals.iter_mut().enumerate() {
        // Evaluate eq(rho, b_val)
        let b_bits= decimal_to_bits_msb_first::<F>(b_val,nu);
        let eq_val= eq_rho.evaluate(&b_bits);
        if eq_val.is_zero() {
            // Q(b_val) = 0 if eq(rho,b_val)=0
            *q_out = F::zero();
        } else {
            let mut sum_val= F::zero();
            for x_val in 0..len_x {
                // Evaluate each f_j(b_val,x_val)
                let bx_point=build_bx_point::<F>(b_val,x_val,nu,l);
                let mut vals=Vec::with_capacity(t);
                for f_j in f_js {
                    vals.push(f_j.evaluate(&bx_point));
                }
                // Apply F to these t values
                sum_val += F(&vals);
            }
            // Multiply by eq(rho,b_val)
            *q_out = eq_val * sum_val;
        }
    }

    // 4) Return Q(b) as a multilinear polynomial in nu variables.
    MultilinearPolynomial::new(q_evals)
}