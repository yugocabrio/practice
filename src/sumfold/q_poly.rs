use ark_ff::{Field, One, Zero};
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
pub fn build_q_polynomial<F: Field>(
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

/// a simple product F
pub fn product_F<F: Field>(vals:&[F])->F{
    let mut acc=F::one();
    for &v in vals {
        acc*=v;
    }
    acc
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{Field, UniformRand};
    use ark_bls12_381::Fr as FF;
    use crate::sumfold::multilinear::MultilinearPolynomial;
    use rand::rngs::StdRng;
    use rand::{SeedableRng, RngCore};

    #[test]
    fn test_build_q_poly() {
        test_build_q_polynomial_simple(2,2,2);
        test_build_q_polynomial_simple(8,4,2);
    }

    fn test_build_q_polynomial_simple(n:usize,x:usize,t:usize){
        let nu=(n as f64).log2() as usize;
        let l =(x as f64).log2() as usize;
        let mut rng=StdRng::seed_from_u64(111);

        
        let mut g_bj=Vec::with_capacity(n);
        for _ in 0..n {
            let mut polys_b=Vec::with_capacity(t);
            for _j in 0..t {
                let evals= (0..x).map(|_| FF::rand(&mut rng)).collect();
                polys_b.push(MultilinearPolynomial::new(evals));
            }
            g_bj.push(polys_b);
        }

        // build f_j
        let size=1<<(nu+l);
        let mut f_js=Vec::with_capacity(t);
        for j in 0..t {
            let mut f_eval=vec![FF::zero();size];
            for b_val in 0..n {
                for x_val in 0..x {
                    let idx=(b_val<<l)+x_val;
                    f_eval[idx]=g_bj[b_val][j].z[x_val];
                }
            }
            f_js.push(MultilinearPolynomial::new(f_eval));
        }

        let random_u64=rng.next_u64();
        let rho=(random_u64 as usize)%n;

        let Q= build_q_polynomial(&f_js,&product_F::<FF>,rho,nu,l);
        let q_r_b= Q.z[rho];

        // check
        let mut sum_val= FF::zero();
        for x_val in 0..x {
            let mut pv=FF::one();
            for j in 0..t {
                pv *= g_bj[rho][j].z[x_val];
            }
            sum_val += pv;
        }
        assert_eq!(q_r_b,sum_val);
    }
}
