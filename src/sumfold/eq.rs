use ark_ff::Field;

/// eq(b_bits, x_bits) = ∏( x_i^(b_i) * (1-x_i)^(1-b_i) )
/// b_bits, x_bits are 0 or 1
#[derive(Clone)]
pub struct EqPolynomial<F: Field> {
    b_bits: Vec<F>,
}

impl<F: Field> EqPolynomial<F> {
    pub fn new(b_bits: Vec<F>) -> Self {
        Self { b_bits }
    }

    ///  eq(b_bits, x_bits) based on x_bits
    pub fn evaluate(&self, x_bits: &[F]) -> F {
        assert_eq!(x_bits.len(), self.b_bits.len());
        let mut acc = F::one();
        for (b, x) in self.b_bits.iter().zip(x_bits.iter()) {
            if b.is_zero() {
                // b=0
                acc *= F::one() - *x;
            } else {
                // b=1
                acc *= *x;
            }
        }
        acc
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_ff::{One, Zero};
    use ark_bls12_381::Fr as F;

    #[test]
    fn test_eq_poly() {
        let eqp = EqPolynomial::new(vec![F::one(), F::zero(), F::one()]);
        for x in 0..8 {
            let b2 = (x>>2) & 1;
            let b1 = (x>>1) & 1;
            let b0 = x & 1;
            let x_bits = [F::from(b2 as u64), F::from(b1 as u64), F::from(b0 as u64)];
            let val = eqp.evaluate(&x_bits);
            if x==5 {
                assert!(val.is_one(),"x=5 => eq=1");
            } else {
                assert!(val.is_zero(),"x={} => eq=0", x);
            }
        }
    }
}
