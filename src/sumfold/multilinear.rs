use ark_ff::Field;

#[derive(Clone, Debug)]
pub struct MultilinearPolynomial<F: Field> {
    pub z: Vec<F>,
    num_vars: usize,
}

impl<F: Field> MultilinearPolynomial<F> {
    pub fn new(z: Vec<F>) -> Self {
        let len = z.len();
        let num_vars = (len as f64).log2() as usize;
        assert_eq!(1 << num_vars, len, "z length must be power-of-2");
        Self { z, num_vars }
    }

    pub fn get_num_vars(&self) -> usize {
        self.num_vars
    }

    pub fn len(&self) -> usize {
        self.z.len()
    }

    pub fn evaluate(&self, point: &[F]) -> F {
        assert_eq!(point.len(), self.num_vars);
        let mut index = 0usize;
        for &bit in point {
            let is_one = !bit.is_zero();
            index = (index << 1) | (is_one as usize);
        }
        self.z[index]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Fr as FF;

    #[test]
    fn test_mlt() {
        let evals = vec![
            FF::from(0u64),FF::from(1u64),FF::from(2u64),FF::from(3u64),
            FF::from(4u64),FF::from(5u64),FF::from(6u64),FF::from(7u64),
        ];
        let mp = MultilinearPolynomial::new(evals.clone());
        assert_eq!(mp.get_num_vars(),3);

        for i in 0..8 {
            let b2=(i>>2)&1;
            let b1=(i>>1)&1;
            let b0=i&1;
            let bits=[FF::from(b2 as u64),FF::from(b1 as u64),FF::from(b0 as u64)];
            let val=mp.evaluate(&bits);
            assert_eq!(val,FF::from(i as u64));
        }
    }
}
