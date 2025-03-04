use zkp_libra::sumfold::{
  fj_poly::{build_fj_polynomial, evaluate_fj_at_decimals},
  multilinear::MultilinearPolynomial,
  prover::{SumfoldInstance, sumfold},
  q_poly::build_Q_polynomial,
};
use zkp_libra::sumcheck::SumCheckProof;
use ark_ff::{Field, PrimeField, UniformRand, One, Zero};
use ark_bls12_381::Fr as FF;
use ark_bls12_381::{Bls12_381 as E, Fr as BLSFr};
use merlin::Transcript;
use rand::{rngs::StdRng, Rng, RngCore, SeedableRng};
use std::sync::Arc;

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

fn test_build_fj_polynomial(n:usize, x:usize){
  println!("test_build_fj_polynomial n={} x={}",n,x);
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

#[allow(non_snake_case)]
fn test_build_Q_polynomial_simple(n:usize,x:usize,t:usize){
  println!("test_build_Q_polynomial_simple n={} x={} t={}",n,x,t);
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

  let Q= build_Q_polynomial(&f_js,&product_f::<FF>,rho,nu,l);
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

pub fn product_f<F: Field>(vals:&[F])->F{
  let mut acc=F::one();
  for &v in vals {
      acc*=v;
  }
  acc
}

#[test]
fn test_build_fj() {
  for _ in 0..10 {
    test_build_fj_polynomial(2,2);
    test_build_fj_polynomial(4,2);
    test_build_fj_polynomial(8,2);
    test_build_fj_polynomial(16,2);
    test_build_fj_polynomial(32,2);
  }
}

#[allow(non_snake_case)]
#[test]
fn test_build_Q_poly() {
  for _ in 0..10 {
    test_build_Q_polynomial_simple(2,2,2);
    test_build_Q_polynomial_simple(8,4,2);
    test_build_Q_polynomial_simple(8,8,2);
    test_build_Q_polynomial_simple(16,8,2);
  }
}

#[test]
fn test_sumfold_correctness() {
  // n: num of instances
  let ns = [2, 4, 8, 16];
  // l: num of bits for g polynomial
  let ls = [2, 4, 8, 16];
  let num_tries = 10;

  for &n in &ns {
    for &l in &ls {
      for t in 0..num_tries {
        println!("test_sumfold n={} l={} t={}",n,l,t);
        let mut rng=StdRng::seed_from_u64(123);

        // define F as product
        let f_arc: Arc<dyn Fn(&[BLSFr])->BLSFr + Send + Sync> =
            Arc::new(|vals: &[BLSFr]| vals.iter().product());

        // build 2 instances
        let mut instances = Vec::with_capacity(n);
        for _ in 0..n {
            // build g0,g1
            let g0= build_random_poly(l, &mut rng);
            let g1= build_random_poly(l, &mut rng);
            // store
            let inst = SumfoldInstance {
                F_func: f_arc.clone(),
                g_vec: vec![g0, g1],
            };
            instances.push(inst);
        }

        // call sumfold
        let (chosen_inst, rho_field, q_b) = sumfold::<E>(instances);

        // 1) compute the actual sum_x F(g0(x), g1(x)) for chosen_inst
        let size = chosen_inst.g_vec[0].len();
        let mut t_val = BLSFr::zero();
        for i in 0..size {
            let val = (chosen_inst.F_func)(&[
                chosen_inst.g_vec[0].z[i],
                chosen_inst.g_vec[1].z[i],
            ]);
            t_val += val;
        }

        // 2) check Q(rho) == that sum
        let rho_usize = (rho_field.into_repr().as_ref()[0]) as usize;
        let qb_rho = q_b.z[rho_usize];
        assert_eq!(
            qb_rho, t_val,
            "q_b(rho) must match sum_x of F(g_vec)"
        );

        // 3) naive sumcheck of Q(b):
        let total_sum: BLSFr = q_b.z.iter().copied().sum();

        let mut q_b_clone = q_b.clone();

        let mut transcript = Transcript::new(b"test_sumcheck");
        let (proof, _challenges)= SumCheckProof::<E>::prove(&mut q_b_clone, total_sum, &mut transcript);

        let mut verify_transcript = Transcript::new(b"test_sumcheck");
        let ok = proof.verify(total_sum, &mut verify_transcript);
        assert!(ok, "Naive sumcheck on q_b should pass as well");
      }
    }
  }
}

#[test]
fn test_sumfold_soundness() {
  // n: num of instances
  let ns = [2, 4, 8, 16];
  // l: num of bits for g polynomial
  let ls = [2, 4, 8, 16];
  let num_tries = 10;

  for &n in &ns {
    for &l in &ls {
      for t in 0..num_tries {
        println!("test_sumfold n={} l={} t={}",n,l,t);
        let mut rng=StdRng::seed_from_u64(123);

        // define F as product
        let f_arc: Arc<dyn Fn(&[BLSFr])->BLSFr + Send + Sync> =
            Arc::new(|vals: &[BLSFr]| vals.iter().product());

        // build 2 instances
        let mut instances = Vec::with_capacity(n);
        for _ in 0..n {
            // build g0,g1
            let g0= build_random_poly(l, &mut rng);
            let g1= build_random_poly(l, &mut rng);
            // store
            let inst = SumfoldInstance {
                F_func: f_arc.clone(),
                g_vec: vec![g0, g1],
            };
            instances.push(inst);
        }

        // call sumfold
        let (chosen_inst, rho_field, q_b) = sumfold::<E>(instances);

        // 1) compute the actual sum_x F(g0(x), g1(x)) for chosen_inst
        let size = chosen_inst.g_vec[0].len();
        let mut t_val = BLSFr::zero();
        for i in 0..size {
            let val = (chosen_inst.F_func)(&[
                chosen_inst.g_vec[0].z[i],
                chosen_inst.g_vec[1].z[i],
            ]);
            t_val += val;
        }

        // 2) check Q(rho) != invalid sum
        let rho_usize = (rho_field.into_repr().as_ref()[0]) as usize;
        let qb_rho = q_b.z[rho_usize];
        for sigma in 1..11 {
          assert!(
            qb_rho != (t_val + BLSFr::from(sigma as u64)),
            "q_b(rho) must not match invalid sum"
          );

          assert!(
            qb_rho != (t_val - BLSFr::from(sigma as u64)),
            "q_b(rho) must not match invalid sum"
          );
        }

        // 3) naive sumcheck of Q(b):
        let total_sum: BLSFr = q_b.z.iter().copied().sum();

        let mut q_b_clone = q_b.clone();

        let mut transcript = Transcript::new(b"test_sumcheck");
        let (proof, _challenges)= SumCheckProof::<E>::prove(&mut q_b_clone, total_sum, &mut transcript);

        let mut verify_transcript = Transcript::new(b"test_sumcheck");
        let ok = proof.verify(BLSFr::zero(), &mut verify_transcript);
        assert!(!ok, "Should not pass with invalid sum");
      }
    }
  }
}