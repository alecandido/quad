use::rayon::prelude::*;

use std::collections::{BinaryHeap, HashMap};
use crate::constants::*;
use crate::qag_par_integrator_result::QagParIntegratorResult;
use crate::qk15::qk15_array_quadrature;
use crate::qk21::qk21_array_quadrature;
use crate::qk31::qk31_array_quadrature;
use crate::qk41::qk41_array_quadrature;
use crate::qk51::qk51_array_quadrature;
use crate::qk61::qk61_array_quadrature;
use crate::result_state::*;




#[derive(Clone)]
pub struct QagParIter {
    pub key : i32,
    pub limit : usize,
}

///           f      : f64
///                     function
///
///           a      : f64
///                    lower limit of integration
///
///           b      : f64
///                    upper limit of integration
///
///           epsabs : f64
///                    absolute accuracy requested
///
///           epsrel : f64
///                    relative accuracy requested
///                    if  epsabs <= 0 && epsrel <= max(50*rel.mach.acc.,0.5d-28),
///                    the fn will return with result_state = Invalid.
///
///            key   : i32
///                    key for choice of local integration rule. A gauss-kronrod pair is used with:
///                          7 - 15 points if key < 2,
///                         10 - 21 points if key = 2,
///                         15 - 31 points if key = 3,
///                         20 - 41 points if key = 4,
///                         25 - 51 points if key = 5,
///                         30 - 61 points if key > 5.
///
///            limit : i32
///                    gives an upperbound on the number of subintervals in the partition
///                    of (a,b), limit >= 1.
///
///
///
///         On return : QagIntegratorResult :
///
///           QagIntegrationResult:
///           result : f64
///                    Approximation to the integral.
///
///           abserr : f64
///                    Estimate of the modulus of the absolute error,
///                    which should equal or exceed abs(i-result).
///
///           neval  : i32
///                    Number of integrand evaluations.
///
///           alist  : Vec<f64>
///                      Vector of dimension at least limit, the elements of which are the left
///                      end points of the subintervals in the partition of the given integration
///                      range (a,b).
///
///           blist  : Vec<f64>
///                      Vector of dimension at least limit, the elements of which are the right
///                      end points of the subintervals in the partition of the given integration
///                      range (a,b).
///
///           rlist  : Vec<f64>
///                      Vector of dimension at least limit, the elements of which are the integral
///                      approximations on the subintervals.
///
///            rlist  : Vec<f64>
///                      Vector of dimension at least limit, the elements of which are the moduli
///                      of the absolute error estimates on the subintervals.
///
///            iord   : Vec<usize>
///                      Vector of dimension at least limit, the elements of which are pointers to
///                      the error estimates over the subintervals, such that
///                      elist(iord(1)), ...,elist(iord(k)) form a decreasing sequence,
///                      with k = last if last <= (limit/2+2), and
///                      k = limit+1-last otherwise.
///
///            last    : usize
///                      number of subintervals actually produced in the
///                      subdivision process
///
///
///
///
///           ResultState =
///           Success :
///                    Normal and reliable termination of the routine. it is assumed that the
///                    requested accuracy has been achieved.
///           MaxIteration :
///                    The maximum number of steps has been executed. the integral is probably too
///                    difficult to be calculated by dqng.
///           Invalid :
///                     The input is invalid, because epsabs <= 0 &&
///                     epsrel < max(50 * rel.mach.acc.,0.5e-28).
///           BadTolerance :
///                     The occurrence of roundoff error is detected, which prevents the requested
///                     tolerance from being achieved.
///           BadFunction :
///                     Extremely bad integrand behaviour occurs at some points of the integration
///                     interval.
///
///
///           If ResultState != Succes =>    It is assumed that the requested accuracy has not
///           been achieved.
///
///
///





impl QagParIter {
    pub fn qintegrate<const N:usize>(&self, fun : &FnVecGen<N>, a : f64, b : f64, epsabs : f64, epsrel : f64)
                                     ->  QagParIntegratorResult<N> {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagParIntegratorResult::new_error(ResultState::Invalid)
        }

        let f = &fun.components;


        //            first approximation to the integral

        let mut neval = 0;
        let mut last= 1 ;
        let mut result = [0.0; N];
        let mut abserr = 0.0;
        let mut rounderr  =  0.0;

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }
        match keyf {
            1 => (result, abserr, rounderr) = qk15_array_quadrature(&**f, a, b),
            2 => (result, abserr, rounderr) = qk21_array_quadrature(&**f, a, b),
            3 => (result, abserr, rounderr) = qk31_array_quadrature(&**f, a, b),
            4 => (result, abserr, rounderr) = qk41_array_quadrature(&**f, a, b),
            5 => (result, abserr, rounderr) = qk51_array_quadrature(&**f, a, b),
            6 => (result, abserr, rounderr) = qk61_array_quadrature(&**f, a, b),
            _ => (),
        }

        //           test on accuracy.

        let mut errbnd = epsabs.max(epsrel * norm_vec(&result));

        let mut interval_cache = HashMap::from([((Myf64{x:a}, Myf64{x:b}), result.clone())]);
        let mut heap = BinaryHeap::new();
        heap.push(HeapItem::new((a,b),abserr));

        if abserr < rounderr {
            return QagParIntegratorResult::new_error(ResultState::BadTolerance)
        }


        if abserr <= errbnd{
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * last as i32 - 1); }
            if keyf == 1 { neval = 30 * last as i32 + 15; }
            return QagParIntegratorResult::new(result, abserr, neval, last)
        }

        if self.limit == 1 {
            return QagParIntegratorResult::new_error(ResultState::MaxIteration)
        }

        //          initialization

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.




        while last < self.limit {
            let mut to_process = vec![];
            let mut err_sum = 0.0;


            while to_process.len() < 128 && heap.len() != 0 {
                let old_interval = heap.pop().unwrap();
                let ((x, y), old_err) = (old_interval.interval, old_interval.err);
                let old_res = interval_cache.remove(&(Myf64 { x }, Myf64 { x: y })).unwrap();
                err_sum += old_err;
                sub_vec(&mut result, &old_res);
                to_process.push((x, y));
                if err_sum > abserr - errbnd / 8.0 { break }
            }

            abserr -= err_sum;
            last += to_process.len();

            let new_result : (Vec<_>,Vec<_>) = to_process.par_iter().map(|comp| {
                let mut result1 = [0.0; N];
                let mut abserr1 = 0.0;
                let mut rounderr1 = 0.0;

                let mut result2 = [0.0; N];
                let mut abserr2 = 0.0;
                let mut rounderr2 = 0.0;

                let a1 = comp.0;
                let b1 = 0.5 * (comp.0 + comp.1);
                let a2 = b1;
                let b2 = comp.1;

                match keyf {
                    1 => {
                        (result1, abserr1, rounderr1) = qk15_array_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk15_array_quadrature(&**f, a2, b2);
                    },
                    2 => {
                        (result1, abserr1, rounderr1) = qk21_array_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk21_array_quadrature(&**f, a2, b2);
                    },
                    3 => {
                        (result1, abserr1, rounderr1) = qk31_array_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk31_array_quadrature(&**f, a2, b2);
                    },
                    4 => {
                        (result1, abserr1, rounderr1) = qk41_array_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk41_array_quadrature(&**f, a2, b2);
                    },
                    5 => {
                        (result1, abserr1, rounderr1) = qk51_array_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk51_array_quadrature(&**f, a2, b2);
                    },
                    6 => {
                        (result1, abserr1, rounderr1) = qk61_array_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk61_array_quadrature(&**f, a2, b2);
                    },
                    _ => (),
                }
                ((a1, b1, result1, abserr1, rounderr1), (a2, b2, result2, abserr2, rounderr2))
            }
            ).collect();

            for k in 0..new_result.0.len() {
                add_vec(&mut result, &new_result.0[k].2);
                add_vec(&mut result, &new_result.1[k].2);
                abserr += new_result.0[k].3 + new_result.1[k].3;
                rounderr += new_result.0[k].4 + new_result.1[k].4;
                interval_cache.insert((Myf64 { x: new_result.0[k].0 }, Myf64 { x: new_result.0[k].1 }), new_result.0[k].2);
                interval_cache.insert((Myf64 { x: new_result.1[k].0 }, Myf64 { x: new_result.1[k].1 }), new_result.1[k].2);
                heap.push(HeapItem::new((new_result.0[k].0, new_result.0[k].1), new_result.0[k].3));
                heap.push(HeapItem::new((new_result.1[k].0, new_result.1[k].1), new_result.1[k].3));
            }




            errbnd = epsabs.max(epsrel * norm_vec(&result));


            if abserr <= errbnd / 8.0{ break;}
            if abserr < rounderr {
                return QagParIntegratorResult::new_error(ResultState::BadTolerance)
            }

            if last >= self.limit {
                return QagParIntegratorResult::new_error(ResultState::MaxIteration)
            }

        }



        if keyf != 1 { neval = (10 * keyf + 1) * (2 * last as i32 - 1); }
        if keyf == 1 { neval = 30 * last as i32 + 15; }


        return QagParIntegratorResult::new(result, abserr, neval, last)

    }
}

fn sub_vec( a : &mut [f64], b : &[f64] ){
    for k in 0..a.len(){
        a[k] -= b[k];
    }
}

fn add_vec( a : &mut [f64], b : &[f64] ){
    for k in 0..a.len(){
        a[k] += b[k];
    }
}


#[cfg(test)]
mod tests {
    use std::sync::Arc;
    use std::time::Instant;
    use crate::constants::FnVecGen;
    use crate::qag_array::QagArray;
    use crate::qag_par_iter::QagParIter;

    #[test]
    fn test() {
        let a = 0.0;
        let b = 5000000.0;
        let epsrel = 0.0;
        let epsabs = 1.0e-2;
        let limit = 100000;
        let key = 6;
        let max = 30;

        let qag1 = QagArray {key,limit, points: vec![], more_info: false };
        let qag2 = QagParIter {key,limit};

        let f1 = |x:f64| [x.cos(),x.sin()];
        let f = FnVecGen{ components : Arc::new(f1.clone())};

        let mut res1;
        let mut res2;

        let (mut t1,mut t2) = (0.0,0.0);


        for k in 0..max{
            let start = Instant::now();
            res1 = qag1.qintegrate(&f1,a,b,epsabs,epsrel);
            if k > 10 { t1 += start.elapsed().as_secs_f64();}
            let start = Instant::now();
            res2 = qag2.qintegrate(&f,a,b,epsabs,epsrel);
            if k > 10 { t2 += start.elapsed().as_secs_f64();}

            if k == max-1{
                println!("{:?}",res1);
                println!("{:?}",res2);
            }
        }

        t1 /= max as f64 - 10.0;
        t2 /= max as f64 - 10.0;

        println!("no parallel time : {t1} ; parallel time : {t2}");



    }
}

