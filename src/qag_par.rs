use::rayon::prelude::*;

use std::collections::{BinaryHeap, HashMap};
use crate::constants::*;
use crate::qag_integrator_result::QagIntegratorResult;
use crate::qk15::qk15_quadrature;
use crate::qk21::qk21_quadrature;
use crate::qk31::qk31_quadrature;
use crate::qk41::qk41_quadrature;
use crate::qk51::qk51_quadrature;
use crate::qk61::qk61_quadrature;
use crate::result_state::*;


#[derive(Clone)]
pub struct QagPar {
    pub key : i32,
    pub limit : usize,
    pub points : Vec<f64>,
    pub more_info : bool
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





impl QagPar {
/*
    pub fn integrate<F>(&self, fun : &FnVec, a : f64, b : f64, epsabs : f64, epsrel : f64)
                        -> QagIntegratorResult
        where F : Fn(f64) -> Vec<f64> + ?Sized {
        if b == f64::INFINITY && a.is_finite() {
            let f = &fun.components;
            let f3 =  |x: f64| semi_infinite_function(&**f, x, a, b);
            let f2 = FnVec{ components : Arc::new(f3.clone())};
            return self.qintegrate(&f2, 0.0, 1.0, epsabs, epsrel)
        }
        if a == f64::NEG_INFINITY && b.is_finite() {
            let f = &fun.components;
            let f3 =  |x: f64| semi_infinite_function(&**f, x, b, a);
            let f2 = FnVec{ components : Arc::new(f3.clone())};
            return self.qintegrate(&f2, 0.0, 1.0, epsabs, epsrel)
        }
        if a == f64::NEG_INFINITY && b == f64::INFINITY {
            let f = &fun.components;
            let f2 = FnVec{ components : Arc::new(|x: f64| double_infinite_function(&**f, x))};
            return self.qintegrate(&f2, -1.0, 1.0, epsabs, epsrel)
        }

        self.qintegrate(&fun,a,b,epsabs,epsrel)
    }

 */




    pub fn qintegrate(&self, fun : &FnVec, a : f64, b : f64, epsabs : f64, epsrel : f64)
                                     ->  QagIntegratorResult {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagIntegratorResult::new_error(ResultState::Invalid)
        }

        let f = &fun.components;
        let n: usize = f(0.0).len();

        let mut initial_intervals = vec![];

        if self.points.is_empty() { initial_intervals.push((a,b)); }
        else {
            let mut prev = a ;
            for p in &self.points{
                initial_intervals.push((prev,*p));
                prev = *p;
            }
            initial_intervals.push((prev,b));
        }




        let mut neval = 0;
        let mut last= 1 ;

        let mut interval_cache = HashMap::new();
        let mut heap = BinaryHeap::new();
        let mut result = vec![0.0; n];
        let mut abserr = 0.0;
        let mut rounderr  =  0.0;

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }




        for comp in initial_intervals{
            let (result_temp, abserr_temp, rounderr_temp) = match keyf {
                1 => qk15_quadrature(&**f, comp.0, comp.1),
                2 => qk21_quadrature(&**f, comp.0, comp.1),
                3 => qk31_quadrature(&**f, comp.0, comp.1),
                4 => qk41_quadrature(&**f, comp.0, comp.1),
                5 => qk51_quadrature(&**f, comp.0, comp.1),
                6 => qk61_quadrature(&**f, comp.0, comp.1),
                _ => (vec![0.0; n], 0.0, 0.0),
            };
            add_res(&mut result,&result_temp);
            abserr += abserr_temp;
            rounderr += rounderr_temp;
            heap.push(HeapItem::new((comp.0,comp.1),abserr_temp));
            interval_cache.insert((Myf64{x:comp.0},Myf64{x:comp.1}),result_temp);
        }

        //           test on accuracy.

        let mut errbnd = epsabs.max(epsrel * norm_vec(&result));



        if abserr + rounderr <= errbnd{
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * last as i32 - 1); }
            if keyf == 1 { neval = 30 * last as i32 + 15; }
            abserr = abserr + rounderr;
            if self.more_info { return QagIntegratorResult::new_more_info(result, abserr, neval, last, interval_cache, heap) }
            else { return QagIntegratorResult::new(result, abserr) }
        }

        if self.limit == 1 {
            return QagIntegratorResult::new_error(ResultState::MaxIteration)
        }

        if abserr < rounderr {
            return QagIntegratorResult::new_error(ResultState::BadTolerance)
        }



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

                let mut result1 = vec![0.0; n];
                let mut abserr1 = 0.0;
                let mut rounderr1 = 0.0;

                let mut result2 = vec![0.0; n];
                let mut abserr2 = 0.0;
                let mut rounderr2 = 0.0;

                let a1 = comp.0;
                let b1 = 0.5 * (comp.0 + comp.1);
                let a2 = b1;
                let b2 = comp.1;

                match keyf {
                    1 => {
                        (result1, abserr1, rounderr1) = qk15_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk15_quadrature(&**f, a2, b2);
                    },
                    2 => {
                        (result1, abserr1, rounderr1) = qk21_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk21_quadrature(&**f, a2, b2);
                    },
                    3 => {
                        (result1, abserr1, rounderr1) = qk31_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk31_quadrature(&**f, a2, b2);
                    },
                    4 => {
                        (result1, abserr1, rounderr1) = qk41_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk41_quadrature(&**f, a2, b2);
                    },
                    5 => {
                        (result1, abserr1, rounderr1) = qk51_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk51_quadrature(&**f, a2, b2);
                    },
                    6 => {
                        (result1, abserr1, rounderr1) = qk61_quadrature(&**f, a1, b1);
                        (result2, abserr2, rounderr2) = qk61_quadrature(&**f, a2, b2);
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
                interval_cache.insert((Myf64 { x: new_result.0[k].0 }, Myf64 { x: new_result.0[k].1 }), new_result.0[k].2.clone());
                interval_cache.insert((Myf64 { x: new_result.1[k].0 }, Myf64 { x: new_result.1[k].1 }), new_result.1[k].2.clone());
                heap.push(HeapItem::new((new_result.0[k].0, new_result.0[k].1), new_result.0[k].3));
                heap.push(HeapItem::new((new_result.1[k].0, new_result.1[k].1), new_result.1[k].3));
            }




            errbnd = epsabs.max(epsrel * norm_vec(&result));


            if abserr <= errbnd / 8.0{ break;}
            if abserr < rounderr {
                return QagIntegratorResult::new_error(ResultState::BadTolerance)
            }

            if last >= self.limit {
                return QagIntegratorResult::new_error(ResultState::MaxIteration)
            }

        }



        if keyf != 1 { neval = (10 * keyf + 1) * (2 * last as i32 - 1); }
        if keyf == 1 { neval = 30 * last as i32 + 15; }

        abserr = abserr + rounderr;

        if self.more_info { return QagIntegratorResult::new_more_info(result, abserr, neval, last, interval_cache, heap) }
        else { return QagIntegratorResult::new(result, abserr) }

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
    use crate::constants::FnVec;
    use crate::qag::Qag;
    use crate::qag_par::QagPar;

    #[test]
    fn test() {
        let a = 0.0;
        let b = 10000.0;
        let epsrel = 0.0;
        let epsabs = 1.0e-2;
        let limit = 10000000;
        let key = 6;
        let max = 30;

        let qag1 = Qag{key,limit, points: vec![], more_info: false };
        let qag2 = QagPar {key,limit, points: vec![], more_info: false};

        let f1 = |x:f64| vec![x.cos(),x.sin()];
        let f = FnVec{ components : Arc::new(f1.clone())};

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




