use std::collections::{BinaryHeap, HashMap};
use crate::constants::*;
use crate::qag_vec_norm_integrator_result::QagVecNormIntegratorResult;
use crate::qk15::qk15_quadrature;
use crate::qk21::qk21_quadrature;
use crate::qk31::qk31_quadrature;
use crate::qk41::qk41_quadrature;
use crate::qk51::qk51_quadrature;
use crate::qk61::qk61_quadrature;
use crate::qk61_vec_norm2::Qk61VecNorm2;
use crate::result_state::*;
use crate::semi_infinite_function::{DoubleInfiniteFunction, SemiInfiniteFunction};


#[derive(Clone)]
pub struct Qag {
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





impl Qag {

    pub fn integrate<const N:usize,F>(&self, f : &F, a : f64, b : f64, epsabs : f64, epsrel : f64)
                                       ->  QagVecNormIntegratorResult<N>
        where F : Fn(f64) -> [f64;N] {
        if b == f64::INFINITY && a.is_finite() {
            let f2 = |x: f64| SemiInfiniteFunction(&f, x, a, b);
            return self.qintegrate(&f2, 0.0, 1.0, epsabs, epsrel)
        }
        if a == f64::NEG_INFINITY && b.is_finite() {
            let f2 = |x: f64| SemiInfiniteFunction(&f, x, b, a);
            return self.qintegrate(&f2, 0.0, 1.0, epsabs, epsrel)
        }
        if a == f64::NEG_INFINITY && b == f64::INFINITY {
            let f2 = |x: f64| DoubleInfiniteFunction(&f, x);
            return self.qintegrate(&f2, -1.0, 1.0, epsabs, epsrel)
        }




        self.qintegrate(&f,a,b,epsabs,epsrel)
    }


    pub fn qintegrate<const N:usize,F>(&self, f : &F, a : f64, b : f64, epsabs : f64, epsrel : f64)
                                     ->  QagVecNormIntegratorResult<N>
    where F : Fn(f64) -> [f64;N]{

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagVecNormIntegratorResult::new_error(ResultState::Invalid)
        }

        let mut neval = 0;
        let mut last= 1 ;

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }

        let (mut result, mut abserr, mut rounderr) = match keyf {
            1 => qk15_quadrature(f, a, b),
            2 => qk21_quadrature(f, a, b),
            3 => qk31_quadrature(f, a, b),
            4 => qk41_quadrature(f, a, b),
            5 => qk51_quadrature(f, a, b),
            6 => qk61_quadrature(f, a, b),
            _ => ([0.0; N], 0.0, 0.0),
        };


        let mut errbnd = epsabs.max(epsrel * norm_vec(&result));

        let mut interval_cache = HashMap::from([((Myf64{x:a},Myf64{x:b}),result.clone())]);
        let mut heap = BinaryHeap::new();
        heap.push(HeapItem::new((a,b),abserr));


        if abserr + rounderr <= errbnd{
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * last as i32 - 1); }
            if keyf == 1 { neval = 30 * last as i32 + 15; }
            abserr = abserr + rounderr;
            return QagVecNormIntegratorResult::new(result,abserr,neval,last)
        }

        if self.limit == 1 {
            return QagVecNormIntegratorResult::new_error(ResultState::MaxIteration)
        }

        if abserr < rounderr {
            return QagVecNormIntegratorResult::new_error(ResultState::BadTolerance)
        }


        while last < self.limit{
            let mut to_process = vec![];
            let mut err_sum = 0.0;


            while to_process.len() < 128 && heap.len() != 0{
                let old_interval = heap.pop().unwrap();
                let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                let old_res = interval_cache.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                err_sum += old_err;
                to_process.push((x,y,old_err,old_res));
                if err_sum > abserr - errbnd / 8.0 { break}
            }


            for comp in to_process{
                last += 1;

                let mut result1 = [0.0; N];
                let mut abserr1 = 0.0;
                let mut rounderr1  = 0.0;

                let mut result2 = [0.0; N];
                let mut abserr2 = 0.0;
                let mut rounderr2  = 0.0;

                let a1 = comp.0;
                let b1 = 0.5 * (comp.0 + comp.1);
                let a2 = b1;
                let b2 = comp.1;


                match keyf {
                    1 => {
                        (result1, abserr1, rounderr1) = qk15_quadrature(f, a1, b1);
                        (result2, abserr2, rounderr2) = qk15_quadrature(f, a2, b2);
                    },
                    2 => {
                        (result1, abserr1, rounderr1) = qk21_quadrature(f, a1, b1);
                        (result2, abserr2, rounderr2) = qk21_quadrature(f, a2, b2);
                    },
                    3 => {
                        (result1, abserr1, rounderr1) = qk31_quadrature(f, a1, b1);
                        (result2, abserr2, rounderr2) = qk31_quadrature(f, a2, b2);
                    },
                    4 => {
                        (result1, abserr1, rounderr1) = qk41_quadrature(f, a1, b1);
                        (result2, abserr2, rounderr2) = qk41_quadrature(f, a2, b2);
                    },
                    5 => {
                        (result1, abserr1, rounderr1) = qk51_quadrature(f, a1, b1);
                        (result2, abserr2, rounderr2) = qk51_quadrature(f, a2, b2);
                    },
                    6 => {
                        (result1, abserr1, rounderr1) = qk61_quadrature(f, a1, b1);
                        (result2, abserr2, rounderr2) = qk61_quadrature(f, a2, b2);
                    },
                    _ => (),
                }

                res_update(&mut result, & result1, & result2, &comp.3);

                interval_cache.insert((Myf64{x:a1},Myf64{x:b1}),result1);
                interval_cache.insert((Myf64{x:a2},Myf64{x:b2}),result2);

                heap.push(HeapItem::new((a1,b1),abserr1));
                heap.push(HeapItem::new((a2,b2),abserr2));


                abserr += -comp.2 + abserr1 + abserr2;
                rounderr += rounderr1 + rounderr2;

            }
            errbnd = epsabs.max(epsrel * norm_vec(&result));


            if abserr <= errbnd / 8.0{ break;}
            if abserr < rounderr {
                return QagVecNormIntegratorResult::new_error(ResultState::BadTolerance)
            }
        }


        if last >= self.limit {
            return QagVecNormIntegratorResult::new_error(ResultState::MaxIteration)
        }

        if keyf != 1 { neval = (10 * keyf + 1) * (2 * last as i32 - 1); }
        if keyf == 1 { neval = 30 * last as i32 + 15; }

        abserr = abserr + rounderr;

        return QagVecNormIntegratorResult::new(result,abserr,neval,last)

    }
}



#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qag::Qag;
    use crate::qage_vec_norm3::QagVecNorm3;

    #[test]
    fn test(){

        let a = 0.0;
        let b = 1.0e7;
        let key = 6;
        let limit = 1000000;
        let epsrel = 0.0;
        let epsabs = 1.0;
        let max = 25;

        let f = |x:f64| [x.cos(),x.sin(),x.cos(),x.sin()];
        let qag3 = QagVecNorm3{key,limit};
        let qag = Qag{key,limit};

        let mut res;
        let mut res2;

        for k in 0..max{
            let start = Instant::now();
            res = qag3.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
            println!("qk struct time : {:?}",start.elapsed());
            let start = Instant::now();
            res2 = qag.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
            println!("qk nostruct time : {:?}",start.elapsed());

            if k == max-1{
                println!("{:?}",res);
                println!("{:?}",res2);
            }

        }
    }


    #[test]
    fn keys() {

        let a = 0.0;
        let b = 1.0;
        let limit = 1000000;
        let epsrel = 0.0;
        let epsabs = 7e-3;
        let max = 2;

        let f = |x:f64| [(1.0/(x*x)).cos(),(1.0/(x*x)).sin(),(1.0/(x*x)).cos(),(1.0/(x*x)).sin()];
        let qag1 = Qag{key : 1,limit};
        let qag2 = Qag{key : 2,limit};
        let qag3 = Qag{key : 3,limit};
        let qag4 = Qag{key : 4,limit};
        let qag5 = Qag{key : 5,limit};
        let qag6 = Qag{key : 6,limit};


        let mut res;
        let mut res2;
        let mut res3;
        let mut res4;
        let mut res5;
        let mut res6;


        for k in 0..max{
            let start = Instant::now();
            res = qag1.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
            println!("qk 15 : {:?}",start.elapsed());
            let start = Instant::now();
            res2 = qag2.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
            println!("qk 21 : {:?}",start.elapsed());
            let start = Instant::now();
            res3 = qag3.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
            println!("qk 31 : {:?}",start.elapsed());
            let start = Instant::now();
            res4 = qag4.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
            println!("qk 41 : {:?}",start.elapsed());
            let start = Instant::now();
            res5 = qag5.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
            println!("qk 51 : {:?}",start.elapsed());
            let start = Instant::now();
            res6 = qag6.qintegrate(&f,a,b,epsabs,epsrel).unwrap();
            println!("qk 61 : {:?}",start.elapsed());

            if k == max-1{
                println!("{:?}",res);
                println!("{:?}",res2);
                println!("{:?}",res3);
                println!("{:?}",res4);
                println!("{:?}",res5);
                println!("{:?}",res6);
            }

        }


    }


    #[test]
    fn SemiInfinite() {
        let a = 0.0;
        let b = f64::INFINITY;
        let limit = 1000000;
        let epsrel = 0.0;
        let epsabs = 1e-2;

        let f = |x:f64| [(1.0/(x*x)).sin()];
        let qag = Qag{key : 1,limit};

        let res = qag.integrate(&f,a,b,epsabs,epsrel);
        println!("{:?}",res.unwrap());
    }


    #[test]
    fn DoubleInfinite() {
        let a = f64::NEG_INFINITY;
        let b = f64::INFINITY;
        let limit = 1000000;
        let epsrel = 0.0;
        let epsabs = 1e-2;

        let f = |x:f64| [(1.0/(x*x)).sin()];
        let qag = Qag{key : 1,limit};

        let res = qag.integrate(&f,a,b,epsabs,epsrel);
        println!("{:?}",res.unwrap());
    }



}


