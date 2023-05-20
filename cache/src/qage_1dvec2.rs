use crate::qk::*;
use crate::qk15::*;
use crate::qk21::*;
use crate::qk31::*;
use crate::qk41::*;
use crate::qk51::*;
use crate::qk61::*;
use crate::qsrt2::*;
use crate::qag_1dvec_integrator_result::Qag1DVecIntegratorResult;
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;
use crate::result_state::*;
use crate::qage_1dvec::*;
use crate::qk21_simd::Qk21Simd;
use crate::qk61_simd2::Qk61Simd2;
use crate::qk61_simd::Qk61Simd;


#[derive(Clone)]
pub struct Qag_1dvec2 {
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





impl Qag_1dvec2 {
    pub fn qintegrate(&self, f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> Qag1DVecIntegratorResult {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return Qag1DVecIntegratorResult::new_error(ResultState::Invalid)
        }
        //            first approximation to the integral

        let mut neval = 0;
        let mut last= 1 ;
        let mut result = 0.0;
        let mut abserr = 0.0;
        let mut defabs = 0.0;
        let mut resabs = 0.0;
        let mut result_list = vec![];


        //let qk15 = Qk15 {};
        let qk21 = Qk21Simd {};
        //let qk31 = Qk31 {};
        //let qk41 = Qk41 {};
        //let qk51 = Qk51 {};
        let qk61 = Qk61Simd2 {};

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }
        match keyf {
            //1 => (result, abserr, defabs, resabs) = qk15.integrate(f, a, b),
            2 => (result, abserr, defabs, resabs) = qk21.integrate(f, a, b),
            //3 => (result, abserr, defabs, resabs) = qk31.integrate(f, a, b),
            //4 => (result, abserr, defabs, resabs) = qk41.integrate(f, a, b),
            //5 => (result, abserr, defabs, resabs) = qk51.integrate(f, a, b),
            6 => (result, abserr, defabs, resabs) = qk61.integrate(f, a, b),
            _ => (),
        }

        result_list.push(Result::new(a,b,result,abserr));


        //           test on accuracy.

        let mut errbnd = epsabs.max(epsrel * result.abs());
        if abserr <= 50.0 * EPMACH * defabs && abserr > errbnd {
            return Qag1DVecIntegratorResult::new_error(ResultState::BadTolerance)
        }
        if (abserr <= errbnd && abserr != resabs) || abserr == 0.0 {
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
            if keyf == 1 { neval = 30 * neval + 15; }
            return Qag1DVecIntegratorResult::new(result, abserr, neval, result_list, last)
        }
        if self.limit == 1 {
            return Qag1DVecIntegratorResult::new_error(ResultState::MaxIteration)
        }

        //          initialization
        let mut area = result;
        let mut errsum = abserr;
        let mut iroff1 = 0;
        let mut iroff2 = 0;

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.


        while last <= self.limit {

            let mut new_interval = 0;

            for k in 0..last {
                if result_list[k].error > epsabs.max(epsrel * area.abs()) / last as f64 {
                    new_interval += 1;
                    let a1 = result_list[k].a;
                    let b1 = 0.5 * (result_list[k].a + result_list[k].b);
                    let a2 = b1;
                    let b2 = result_list[k].b;


                    let area1: f64;
                    let error1: f64;
                    let area2: f64;
                    let error2: f64;
                    let defab1: f64;
                    let defab2: f64;


                    match keyf {
                        //1 => {
                        //    (area1, error1, _, defab1) = qk15.integrate(f, a1, b1);
                        //    (area2, error2, _, defab2) = qk15.integrate(f, a2, b2);
                        //},
                        2 => {
                            (area1, error1, _, defab1) = qk21.integrate(f, a1, b1);
                            (area2, error2, _, defab2) = qk21.integrate(f, a2, b2);
                        },
                        //3 => {
                        //    (area1, error1, _, defab1) = qk31.integrate(f, a1, b1);
                        //    (area2, error2, _, defab2) = qk31.integrate(f, a2, b2);
                        //},
                        //4 => {
                        //    (area1, error1, _, defab1) = qk41.integrate(f, a1, b1);
                        //    (area2, error2, _, defab2) = qk41.integrate(f, a2, b2);
                        //},
                        //5 => {
                        //    (area1, error1, _, defab1) = qk51.integrate(f, a1, b1);
                        //    (area2, error2, _, defab2) = qk51.integrate(f, a2, b2);
                        //},
                        6 => {
                            (area1, error1, _, defab1) = qk61.integrate(f, a1, b1);
                            (area2, error2, _, defab2) = qk61.integrate(f, a2, b2);
                        },
                        _ => {
                            (area1, error1, _, defab1) = (0.0, 0.0, 0.0, 0.0);
                            (area2, error2, _, defab2) = (0.0, 0.0, 0.0, 0.0);
                        },
                    }


                    //           improve previous approximations to integral
                    //           and error and test for accuracy.

                    neval += 1;
                    let area12 = area1 + area2;
                    let erro12 = error1 + error2;
                    errsum += erro12 - result_list[k].error;
                    area += area12 - result_list[k].result;


                    if defab1 == error1 || defab2 == error2 {} else {
                        if (result_list[k].result - area12).abs() <= 0.00001 * area12.abs() && erro12 >= 0.99 * result_list[k].error {
                            iroff1 += 1;
                        }
                        if last > 10 && erro12 > result_list[k].error { iroff2 += 1; }
                    }
                    errbnd = epsabs.max(epsrel * area.abs());


                    if errsum > errbnd {

                        //           test for roundoff error.

                        if iroff1 >= 6 || iroff2 >= 20 {
                            return Qag1DVecIntegratorResult::new_error(ResultState::BadTolerance)
                        }

                        //           set error flag in the case that the number of subintervals
                        //           equals limit.

                        if last == self.limit {
                            return Qag1DVecIntegratorResult::new_error(ResultState::MaxIteration)
                        }

                        //           set error flag in the case of bad integrand behaviour
                        //          at a point of the integration range.

                        if a1.abs().max(b2.abs()) <= (1.0 + 100.0 * EPMACH) * (a2.abs() + 1000.0 * UFLOW) {
                            return Qag1DVecIntegratorResult::new_error(ResultState::BadFunction)
                        }

                        //           append the newly-created intervals to the list.
                    }


                    result_list[k] = Result::new(a1, b1, area1, error1);
                    result_list.push(Result::new(a2, b2, area2, error2));

                    if errsum <= errbnd {
                        break
                    }
                }
            }

            //  result_list.sort_by(|a,b| b.error.total_cmp(&a.error));
            last += new_interval;

            if errsum <= errbnd {
                break
            }
        }
        //           compute final result.

        result = 0.0;
        for k in &result_list {
            result += k.result;
        }

        abserr = errsum;


        if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
        if keyf == 1 { neval = 30 * neval + 15; }


        Qag1DVecIntegratorResult::new(result, abserr, neval, result_list, last)
    }
}


impl QuadIntegralMethod for Qag_1dvec2 {
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64) -> QuadIntegratorResult{
        QuadIntegratorResult::new_qag_vec( self.qintegrate(f,a,b,epsabs,epsrel))
    }
}


#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qage2::Qag2;
    use crate::qage::Qag;
    use crate::qage_1dvec2::Qag_1dvec2;
    use crate::qage_1dvec::Qag_1dvec;
    use crate::qage_nosort::QagNosort;
    use crate::qage_nosort_findmax::QagNosortFindmax;

    #[test]
    fn test() {
        let test: bool = true; // set to true if you want to have the results printed out

        let f = |x: f64| (1.0/x.powi(2)).sin();
        let a = 0.0;
        let b = 1.0;
        let epsabs = 5.0e-3;
        let epsrel = 0.0;
        let key = 6;
        let limit = 10000000;
        let max = 13;
        let qag = Qag { key, limit };
        let qag2 = Qag2 { key, limit };
        let qag_1dvec = Qag_1dvec { key, limit };
        let qag_1dvec2 = Qag_1dvec2 { key, limit };
        let qag_nosort = QagNosort{key,limit};
        let qag_nosort_findmax = QagNosortFindmax{key,limit};


        for _k in 0..10{

        let (mut t1, mut t2, mut t3, mut t4, mut t5, mut t6) =
            (0.0, 0.0, 0.0, 0.0, 0.0,0.0);


        for i in 0..max {
            //let start = Instant::now();
            //let res = qag.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
            //if i > 10 { t1 += start.elapsed().as_secs_f64(); }
            let start = Instant::now();
            let res3 = qag_1dvec.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
            if i > 10 { t3 += start.elapsed().as_secs_f64(); }
            //let start = Instant::now();
            //let res4 = qag_1dvec2.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
            //if i > 10 { t4 += start.elapsed().as_secs_f64(); }
            let start = Instant::now();
            let res2 = qag2.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
            if i > 10 { t2 += start.elapsed().as_secs_f64(); }
            let start = Instant::now();
            let res5 = qag_nosort.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
            if i > 10 { t5 += start.elapsed().as_secs_f64(); }
            let start = Instant::now();
            let res6 = qag_nosort_findmax.qintegrate(&f, a, b, epsabs, epsrel).unwrap();
            if i > 10 { t6 += start.elapsed().as_secs_f64(); }


            if test && i == max-1 {

                //println!("{:?}",res2.result);
                //println!("{:?}",res3.result);
                println!("{:?}",res5.result);
                println!("{:?}",res6.result);
                println!("{:?}",res5.last);
                println!("{:?}",res6.last);
            }

        }
            t1 = t1 / (max as f64 - 10.0);
            t2 = t2 / (max as f64 - 10.0);
            t3 = t3 / (max as f64 - 10.0);
            t4 = t4 / (max as f64 - 10.0);
            t5 = t5 / (max as f64 - 10.0);
            t6 = t6 / (max as f64 - 10.0);


            println!("qag : {:?} ; qag2 : {:?}; qag_1dvec: {:?}; qag_1dvec2 : {:?}; \
            qag_nosort : {:?}; qag_nosort_findmax : {:?}", t1, t2, t3, t4,t5,t6);
    }

    }

}