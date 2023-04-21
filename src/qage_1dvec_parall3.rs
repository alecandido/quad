use std::sync::{Arc, Mutex};
use na::max;
use rayon::{current_thread_index, ThreadPool};
use crate::funct_vector::{FnVecP, FnVecPa};
use crate::qk::*;
use crate::qsrt2::*;
use crate::qag_1dvec_integrator_result::Qag1DVecIntegratorResult;
use crate::qag_1dvec_parall_integration_result::Qag1DVecParIntegrationResult;
use crate::qag_1dvec_parall_integrator_result::Qag1DVecParIntegratorResult;
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;
use crate::result_state::*;
use crate::qage_1dvec::*;
use crate::qk21_simd::Qk21Simd;
use crate::qk61_simd::Qk61Simd;


#[derive(Clone)]
pub struct Qag_1dvec_parall3 {
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





impl Qag_1dvec_parall3 {
    pub fn qintegrate(&self, fun : &FnVecPa , a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> Qag1DVecParIntegratorResult {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return Qag1DVecParIntegratorResult::new_error(ResultState::Invalid)
        }
        //            first approximation to the integral

        let mut neval = 0;
        let mut last= 1 ;
        let mut result = 0.0;
        let mut abserr = 0.0;
        let mut defabs = 0.0;
        let mut resabs = 0.0;
        let mut result_list = vec![];
        let f = fun.components[0].clone();


        //let qk15 = Qk15 {};
        let qk21 = Qk21Simd {};
        //let qk31 = Qk31 {};
        //let qk41 = Qk41 {};
        //let qk51 = Qk51 {};
        let qk61 = Qk61Simd {};

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }
        match keyf {
            //1 => (result, abserr, defabs, resabs) = qk15.integrate(&*f, a, b),
            2 => (result, abserr, defabs, resabs) = qk21.integrate(&*f, a, b),
            //3 => (result, abserr, defabs, resabs) = qk31.integrate(&*f, a, b),
            //4 => (result, abserr, defabs, resabs) = qk41.integrate(&*f, a, b),
            //5 => (result, abserr, defabs, resabs) = qk51.integrate(&*f, a, b),
            6 => (result, abserr, defabs, resabs) = qk61.integrate(&*f, a, b),
            _ => (),
        }

        result_list.push(ArcResult::new(a,b,result,abserr));


        //           test on accuracy.

        let mut errbnd = epsabs.max(epsrel * result.abs());
        if abserr <= 50.0 * EPMACH * defabs && abserr > errbnd {
            return Qag1DVecParIntegratorResult::new_error(ResultState::BadTolerance)
        }
        if (abserr <= errbnd && abserr != resabs) || abserr == 0.0 {
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
            if keyf == 1 { neval = 30 * neval + 15; }
            return Qag1DVecParIntegratorResult::new(result, abserr, neval, result_list, last)
        }
        if self.limit == 1 {
            return Qag1DVecParIntegratorResult::new_error(ResultState::MaxIteration)
        }


        let mut area = result;
        let mut errsum = abserr;
        let mut errbnd = errbnd;
        let mut iroff1 = Arc::new(Mutex::new(0));
        let mut iroff2 = Arc::new(Mutex::new(0));
        let mut neval = Arc::new(Mutex::new(0));

        let mut bad_tolerance = Arc::new(Mutex::new(0));
        let mut max_iteration = Arc::new(Mutex::new(0));
        let mut bad_function = Arc::new(Mutex::new(0));
        let mut break_flag = Arc::new(Mutex::new(0));


        //          main do-loop
        //           bisect the subinterval with the largest error estimate.



        let pool = threadpool::Builder::new().build();

        while last <= self.limit {

            let mut new_interval = Arc::new(Mutex::new(0));
            let mut new_contributes = Arc::new(Mutex::new(vec![]));
            let limit = self.limit;



                let qk21 = Arc::new(Qk21Simd {});
                let qk61 = Arc::new(Qk61Simd {});

                for component in &result_list {
                    let mut comp = component.inner.lock().unwrap();

                    if comp.error > epsabs.max(epsrel * area.abs()) / last as f64 {
                        let mut new_interval = Arc::clone(&new_interval);
                        let mut area = area.clone();
                        let mut errsum = errsum.clone();
                        let mut iroff1 = Arc::clone(&iroff1);
                        let mut iroff2 = Arc::clone(&iroff2);
                        let mut neval = Arc::clone(&neval);
                        let mut bad_tolerance = Arc::clone(&bad_tolerance);
                        let mut max_iteration = Arc::clone(&max_iteration);
                        let mut bad_function = Arc::clone(&bad_function);
                        let mut break_flag = Arc::clone(&break_flag);
                        let mut component = component.clone();
                        let f = fun.components[0].clone();
                        //let func = fun;


                        let qk21 = Arc::clone(&qk21);
                        let qk61 = Arc::clone(&qk61);


                        let mut new_contributes = Arc::clone(&new_contributes);
                        let limit = limit.clone();

                        pool.execute(move || {
                            let mut new_interval = new_interval.lock().unwrap();
                            *new_interval += 1;
                            drop(new_interval);
                            //let f = fun.components[0].clone();
                            //   !!    let f = &fun.components[0];
                            let mut component = component.inner.lock().unwrap();

                            let a1 = component.a;
                            let b1 = 0.5 * (component.a + component.b);
                            let a2 = b1;
                            let b2 = component.b;


                            let area1: f64;
                            let error1: f64;
                            let area2: f64;
                            let error2: f64;
                            let defab1: f64;
                            let defab2: f64;

                            //let qk15 = Qk15 {};
                            //let qk21 = Qk21Simd2{};
                            //let qk31 = Qk31 {};
                            //let qk41 = Qk41 {};
                            //let qk51 = Qk51 {};
                            //let qk61 = Qk611DVec_Simd{};


                            match keyf {
                                //1 => {
                                //    (area1, error1, _, defab1) = qk15.integrate(&*f, a1, b1);
                                //    (area2, error2, _, defab2) = qk15.integrate(&*f, a2, b2);
                                //},
                                2 => {
                                    (area1, error1, _, defab1) = qk21.integrate(&*f, a1, b1);
                                    (area2, error2, _, defab2) = qk21.integrate(&*f, a2, b2);
                                },
                                //3 => {
                                //    (area1, error1, _, defab1) = qk31.integrate(&*f, a1, b1);
                                //    (area2, error2, _, defab2) = qk31.integrate(&*f, a2, b2);
                                //},
                                //4 => {
                                //    (area1, error1, _, defab1) = qk41.integrate(&*f, a1, b1);
                                //    (area2, error2, _, defab2) = qk41.integrate(&*f, a2, b2);
                                //},
                                //5 => {
                                //    (area1, error1, _, defab1) = qk51.integrate(&*f, a1, b1);
                                //    (area2, error2, _, defab2) = qk51.integrate(&*f, a2, b2);
                                //},
                                6 => {
                                    (area1, error1, _, defab1) = qk61.integrate(&*f, a1, b1);
                                    (area2, error2, _, defab2) = qk61.integrate(&*f, a2, b2);
                                },
                                _ => {
                                    (area1, error1, _, defab1) = (0.0, 0.0, 0.0, 0.0);
                                    (area2, error2, _, defab2) = (0.0, 0.0, 0.0, 0.0);
                                },
                            }


                            //           improve previous approximations to integral
                            //           and error and test for accuracy.

                            let area12 = area1 + area2;
                            let erro12 = error1 + error2;

                            *(neval.lock().unwrap()) += 1;
                            drop(neval);


                            errsum += erro12 - component.error;
                            area += area12 - component.result;
                            let errbnd = epsabs.max(epsrel * area.abs());



                            //  let mut iroff11 = iroff1.lock().unwrap();
                            //  let mut iroff1_temp = *iroff11;
                            //  drop(iroff11);
                            //  let mut iroff22 = iroff2.lock().unwrap();
                            //  let mut iroff2_temp = *iroff22;
                            //  drop(iroff22);

                            if defab1 == error1 || defab2 == error2 {} else {


                                if (component.result - area12).abs() <= 0.00001 * area12.abs() && erro12 >= 0.99 * component.error {
                                    let iroff11 = iroff1.clone();
                                    *(iroff11.lock().unwrap()) += 1;
                                    drop(iroff11);
                                    //iroff1_temp += 1;

                                }
                                if last > 10 && erro12 > component.error {
                                    let iroff22 = iroff2.clone();
                                    *(iroff22.lock().unwrap()) += 1;
                                    drop(iroff22);
                                    //iroff2_temp += 1;
                                }
                            }


                            if errsum > errbnd {

                                //           test for roundoff error.

                                if *iroff1.lock().unwrap() >= 6 || *iroff2.lock().unwrap() >= 20 {
                                    let mut bad_tolerance = bad_tolerance.lock().unwrap();
                                    *bad_tolerance += 1;
                                    return ()
                                    //  QagVecParIntegratorResult::new_error(ResultState::BadTolerance)
                                }

                                drop(iroff1);
                                drop(iroff2);

                                //           set error flag in the case that the number of subintervals
                                //           equals limit.

                                if last == limit {
                                    let mut max_iteration = max_iteration.lock().unwrap();
                                    *max_iteration += 1;
                                    return ()
                                    //QagVecIntegratorResult::new_error(ResultState::MaxIteration)
                                }

                                //           set error flag in the case of bad integrand behaviour
                                //          at a point of the integration range.

                                if a1.abs().max(b2.abs()) <= (1.0 + 100.0 * EPMACH) * (a2.abs() + 1000.0 * UFLOW) {
                                    let mut bad_function = bad_function.lock().unwrap();
                                    *bad_function += 1;
                                    return ()
                                    //QagVecIntegratorResult::new_error(ResultState::BadFunction)
                                }

                                //           append the newly-created intervals to the list.
                            }

                            *component = Result::new(a1, b1, area1, error1);
                            let mut new_contributes = new_contributes.lock().unwrap();
                            new_contributes.push(ArcResult::new(a2, b2, area2, error2));

                            if errsum <= errbnd {
                                *break_flag.lock().unwrap() += 1;
                            }
                        });
                    }


                    let break_fl = break_flag.lock().unwrap();
                    if *break_fl == 1 {
                        break;
                    }
                    drop(break_fl);

                }

            pool.join();
            let mut new_contributes = new_contributes.lock().unwrap();
            let mut new_interval = new_interval.lock().unwrap();
            result_list.append(&mut *new_contributes);
            //result_list.sort_by(|a,b| b.inner.lock().unwrap().error.total_cmp(&a.inner.lock().unwrap().error));
            last += *new_interval;

            (area, errsum, errbnd) = (0.0,0.0,0.0);
            for k in &result_list{
                area += k.inner.lock().unwrap().result;
                errsum += k.inner.lock().unwrap().error;
            }

            errbnd = epsabs.max(epsrel * area.abs());



            if errsum <= errbnd {
                break
            }
        }










        //           compute final result.

        let mut result = area;

        let mut neval = neval.lock().unwrap();
        abserr = errsum;


        if keyf != 1 { *neval = (10 * keyf + 1) * (2 * *neval + 1); }
        if keyf == 1 { *neval = 30 * *neval + 15; }


        Qag1DVecParIntegratorResult::new(result, abserr, *neval, result_list, last)
    }
}


//impl QuadIntegralMethod for Qag_vec_parall {
//    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64) -> QuadIntegratorResult{
//        QuadIntegratorResult::new_qag_vec( self.qintegrate(f,a,b,epsabs,epsrel))
//    }
//}

