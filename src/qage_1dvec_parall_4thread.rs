use std::sync::{Arc, Mutex};
use na::max;
use rayon::{current_thread_index, join, ThreadPool};
use crate::funct_vector::{FnPa, FnVecP, FnVecPa};
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
pub struct Qag_1dvec_parall_4thread {
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





impl Qag_1dvec_parall_4thread {
    pub fn qintegrate(&self, fun : &FnPa , a : f64, b : f64, epsabs : f64, epsrel : f64)
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
        let f = fun.components.clone();


        let qk21 = Qk21Simd {};
        let qk61 = Qk61Simd {};

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }
        match keyf {
            2 => (result, abserr, defabs, resabs) = qk21.integrate(&*f, a, b),
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

        while last <= self.limit {

            let len = result_list.len();
            let (res1,res3) = result_list.as_mut_slice().split_at_mut(len/2);
            let (res1, res2) = res1.split_at_mut(len/4);
            let (res3, res4) = res3.split_at_mut(len/4);


            let (((mut vec1,area1,error1),(mut vec2, area2, error2)),
                ((mut vec3,area3,error3),(mut vec4, area4, error4))) =
                join(||
                         {join( ||{
                let mut new_contributes = vec![];
                let mut areat1 = 0.0;
                let mut errort1 = 0.0;
                for component in res1 {
                    let mut component = component.inner.lock().unwrap();

                    if component.error > epsabs.max(epsrel * area.abs()) / last as f64 {

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

                        let qk21 = Qk21Simd {};
                        let qk61 = Qk61Simd {};


                        match keyf {
                            2 => {
                                (area1, error1, _, defab1) = qk21.integrate(&*f, a1, b1);
                                (area2, error2, _, defab2) = qk21.integrate(&*f, a2, b2);
                            },
                            6 => {
                                (area1, error1, _, defab1) = qk61.integrate(&*f, a1, b1);
                                (area2, error2, _, defab2) = qk61.integrate(&*f, a2, b2);
                            },
                            _ => {
                                (area1, error1, _, defab1) = (0.0, 0.0, 0.0, 0.0);
                                (area2, error2, _, defab2) = (0.0, 0.0, 0.0, 0.0);
                            },
                        }

                        areat1 += area1 + area2;
                        errort1 += error1 + error2;

                        *component = Result::new(a1, b1, area1, error1);
                        new_contributes.push(ArcResult::new(a2, b2, area2, error2));
                    }
                    else{
                        areat1 += component.result;
                        errort1 += component.error;
                    }
                }
                (new_contributes,areat1,errort1)
            }
                                                                     ,||{
                    let mut new_contributes = vec![];
                    let mut areat2 = 0.0;
                    let mut errort2 = 0.0;
                    for component in res2 {
                        let mut component = component.inner.lock().unwrap();


                        if component.error > epsabs.max(epsrel * area.abs()) / last as f64 {

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

                            let qk21 = Qk21Simd {};
                            let qk61 = Qk61Simd {};


                            match keyf {
                                2 => {
                                    (area1, error1, _, defab1) = qk21.integrate(&*f, a1, b1);
                                    (area2, error2, _, defab2) = qk21.integrate(&*f, a2, b2);
                                },
                                6 => {
                                    (area1, error1, _, defab1) = qk61.integrate(&*f, a1, b1);
                                    (area2, error2, _, defab2) = qk61.integrate(&*f, a2, b2);
                                },
                                _ => {
                                    (area1, error1, _, defab1) = (0.0, 0.0, 0.0, 0.0);
                                    (area2, error2, _, defab2) = (0.0, 0.0, 0.0, 0.0);
                                },
                            }

                            areat2 += area1 + area2;
                            errort2 += error1 + error2;

                            *component = Result::new(a1, b1, area1, error1);
                            new_contributes.push(ArcResult::new(a2, b2, area2, error2));
                        }
                        else{
                            areat2 += component.result;
                            errort2 += component.error;
                        }
                    }
                    (new_contributes,areat2,errort2)

                             }
                         )
                         },
                     ||
                         {join( ||{
                             let mut new_contributes = vec![];
                             let mut areat3 = 0.0;
                             let mut errort3 = 0.0;
                             for component in res3 {
                                 let mut component = component.inner.lock().unwrap();

                                 if component.error > epsabs.max(epsrel * area.abs()) / last as f64 {

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

                                     let qk21 = Qk21Simd {};
                                     let qk61 = Qk61Simd {};


                                     match keyf {
                                         2 => {
                                             (area1, error1, _, defab1) = qk21.integrate(&*f, a1, b1);
                                             (area2, error2, _, defab2) = qk21.integrate(&*f, a2, b2);
                                         },
                                         6 => {
                                             (area1, error1, _, defab1) = qk61.integrate(&*f, a1, b1);
                                             (area2, error2, _, defab2) = qk61.integrate(&*f, a2, b2);
                                         },
                                         _ => {
                                             (area1, error1, _, defab1) = (0.0, 0.0, 0.0, 0.0);
                                             (area2, error2, _, defab2) = (0.0, 0.0, 0.0, 0.0);
                                         },
                                     }

                                     areat3 += area1 + area2;
                                     errort3 += error1 + error2;

                                     *component = Result::new(a1, b1, area1, error1);
                                     new_contributes.push(ArcResult::new(a2, b2, area2, error2));
                                 }
                                 else{
                                     areat3 += component.result;
                                     errort3 += component.error;
                                 }
                             }
                             (new_contributes, areat3, errort3)
                         }
                                ,||{
                                 let mut new_contributes = vec![];
                                 let mut areat4 = 0.0;
                                 let mut errort4 = 0.0;
                                 for component in res4 {
                                     let mut component = component.inner.lock().unwrap();


                                     if component.error > epsabs.max(epsrel * area.abs()) / last as f64 {

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

                                         let qk21 = Qk21Simd {};
                                         let qk61 = Qk61Simd {};


                                         match keyf {
                                             2 => {
                                                 (area1, error1, _, defab1) = qk21.integrate(&*f, a1, b1);
                                                 (area2, error2, _, defab2) = qk21.integrate(&*f, a2, b2);
                                             },
                                             6 => {
                                                 (area1, error1, _, defab1) = qk61.integrate(&*f, a1, b1);
                                                 (area2, error2, _, defab2) = qk61.integrate(&*f, a2, b2);
                                             },
                                             _ => {
                                                 (area1, error1, _, defab1) = (0.0, 0.0, 0.0, 0.0);
                                                 (area2, error2, _, defab2) = (0.0, 0.0, 0.0, 0.0);
                                             },
                                         }

                                         areat4 += area1 + area2;
                                         errort4 += error1 + error2;

                                         *component = Result::new(a1, b1, area1, error1);
                                         new_contributes.push(ArcResult::new(a2, b2, area2, error2));
                                     }
                                     else{
                                         areat4 += component.result;
                                         errort4 += component.error;
                                     }
                                 }
                                 (new_contributes, areat4, errort4)

                             }
                         )
                         }





                );

            last += vec1.len() + vec2.len() + vec3.len() + vec4.len();
            result_list.append(&mut vec1);
            result_list.append(&mut vec2);
            result_list.append(&mut vec3);
            result_list.append(&mut vec4);

            //  (area, errsum, errbnd) = (0.0,0.0,0.0);
            //  for k in &result_list{
            //      area += k.inner.lock().unwrap().result;
            //      errsum += k.inner.lock().unwrap().error;
            //  }
            area = area1 + area2 + area3 + area4;
            errsum = error1 + error2 + error3 + error4;

            errbnd = epsabs.max(epsrel * area.abs());



            if errsum <= errbnd {
                break
            }
        }










        //           compute final result.

        let mut result = area;

        let mut neval = last as i32;
        abserr = errsum;


        if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
        if keyf == 1 { neval = 30 * neval + 15; }


        Qag1DVecParIntegratorResult::new(result, abserr, neval, result_list, last)
    }
}


//impl QuadIntegralMethod for Qag_vec_parall {
//    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64) -> QuadIntegratorResult{
//        QuadIntegratorResult::new_qag_vec( self.qintegrate(f,a,b,epsabs,epsrel))
//    }
//}

