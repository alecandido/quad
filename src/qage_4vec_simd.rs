/*
use std::iter::zip;
use std::simd::{f64x4, Simd, SimdFloat};
use crate::funct_vector::{FnVec, FnVec4};
use crate::qk::*;
use crate::qk15::*;
use crate::qk21::*;
use crate::qk31::*;
use crate::qk41::*;
use crate::qk51::*;
use crate::qk61::*;
use crate::qsrt2::*;
use crate::qag_1dvec_integrator_result::Qag1DVecIntegratorResult;
use crate::qag_vec4_integration_result::ResultVec4;
use crate::qag_vec4_integrator_result::QagVec4IntegratorResult;
use crate::qag_vec_integrator_result::QagVecIntegratorResult;
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;
use crate::result_state::*;
use crate::qage_1dvec::*;
use crate::qage_vec::ResultVec;
use crate::qk21_4vec_simd::Qk21Vec4Simd;
use crate::qk61_4vec_simd2::Qk61Vec4Simd2;
use crate::qk61_4vec_simd::Qk61Vec4Simd;
use crate::qk61_vec::Qk61Vec;


#[derive(Clone)]
pub struct QagVec4Simd {
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





impl QagVec4Simd {
    pub fn qintegrate(&self, f : &FnVec4, a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> QagVec4IntegratorResult {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagVec4IntegratorResult::new_error(ResultState::Invalid)
        }
        //            first approximation to the integral

        const  len : usize = 4;
        let mut neval = 0;
        let mut last= 1 ;
        let mut result = f64x4::splat(0.0);
        let mut abserr = f64x4::splat(0.0);
        let mut defabs = f64x4::splat(0.0);
        let mut resabs = f64x4::splat(0.0);
        let mut result_list = vec![];

        let qk21 = Qk21Vec4Simd {};
        let qk61 = Qk61Vec4Simd2 {};

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }
        match keyf {
            2 => (result, abserr, defabs, resabs) = qk21.integrate(f, a, b),
            6 => (result, abserr, defabs, resabs) = qk61.integrate(f, a, b),
            _ => (),
        }

        result_list.push(ResultVec4::new(a,b,result,abserr));


        //           test on accuracy.

        let mut errbnd = [0.0;len];
        for k in 0..len{
            errbnd[k] = ( epsabs.max( epsrel * result_list[0].result[k].abs()));
            if result_list[0].error[k] <= 50.0 * EPMACH * defabs[k] && result_list[0].error[k] > errbnd[k] {
                return QagVec4IntegratorResult::new_error(ResultState::BadTolerance)
            }
        }
        let mut errbnd = Simd::from_array(errbnd);

        if condition1(&result_list[0].error,&errbnd,&resabs) {
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
            if keyf == 1 { neval = 30 * neval + 15; }
            return QagVec4IntegratorResult::new(result_list[0].result.clone(), result_list[0].error.clone(),
                                               neval, result_list, last)
        }

        if self.limit == 1 {
            return QagVec4IntegratorResult::new_error(ResultState::MaxIteration)
        }

        //          initialization
        let mut area = result_list[0].result.clone();
        let mut errsum = result_list[0].error.clone();
        //let mut iroff1 = [0;len];
        //let mut iroff2 = [0;len];

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.



        let mut area1 = f64x4::splat(0.0);
        let mut error1 = f64x4::splat(0.0);
        let mut area2 = f64x4::splat(0.0);
        let mut error2 = f64x4::splat(0.0);
        let mut defab1 = f64x4::splat(0.0);
        let mut defab2 = f64x4::splat(0.0);

        let mut area12 = f64x4::splat(0.0);
        let mut erro12 = f64x4::splat(0.0);

        let epsabs_simd = f64x4::splat(epsabs);
        let epsrel_simd = f64x4::splat(epsabs);


        while last <= self.limit {

            let last_simd = f64x4::splat(last as f64);

            for k in 0..last {
                if result_list[k].error > epsabs_simd.simd_max(epsrel_simd * area.abs()) / last_simd{
                //condition2(&result_list[k].error,&area,epsabs,epsrel,last as f64){




                    let a1 = result_list[k].a;
                    let b1 = 0.5 * (result_list[k].a + result_list[k].b);
                    let a2 = b1;
                    let b2 = result_list[k].b;



                    match keyf {
                        2 => {
                            (area1, error1, _, defab1) = qk21.integrate(f, a1, b1);
                            (area2, error2, _, defab2) = qk21.integrate(f, a2, b2);
                        },
                        6 => {
                            (area1, error1, _, defab1) = qk61.integrate(f, a1, b1);
                            (area2, error2, _, defab2) = qk61.integrate(f, a2, b2);
                        },
                        _ => (),
                    }



                    //           improve previous approximations to integral
                    //           and error and test for accuracy.



                    area12 = area1 + area2;
                    erro12 = error1 + error2;
                    errsum += erro12 - result_list[k].error;
                    area += area12 - result_list[k].result;

/*
                    if defab1 == error1 || defab2 == error2 {} else {
                        for j in 0..len {
                            if (result_list[k].result[j] - area12[j]).abs() <= 0.00001 * area12[j].abs() && erro12[j] >= 0.99 * result_list[k].error[j] {
                                iroff1[j] += 1;
                            }
                            if last > 10 && erro12[j] > result_list[k].error[j] { iroff2[j] += 1; }
                        }
                    }

 */
                   // for k in 0..len {
                   //     errbnd[k] = epsabs.max(epsrel * area[k].abs());
                   // }

                    errbnd = epsabs_simd.simd_max(epsrel_simd * area.abs());

                    if errsum > errbnd {

                        //           set error flag in the case that the number of subintervals
                        //           equals limit.

                        if last == self.limit {
                            return QagVec4IntegratorResult::new_error(ResultState::MaxIteration)
                        }

                        //           set error flag in the case of bad integrand behaviour
                        //          at a point of the integration range.

                        if a1.abs().max(b2.abs()) <= (1.0 + 100.0 * EPMACH) * (a2.abs() + 1000.0 * UFLOW) {
                            return QagVec4IntegratorResult::new_error(ResultState::BadFunction)
                        }

                        /*
                        for k in 0..len {
                            //           test for roundoff error.

                            if iroff1[k] >= 6 || iroff2[k] >= 20 {
                                return QagVec4IntegratorResult::new_error(ResultState::BadTolerance)
                            }

                        }

                         */
                        //           append the newly-created intervals to the list.
                    }




                    result_list[k] = ResultVec4::new(a1, b1, area1, error1);
                    result_list.push(ResultVec4::new(a2, b2, area2, error2));
                }
                if ! (errsum > errbnd) {
                    break
                }
            }

            /*
            result_list.sort_by(|a,b|
                                    {
                                        let norm1 = (b.error * b.error).reduce_sum();
                                        let norm2 = (a.error * a.error).reduce_sum();
                                        norm1.total_cmp(&norm2)
                                    });

             */


            last = result_list.len();

            if ! (errsum > errbnd) {
                break
            }
        }
        //           compute final result.

        let mut result = f64x4::splat(0.0);
        for k in &result_list {
            result += k.result;
        }


        if keyf != 1 { neval = (10 * keyf + 1) * (2 * last as i32 + 1); }
        if keyf == 1 { neval = 30 * last as i32 + 15; }


        QagVec4IntegratorResult::new(result, errsum, neval, result_list, last)
    }
}





pub fn condition1(a : &Simd<f64,4>, b : &Simd<f64,4>, c : &Simd<f64,4>) -> bool{
    for k in 0..a.lanes(){
        if !((a[k] <= b[k] && a[k] != c[k]) || a[k] == 0.0) {
            return false;
        }
    }
    true
}

pub fn condition2(a : &Simd<f64,4>, b : &Simd<f64,4>, c : f64, d : f64, e : f64) -> bool{
    for k in 0..a.lanes(){
        if a[k] > c.max(d * b[k].abs()) / e { return true; }
    }
    false
}





pub fn vec_norm(a: &Simd<f64,4>) -> f64{
    let mut sum = 0.0;
    for k in 0..4 {
        sum += a[k].powi(2);
    }
    sum.sqrt()
}






//  impl QuadIntegralMethod for QagVec {
//      fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64) -> QuadIntegratorResult{
//          QuadIntegratorResult::new_qag_vec( self.qintegrate(f,a,b,epsabs,epsrel))
//      }
//  }










#[cfg(test)]
mod tests {
    use std::simd::{f64x4, Simd};
    use std::time::Instant;
    use crate::funct_vector::{FnVec, FnVec4};
    use crate::qage_1dvec2::Qag_1dvec2;
    use crate::qage_4vec_simd::QagVec4Simd;
    use crate::qage_vec::QagVec;

    #[test]
    fn test(){
        let key = 6;
        let limit = 1000000;
        let qag_vec = QagVec{key : 6,limit};
        let qag_vec2 = Qag_1dvec2 {key,limit};
        let qag_vec3 = QagVec4Simd{key,limit};

        let f1 = |x:f64|  x.cos();
        let f2 = |x:f64| x.cos();
        let f3 = |x:f64| x.cos() ;
        let f4 = |x:f64| x.cos() ;
        let fun = FnVec{components : vec![Box::new(f1),Box::new(f2),Box::new(f3),Box::new(f4)]};
        let fun2 = FnVec4{components : [Box::new(f1),Box::new(f2),Box::new(f3),Box::new(f4)]};

        let a = 0.0;
        let b = 100000.0;
        let epsabs = 1.0e-5;
        let epsrel = 0.0;
        let max = 10;


        let mut res1;
        let mut res2;
        let mut res3;
        let mut res4;
        let mut res_simd;

        for k in 0..max{
            //let start = Instant::now();
            //let res = qag_vec.qintegrate(&fun,a,b,epsabs,epsrel);
            //println!("vec time : {:?} ",start.elapsed());

            let qag_vec2 = Qag_1dvec2 {key,limit};
            let start = Instant::now();
            res1 = qag_vec2.qintegrate(&f1,a,b,epsabs,epsrel);
            res2 = qag_vec2.qintegrate(&f2,a,b,epsabs,epsrel);
            res3 = qag_vec2.qintegrate(&f3,a,b,epsabs,epsrel);
            res4 = qag_vec2.qintegrate(&f4,a,b,epsabs,epsrel);

            println!(" 1D time : {:?} ",start.elapsed());

            let start = Instant::now();
            res_simd = qag_vec3.qintegrate(&fun2,a,b,epsabs,epsrel);
            println!("simd time : {:?} ",start.elapsed());

            if k == max-1 {
                println!("1d last : {:?}", res1.integration_result.last);
                println!("simd last : {:?}", res_simd.integration_result.last);
            }

        }




        let a = Simd::from_array([0.0,99.0,100.0,101.0]);
        let b = Simd::from_array([0.0,100.0,100.0,100.0]);
        let boo = a.partial_cmp(&b);
        println!("{:?}",boo);

    }
}

 */