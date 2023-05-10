use std::iter::zip;
use crate::funct_vector::FnVec;
use crate::qk::*;
use crate::qk15::*;
use crate::qk21::*;
use crate::qk31::*;
use crate::qk41::*;
use crate::qk51::*;
use crate::qk61::*;
use crate::qsrt2::*;
use crate::qag_1dvec_integrator_result::Qag1DVecIntegratorResult;
use crate::qag_vec_integrator_result_pre::QagVecIntegratorResultPre;
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;
use crate::result_state::*;
use crate::qage_1dvec::*;
use crate::qk61_vec_pre::Qk61VecPre;

#[derive(Clone,Debug)]
pub struct ResultVecPre {
    pub a : f64,
    pub b : f64,
    pub result : Vec<f64>,
    pub error : Vec<f64>,
}

impl ResultVecPre {
    pub fn new(a : f64, b : f64, result : Vec<f64>, error : Vec<f64>) -> Self{
        Self{
            a,b,result,error
        }
    }
}

#[derive(Clone)]
pub struct QagVecPre {
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





impl QagVecPre {
    pub fn qintegrate(&self, f : &FnVec, a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> QagVecIntegratorResultPre {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagVecIntegratorResultPre::new_error(ResultState::Invalid)
        }
        //            first approximation to the integral

        let len = f.components.len();
        let mut neval = 0;
        let mut last= 1 ;
        let mut result = vec![];
        let mut abserr = vec![];
        let mut defabs = vec![];
        let mut resabs = vec![];
        let mut result_list = vec![];


        //let qk15 = Qk15 {};
        //let qk21 = Qk21 {};
        //let qk31 = Qk31 {};
        //let qk41 = Qk41 {};
        //let qk51 = Qk51 {};
        let qk61 = Qk61VecPre {};

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }
        match keyf {
            //1 => (result, abserr, defabs, resabs) = qk15.integrate(f, a, b),
            //2 => (result, abserr, defabs, resabs) = qk21.integrate(f, a, b),
            //3 => (result, abserr, defabs, resabs) = qk31.integrate(f, a, b),
            //4 => (result, abserr, defabs, resabs) = qk41.integrate(f, a, b),
            //5 => (result, abserr, defabs, resabs) = qk51.integrate(f, a, b),
            6 => (result, abserr, defabs, resabs) = qk61.integrate(f, a, b),
            _ => (),
        }

        result_list.push(ResultVecPre::new(a, b, result, abserr));


        //           test on accuracy.

        let mut errbnd = vec![];
        for k in 0..len{
            errbnd.push( epsabs.max( epsrel * result_list[0].result[k].abs()));
            if result_list[0].error[k] <= 50.0 * EPMACH * defabs[k] && result_list[0].error[k] > errbnd[k] {
                return QagVecIntegratorResultPre::new_error(ResultState::BadTolerance)
            }
        }
        if condition1(&result_list[0].error,&errbnd,&resabs) {
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
            if keyf == 1 { neval = 30 * neval + 15; }
            return QagVecIntegratorResultPre::new(result_list[0].result.clone(), result_list[0].error.clone(),
                                                  neval, result_list, last)
        }

        if self.limit == 1 {
            return QagVecIntegratorResultPre::new_error(ResultState::MaxIteration)
        }

        //          initialization
        let mut area = result_list[0].result.clone();
        let mut errsum = result_list[0].error.clone();
        let mut iroff1 = vec![0;len];
        let mut iroff2 = vec![0;len];

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.


        while last <= self.limit {

            let mut new_interval = 0;

            for k in 0..last {
                if condition2(&result_list[k].error,&area,epsabs,epsrel,last as f64){
                    new_interval += 1;
                    let a1 = result_list[k].a;
                    let b1 = 0.5 * (result_list[k].a + result_list[k].b);
                    let a2 = b1;
                    let b2 = result_list[k].b;


                    let mut area1 : Vec<f64> = vec![];
                    let mut error1 : Vec<f64> = vec![];
                    let mut area2 : Vec<f64> = vec![];
                    let mut error2 : Vec<f64> = vec![];
                    let mut defab1 : Vec<f64> = vec![];
                    let mut defab2 : Vec<f64> = vec![];


                    match keyf {
                        //1 => {
                        //    (area1, error1, _, defab1) = qk15.integrate(f, a1, b1);
                        //    (area2, error2, _, defab2) = qk15.integrate(f, a2, b2);
                        //},
                        //2 => {
                        //    (area1, error1, _, defab1) = qk21.integrate(f, a1, b1);
                        //    (area2, error2, _, defab2) = qk21.integrate(f, a2, b2);
                        //},
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
                            (area1, error1, _, defab1) = (vec![0.0], vec![0.0], vec![0.0], vec![0.0]);
                            (area2, error2, _, defab2) = (vec![0.0], vec![0.0], vec![0.0], vec![0.0]);
                        },
                    }


                    //           improve previous approximations to integral
                    //           and error and test for accuracy.

                    neval += 1;


                    let area12 = vec_sum(&area1,&area2);
                    let erro12 = vec_sum(&error1, &error2);
                    errsum_update(&mut errsum,&erro12,&result_list[k].error);
                    area_update(&mut area,&area12,&result_list[k].result);


                    if defab1 == error1 || defab2 == error2 {} else {
                        for j in 0..len {
                            if (result_list[k].result[j] - area12[j]).abs() <= 0.00001 * area12[j].abs() && erro12[j] >= 0.99 * result_list[k].error[j] {
                                iroff1[j] += 1;
                            }
                            if last > 10 && erro12[j] > result_list[k].error[j] { iroff2[j] += 1; }
                        }
                    }
                    for k in 0..len {
                        errbnd[k] = epsabs.max(epsrel * area[k].abs());
                    }

                    if gt(&errsum ,&errbnd) {

                        //           set error flag in the case that the number of subintervals
                        //           equals limit.

                        if last == self.limit {
                            return QagVecIntegratorResultPre::new_error(ResultState::MaxIteration)
                        }

                        //           set error flag in the case of bad integrand behaviour
                        //          at a point of the integration range.

                        if a1.abs().max(b2.abs()) <= (1.0 + 100.0 * EPMACH) * (a2.abs() + 1000.0 * UFLOW) {
                            return QagVecIntegratorResultPre::new_error(ResultState::BadFunction)
                        }

                        for k in 0..len {
                            //           test for roundoff error.

                            if iroff1[k] >= 6 || iroff2[k] >= 20 {
                                return QagVecIntegratorResultPre::new_error(ResultState::BadTolerance)
                            }

                        }
                        //           append the newly-created intervals to the list.
                    }


                    result_list[k] = ResultVecPre::new(a1, b1, area1, error1);
                    result_list.push(ResultVecPre::new(a2, b2, area2, error2));
                }
            }

            result_list.sort_by(|a,b| vec_norm(&b.error).total_cmp(&vec_norm(&a.error)));
            last += new_interval;

            if ! gt(&errsum, &errbnd) {
                break
            }
        }
        //           compute final result.

        let mut result = vec![];
        for j in 0..len {
            let mut comp_result = 0.0;
            for k in &result_list {
                comp_result += k.result[j];
            }
            result.push(comp_result);
        }


        if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
        if keyf == 1 { neval = 30 * neval + 15; }


        QagVecIntegratorResultPre::new(result, errsum, neval, result_list, last)
    }
}


pub fn le( a : &Vec<f64>, b : &Vec<f64> ) -> bool{
    for (i,j) in zip(a,b){
        if i <= j { return true;}
    }
    false
}

pub fn gt( a : &Vec<f64>, b : &Vec<f64> ) -> bool{
    for (i,j) in zip(a,b){
        if i > j { return true;}
    }
    false
}

pub fn condition1(a : &Vec<f64>, b : &Vec<f64>, c : &Vec<f64>) -> bool{
    for k in 0..a.len(){
        if !((a[k] <= b[k] && a[k] != c[k]) || a[k] == 0.0) {
            return false;
        }
    }
    true
}

pub fn condition2(a : &Vec<f64>, b : &Vec<f64>, c : f64, d : f64, e : f64) -> bool{
    for k in 0..a.len(){
        if a[k] > c.max(d * b[k].abs()) / e { return true; }
    }
    false
}

pub fn vec_sum(a : &Vec<f64>, b : &Vec<f64>) -> Vec<f64>{
    let mut sum = vec![];
    for (i,j) in zip(a,b){
        sum.push(i+j);
    }
    sum
}

pub fn vec_diff(a : &Vec<f64>, b : &Vec<f64>) -> Vec<f64>{
    let mut diff = vec![];
    for (i,j) in zip(a,b){
        diff.push(i-j);
    }
    diff
}

pub fn vec_norm(a: &Vec<f64>) -> f64{
    let mut sum = 0.0;
    for comp in a {
        sum += comp.powi(2);
    }
    sum.sqrt()
}


pub fn errsum_update(a : &mut Vec<f64>, b : & Vec<f64>, c : & Vec<f64>){
    for k in 0..a.len(){
        a[k] += b[k] - c[k];
    }
}

pub fn area_update(a : &mut Vec<f64>, b : & Vec<f64>, c : & Vec<f64>){
    for k in 0..a.len(){
        a[k] += b[k] - c[k];
    }
}



//  impl QuadIntegralMethod for QagVec {
//      fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64) -> QuadIntegratorResult{
//          QuadIntegratorResult::new_qag_vec( self.qintegrate(f,a,b,epsabs,epsrel))
//      }
//  }










#[cfg(test)]
mod tests {
    use std::simd::Simd;
    use std::time::Instant;
    use crate::funct_vector::FnVec;
    use crate::qage_1dvec2::Qag_1dvec2;
    use crate::qage_vec_pre::QagVecPre;

    #[test]
    fn test(){
        let key = 6;
        let limit = 100000;
        let qag_vec = QagVecPre {key,limit};

        let f1 = |x:f64|  x.cos();
        let f2 = |x:f64| x.cos();
        let f3 = |x:f64| x.cos() ;
        let fun = FnVec{components : vec![Box::new(f1),Box::new(f2),Box::new(f3)]};

        let a = 0.0;
        let b = 100.0;
        let epsabs = 1.0;
        let epsrel = 0.0;

        let start = Instant::now();
        let res = qag_vec.qintegrate(&fun,a,b,epsabs,epsrel);
        println!("time : {:?} ... {:?}",start.elapsed(),res);



        let qag_vec2 = Qag_1dvec2 {key,limit};
        let start = Instant::now();
        let res1 = qag_vec2.qintegrate(&f1,a,b,epsabs,epsrel);
        let res2 = qag_vec2.qintegrate(&f2,a,b,epsabs,epsrel);
        let res3 = qag_vec2.qintegrate(&f3,a,b,epsabs,epsrel);
        println!("time : {:?} ... {:?}",start.elapsed(),res1);
        println!("{:?}",res2);
        println!("{:?}",res3);

        for k in 0..100{
            let start = Instant::now();
            let res = qag_vec.qintegrate(&fun,a,b,epsabs,epsrel);
            println!("vec time : {:?} ... {:?}",start.elapsed(),res);



            let qag_vec2 = Qag_1dvec2 {key,limit};
            let start = Instant::now();
            let res1 = qag_vec2.qintegrate(&f1,a,b,epsabs,epsrel);
            let res2 = qag_vec2.qintegrate(&f2,a,b,epsabs,epsrel);
            let res3 = qag_vec2.qintegrate(&f3,a,b,epsabs,epsrel);
            println!(" 1D time : {:?} ... {:?}",start.elapsed(),res1);
        }
        let a0 = [1,2,3,4];
        let a = Simd::from(a0);
        println!("{:?}",a);

    }
}