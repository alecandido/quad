use crate::qk::*;
use crate::qk15::*;
use crate::qk21::*;
use crate::qk31::*;
use crate::qk41::*;
use crate::qk51::*;
use crate::qk61::*;
use crate::qsrt2::*;
use crate::qag_integrator_result::*;
use crate::qag_vec_integrator_result::QagVecIntegratorResult;
use crate::qk61_4vec::Qk614Vec;
use crate::qk61_vec_earlystop::Qk61VecES;
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;
use crate::result_state::*;



#[derive(Clone)]
pub struct QagVecNosortEs {
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





impl QagVecNosortEs {
    pub fn qintegrate(&self, f : &dyn Fn(f64)->[f64;4], a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> QagVecIntegratorResult {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagVecIntegratorResult::new_error(ResultState::Invalid)
        }
        //            first approximation to the integral

        let mut neval = 0;
        let mut last= 1 ;
        let mut result = [0.0;4];
        let mut abserr = [0.0;4];
        let mut defabs = [0.0;4];
        let mut resabs = [0.0;4];
        let mut alist = vec![];
        let mut blist = vec![];
        let mut rlist = vec![];
        let mut elist = vec![];
        let mut flag = [true;4];
        alist.push(a);
        blist.push(b);

        //let qk15 = Qk15 {};
        //let qk21 = Qk21 {};
        //let qk31 = Qk31 {};
        //let qk41 = Qk41 {};
        //let qk51 = Qk51 {};
        let qk61 = Qk61VecES {};

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }
        match keyf {
            //1 => (result, abserr, defabs, resabs) = qk15.integrate(f, a, b),
            //2 => (result, abserr, defabs, resabs) = qk21.integrate(f, a, b),
            //3 => (result, abserr, defabs, resabs) = qk31.integrate(f, a, b),
            //4 => (result, abserr, defabs, resabs) = qk41.integrate(f, a, b),
            //5 => (result, abserr, defabs, resabs) = qk51.integrate(f, a, b),
            6 => (result, abserr, defabs, resabs) = qk61.integrate(f, a, b,&flag),
            _ => (),
        }
        rlist.push(result);
        elist.push(abserr);

        //           test on accuracy.

        let mut errbnd = [epsabs.max(epsrel * result[0].abs()),
            epsabs.max(epsrel * result[1].abs()),epsabs.max(epsrel * result[2].abs()),
            epsabs.max(epsrel * result[3].abs())];
        for k in 0..4 {
            if abserr[k] <= 50.0 * EPMACH * defabs[k] && abserr[k] > errbnd[k] {
                return QagVecIntegratorResult::new_error(ResultState::BadTolerance)
            }
        }

        if ((abserr[0] <= errbnd[0] && abserr[0] != resabs[0]) || abserr[0] == 0.0)  &&
            ((abserr[1] <= errbnd[1] && abserr[1] != resabs[1]) || abserr[1] == 0.0) &&
            ((abserr[2] <= errbnd[2] && abserr[2] != resabs[2]) || abserr[2] == 0.0) &&
            ((abserr[3] <= errbnd[3] && abserr[3] != resabs[3]) || abserr[3] == 0.0) {
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
            if keyf == 1 { neval = 30 * neval + 15; }
            return QagVecIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, last)
        }



        if self.limit == 1 {
            return QagVecIntegratorResult::new_error(ResultState::MaxIteration)
        }

        //          initialization
        let mut errmax = abserr;
        let mut area = result;
        let mut errsum = abserr;
        let mut iroff1 = 0;
        let mut iroff2 = 0;


        for k in 0..4{
            if flag[k] &&  errsum[k] <= errbnd[k] { flag[k] = false;}
        }


        //          main do-loop
        //           bisect the subinterval with the largest error estimate.

        while last  < self.limit + 1 {
            for i in 0..last{
                let mut comp = 5;
                if flag[0] && elist[i][0] > errbnd[0] / last as f64 && elist[i][0] > 0.1 * errmax[0] { comp = 0;}
                else if flag[1] && elist[i][1] > errbnd[1] / last as f64 && elist[i][1] > 0.1 * errmax[1] { comp = 1;}
                else if flag[2] && elist[i][2] > errbnd[2] / last as f64 && elist[i][2] > 0.1 * errmax[2] { comp = 2;}
                else if flag[3] && elist[i][3] > errbnd[3] / last as f64 && elist[i][3] > 0.1 * errmax[3] { comp = 3;}
                if comp != 5{
                    last += 1;
                    let a1 = alist[i];
                    let b1 = 0.5 * (alist[i] + blist[i]);
                    let a2 = b1;
                    let b2 = blist[i];

                    let area1: [f64;4];
                    let error1: [f64;4];
                    let area2: [f64;4];
                    let error2: [f64;4];
                    let defab1: [f64;4];
                    let defab2: [f64;4];



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
                            (area1, error1, _, defab1) = qk61.integrate(f, a1, b1, &flag);
                            (area2, error2, _, defab2) = qk61.integrate(f, a2, b2, &flag);
                        },
                        _ => {
                            (area1, error1, _, defab1) = ([0.0;4], [0.0;4], [0.0;4], [0.0;4]);
                            (area2, error2, _, defab2) = ([0.0;4], [0.0;4], [0.0;4], [0.0;4]);
                        },
                    }


                    //           improve previous approximations to integral
                    //           and error and test for accuracy.

                    neval += 1;
                    let area12 = [area1[0] + area2[0], area1[1] + area2[1], area1[2] + area2[2],
                        area1[3] + area2[3]];
                    let erro12 = [error1[0] + error2[0], error1[1] + error2[1],
                        error1[2] + error2[2], error1[3] + error2[3]];
                    for k in 0..4{
                        if flag[k] == true{
                            errsum[k] += erro12[k] - elist[i][k];
                            area[k] += area12[k] - rlist[i][k];
                        }
                    }


                    if defab1[comp] == error1[comp] || defab2[comp] == error2[comp] {} else {
                        if (rlist[i][comp] - area12[comp]).abs() <= 0.00001 * area12[comp].abs() && erro12[comp] >= 0.99 * errmax[comp] {
                            iroff1 += 1;
                        }
                        if last > 10 && erro12[comp] > errmax[comp] { iroff2 += 1; }
                    }
                    for k in 0..4{
                        if flag[k] == true {errbnd[k] = epsabs.max(epsrel * area[k].abs());}
                    }


                    if errsum[0] > errbnd[0] || errsum[1] > errbnd[1] || errsum[2] > errbnd[2] ||
                        errsum[3] > errbnd[3]{

                        //           test for roundoff error.

                        if iroff1 >= 6 || iroff2 >= 20 {
                            return QagVecIntegratorResult::new_error(ResultState::BadTolerance)
                        }

                        //           set error flag in the case that the number of subintervals
                        //           equals limit.

                        if last == self.limit {
                            return QagVecIntegratorResult::new_error(ResultState::MaxIteration)
                        }

                        //           set error flag in the case of bad integrand behaviour
                        //          at a point of the integration range.

                        if a1.abs().max(b2.abs()) <= (1.0 + 100.0 * EPMACH) * (a2.abs() + 1000.0 * UFLOW) {
                            return QagVecIntegratorResult::new_error(ResultState::BadFunction)
                        }
                    }




                    alist[i] = a2;
                    alist.push(a1);
                    blist.push(b1);
                    rlist[i] = area2;
                    rlist.push(area1);
                    elist[i] = error2;
                    elist.push(error1);
                    //  println!("{last} error2 > error1, rlis {:?}, elist{:?}, alist {:?}, blist {:?}",
                    //          rlist, elist, alist, blist);




                    for k in 0..4{
                        if flag[k] &&  errsum[k] <= errbnd[k] {
                            flag[k] = false;
                            //println!("{k} at {last}");
                        }
                    }

                    if !flag[0] && !flag[1] && !flag[2] && !flag[3]{
                        break
                    }
                }
            }

            if !flag[0] && !flag[1] && !flag[2] && !flag[3]{
                break
            }

            let mut maxerr_list = [0;4];
            for k in 0..4{
                if flag[k] == true{
                    maxerr_list[k] = elist.iter()
                        .enumerate()
                        .max_by(|(_, a), (_, b)| a[k].total_cmp(&b[k]))
                        .map(|(index, _)| index).unwrap();
                }
            }

            errmax = [elist[maxerr_list[0]][0], elist[maxerr_list[1]][1],
                elist[maxerr_list[2]][2], elist[maxerr_list[3]][3]];

        }
        //           compute final result.

        //result = [0.0;4];
        //for k in 1..last+1 {
        //    for j in 0..4{
        //        result[j] += rlist[k-1][j];
        //    }
        //}
        result = area;
        abserr = errsum;


        if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
        if keyf == 1 { neval = 30 * neval + 15; }


        QagVecIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist,  last)
    }
}




#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qage::Qag;
    use crate::qage_vec_nosort::QagVecNosort;
    use crate::qage_vec_nosort_earlystop::QagVecNosortEs;
    use crate::qage_vec_nosort_findmax::QagVecNosortFindmax;
    use crate::qage_vec_nosort_findmax_iroff::QagVecNosortFindmaxIroff;
    use crate::qk61::Qk61;
    use crate::qk61_4vec::Qk614Vec;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f1 = |x:f64| (1.0/(x*x)).cos();
        let f2 = |x:f64| (1.0/(x*x)).sin();
        let f3 = |x:f64| (1.0/(x*x)).sin();
        let f4 = |x:f64| (1.0/(x*x)).sin();

        let a = 0.0;
        let b = 1.0;
        let key = 6;
        let epsabs = 1e-3;
        let epsrel = 0.0;
        let limit = 1000000;
        let max = 3;
        let qag = Qag {key,limit};
        let qag_vec = QagVecNosortFindmax{key,limit};
        let qag_vec_iroff = QagVecNosortFindmaxIroff{key,limit};
        let qag_vec_nosort = QagVecNosort{key,limit};
        let qag_vec_nosort_es = QagVecNosortEs{key,limit};
        let f = |x:f64| [f1(x),f2(x),f3(x),f4(x)];


        let mut res1;
        let mut res2;
        let mut res3;
        let mut res4;
        let mut res_vec;
        let mut res_vec_iroff ;
        let mut res_vec_nosort;
        let mut res_vec_nosort_es;

        for k in 0..max {
            let start = Instant::now();
            res1 = qag.qintegrate(&f1, a, b,epsabs,epsrel);
            res2 = qag.qintegrate(&f2,a,b,epsabs,epsrel);
            res3 = qag.qintegrate(&f3,a,b,epsabs,epsrel);
            res4 = qag.qintegrate(&f4,a,b,epsabs,epsrel);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            res_vec = qag_vec.qintegrate(&f, a, b,epsabs,epsrel);
            println!("vec {:?}", start.elapsed());
            let start = Instant::now();
            res_vec_iroff = qag_vec_iroff.qintegrate(&f, a, b,epsabs,epsrel);
            println!("vec iroff {:?}", start.elapsed());
            let start = Instant::now();
            res_vec_nosort = qag_vec_nosort.qintegrate(&f, a, b,epsabs,epsrel);
            println!("vec nosort {:?}", start.elapsed());
            let start = Instant::now();
            res_vec_nosort_es = qag_vec_nosort_es.qintegrate(&f, a, b,epsabs,epsrel);
            println!("vec nosort es {:?}", start.elapsed());
            if k == max-1{
                println!("normal {:?},{:?},{:?},{:?}",res1.integration_result.result,
                         res2.integration_result.result,res3.integration_result.result,
                         res4.integration_result.result);
                println!("vec {:?}",res_vec.integration_result.result);
                println!("vec iroff {:?}",res_vec_iroff.integration_result.result);
                println!("vec nosort {:?}",res_vec_nosort.integration_result.result);
                println!("vec nosort es {:?}",res_vec_nosort_es.integration_result.result);
                println!("normal {:?},{:?},{:?},{:?}",res1.integration_result.last,
                         res2.integration_result.last,res3.integration_result.last,
                         res4.integration_result.last);
                //println!("vec {:?}",res_vec.integration_result.last);
                //println!("vec iroff{:?}",res_vec_iroff.integration_result.last);
                println!("vec nosort{:?}",res_vec_nosort.integration_result.last);
                println!("vec nosort es{:?}",res_vec_nosort_es.integration_result.last);

                println!("normal {:?},{:?},{:?},{:?}",res1.result_state,
                         res2.result_state,res3.result_state,
                         res4.result_state);
                //println!("vec {:?}",res_vec.integration_result.last);
                //println!("vec iroff{:?}",res_vec_iroff.integration_result.last);
                println!("vec nosort{:?}",res_vec_nosort.result_state);
                println!("vec nosort es{:?}",res_vec_nosort_es.result_state);



            }
        }

    }
}


