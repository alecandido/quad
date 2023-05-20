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
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;
use crate::result_state::*;



#[derive(Clone)]
pub struct QagVecNosortFindmaxIroff {
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





impl QagVecNosortFindmaxIroff {
    pub fn qintegrate(&self, f : &dyn Fn(f64)->[f64;4], a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> QagVecIntegratorResult {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagVecIntegratorResult::new_error(ResultState::Invalid)
        }
        //            first approximation to the integral

        let mut neval = 0;
        let mut lastt= 1 ;
        let mut result = [0.0;4];
        let mut abserr = [0.0;4];
        let mut defabs = [0.0;4];
        let mut resabs = [0.0;4];
        let mut alist = vec![];
        let mut blist = vec![];
        let mut rlist = vec![];
        let mut elist = vec![];
        alist.push(a);
        blist.push(b);

        //let qk15 = Qk15 {};
        //let qk21 = Qk21 {};
        //let qk31 = Qk31 {};
        //let qk41 = Qk41 {};
        //let qk51 = Qk51 {};
        let qk61 = Qk614Vec {};

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
        rlist.push(result);
        elist.push(abserr);

        //           test on accuracy.

        let mut errbnd = [epsabs.max(epsrel * result[0].abs()),
            epsabs.max(epsrel * result[1].abs()),epsabs.max(epsrel * result[2].abs()),
            epsabs.max(epsrel * result[3].abs())];
        for k in 0..4{
            if abserr[k] <= 50.0 * EPMACH * defabs[k] && abserr[k] > errbnd[k] {
                return QagVecIntegratorResult::new_error(ResultState::BadTolerance)
            }
            if (abserr[k] <= errbnd[k] && abserr[k] != resabs[k]) || abserr[k] == 0.0 {
                if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
                if keyf == 1 { neval = 30 * neval + 15; }
                return QagVecIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, lastt)
            }
        }

        if self.limit == 1 {
            return QagVecIntegratorResult::new_error(ResultState::MaxIteration)
        }

        //          initialization
        let mut errmax = abserr;
        let mut maxerr = 0;
        let mut maxerr2 = 0;
        let mut area = result;
        let mut errsum = abserr;
        let mut iroff1 = 0;
        let mut iroff2 = 0;

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.

        for last  in 2..self.limit + 1 {
            let a1 = alist[maxerr];
            let b1 = 0.5 * (alist[maxerr] + blist[maxerr]);
            let a2 = b1;
            let b2 = blist[maxerr];

            let area1: [f64;4];
            let error1: [f64;4];
            let area2: [f64;4];
            let error2: [f64;4];
            let defab1: [f64;4];
            let defab2: [f64;4];

            lastt = last;

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
                errsum[k] += erro12[k] - errmax[k];
                area[k] += area12[k] - rlist[maxerr][k];
            }

            if defab1[maxerr2] == error1[maxerr2] || defab2[maxerr2] == error2[maxerr2] {} else {
                if (rlist[maxerr][maxerr2] - area12[maxerr2]).abs() <= 0.00001 * area12[maxerr2].abs() && erro12[maxerr2] >= 0.99 * errmax[maxerr2] {
                    iroff1 += 1;
                }
                if last > 10 && erro12[maxerr2] > errmax[maxerr2] { iroff2 += 1; }
            }
            for k in 0..4{
                errbnd[k] = epsabs.max(epsrel * area[k].abs());
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

                //           append the newly-created intervals to the list.

            }




            alist[maxerr] = a2;
            alist.push(a1);
            blist.push(b1);
            rlist[maxerr] = area2;
            rlist.push(area1);
            elist[maxerr] = error2;
            elist.push(error1);
            //  println!("{last} error2 > error1, rlis {:?}, elist{:?}, alist {:?}, blist {:?}",
            //          rlist, elist, alist, blist);




            if errsum[0] <= errbnd[0] && errsum[1] <= errbnd[1] && errsum[2] <= errbnd[2] &&
                errsum[3] <= errbnd[3]{
                break
            }

            let maxerr_list = [elist.iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a[0].total_cmp(&b[0]))
                .map(|(index, _)| index).unwrap(),
                elist.iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a[1].total_cmp(&b[1]))
                    .map(|(index, _)| index).unwrap(),
                elist.iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a[2].total_cmp(&b[2]))
                    .map(|(index, _)| index).unwrap(),
                elist.iter()
                    .enumerate()
                    .max_by(|(_, a), (_, b)| a[3].total_cmp(&b[3]))
                    .map(|(index, _)| index).unwrap()];
            let errmax_list = [elist[maxerr_list[0]][0], elist[maxerr_list[1]][1],
                elist[maxerr_list[2]][2], elist[maxerr_list[3]][3]];

            maxerr2 = errmax_list.iter()
                .enumerate()
                .max_by(|(_, a), (_, b)| a.total_cmp(b))
                .map(|(index, _)| index).unwrap();

            maxerr = maxerr_list[maxerr2];
            errmax = elist[maxerr];
        }
        //           compute final result.

        result = [0.0;4];
        for k in 1..lastt+1 {
            for j in 0..4{
                result[j] += rlist[k-1][j];
            }
        }
        abserr = errsum;


        if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
        if keyf == 1 { neval = 30 * neval + 15; }


        QagVecIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist,  lastt)
    }
}




#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qage::Qag;
    use crate::qage_vec_nosort_findmax::QagVecNosortFindmax;
    use crate::qage_vec_nosort_findmax_iroff::QagVecNosortFindmaxIroff;
    use crate::qk61::Qk61;
    use crate::qk61_4vec::Qk614Vec;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f1 = |x:f64| x.cos();
        let f2 = |x:f64| x.sin();
        let f3 = |x:f64| f1(x) + f2(x);
        let f4 = |x:f64| - 2.0 * f3(x);

        let a = 0.0;
        let b = 1000000.0;
        let key = 6;
        let epsabs = 1e-4;
        let epsrel = 0.0;
        let limit = 1000000;
        let qag = Qag {key,limit};
        let qag_vec = QagVecNosortFindmax{key,limit};
        let qag_vec_iroff = QagVecNosortFindmaxIroff{key,limit};
        let f = |x:f64| [f1(x),f2(x),f3(x),f4(x)];


        let mut res1;
        let mut res2;
        let mut res3;
        let mut res4;
        let mut res_vec;
        let mut res_vec_iroff ;

        for k in 0..5 {
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
            println!("vec {:?}", start.elapsed());
            if k == 0{
                println!("normal {:?},{:?},{:?},{:?}",res1.integration_result.result,
                         res2.integration_result.result,res3.integration_result.result,
                         res4.integration_result.result);
                println!("vec {:?}",res_vec.integration_result.result);
                println!("vec iroff {:?}",res_vec_iroff.integration_result.result);
                println!("normal {:?},{:?},{:?},{:?}",res1.integration_result.last,
                         res2.integration_result.last,res3.integration_result.last,
                         res4.integration_result.last);
                println!("vec {:?}",res_vec.integration_result.last);
                println!("vec {:?}",res_vec_iroff.integration_result.last);
            }
        }

    }
}


