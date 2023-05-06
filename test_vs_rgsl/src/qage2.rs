use crate::qk::*;
use crate::qk21_simd::Qk211DVec_Simd;
use crate::qk61_simd::Qk611DVec_Simd;
use crate::qsrt2::*;
use crate::qag_integrator_result::*;
use crate::result_state::*;


#[derive(Clone)]
pub struct Qag2 {
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





impl Qag2 {
    pub fn qintegrate(&self, f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> QagIntegratorResult {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagIntegratorResult::new_error(ResultState::Invalid)
        }
        //            first approximation to the integral

        let mut neval = 0;
        let mut lastt= 1 ;
        let mut result = 0.0;
        let mut abserr = 0.0;
        let mut defabs = 0.0;
        let mut resabs = 0.0;
        let mut alist = vec![];
        let mut blist = vec![];
        let mut rlist = vec![];
        let mut elist = vec![];
        let mut iord = vec![];
        alist.push(a);
        blist.push(b);

        //let qk15 = Qk15 {};
        let qk21 = Qk211DVec_Simd {};
        //let qk31 = Qk31 {};
        //let qk41 = Qk41 {};
        //let qk51 = Qk51 {};
        let qk61 = Qk611DVec_Simd {};

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
        rlist.push(result);
        elist.push(abserr);
        iord.push(1);

        //           test on accuracy.

        let mut errbnd = epsabs.max(epsrel * result.abs());
        if abserr <= 50.0 * EPMACH * defabs && abserr > errbnd {
            return QagIntegratorResult::new_error(ResultState::BadTolerance)
        }
        if (abserr <= errbnd && abserr != resabs) || abserr == 0.0 {
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
            if keyf == 1 { neval = 30 * neval + 15; }
            return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt)
        }
        if self.limit == 1 {
            return QagIntegratorResult::new_error(ResultState::MaxIteration)
        }

        //          initialization
        let mut errmax = abserr;
        let mut maxerr = 1;
        let mut area = result;
        let mut errsum = abserr;
        let mut iroff1 = 0;
        let mut iroff2 = 0;

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.

        for last  in 2..self.limit + 1 {
            let a1 = alist[maxerr - 1];
            let b1 = 0.5 * (alist[maxerr - 1] + blist[maxerr - 1]);
            let a2 = b1;
            let b2 = blist[maxerr - 1];

            let area1: f64;
            let error1: f64;
            let area2: f64;
            let error2: f64;
            let defab1: f64;
            let defab2: f64;

            lastt = last;

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
            errsum += erro12 - errmax;
            area += area12 - rlist[maxerr - 1];


            if defab1 == error1 || defab2 == error2 {} else {
                if (rlist[maxerr - 1] - area12).abs() <= 0.00001 * area12.abs() && erro12 >= 0.99 * errmax {
                    iroff1 += 1;
                }
                if last > 10 && erro12 > errmax { iroff2 += 1; }
            }
            errbnd = epsabs.max(epsrel * area.abs());


            if errsum > errbnd {

                //           test for roundoff error.

                if iroff1 >= 6 || iroff2 >= 20 {
                    return QagIntegratorResult::new_error(ResultState::BadTolerance)
                }

                //           set error flag in the case that the number of subintervals
                //           equals limit.

                if last == self.limit {
                    return QagIntegratorResult::new_error(ResultState::MaxIteration)
                }

                //           set error flag in the case of bad integrand behaviour
                //          at a point of the integration range.

                if a1.abs().max(b2.abs()) <= (1.0 + 100.0 * EPMACH) * (a2.abs() + 1000.0 * UFLOW) {
                    return QagIntegratorResult::new_error(ResultState::BadFunction)
                }

                //           append the newly-created intervals to the list.
            }

            if error2 > error1
            {
                alist[maxerr - 1] = a2;
                alist.push(a1);
                blist.push(b1);
                rlist[maxerr - 1] = area2;
                rlist.push(area1);
                elist[maxerr - 1] = error2;
                elist.push(error1);
                //  println!("{last} error2 > error1, rlis {:?}, elist{:?}, alist {:?}, blist {:?}",
                //          rlist, elist, alist, blist);
            } else {
                rlist[maxerr - 1] = area1;
                rlist.push(area2);
                alist.push(a2);
                blist[maxerr - 1] = b1;
                blist.push(b2);
                elist[maxerr - 1] = error1;
                elist.push(error2);
                //  println!("{last} error1 > error2, rlis {:?}, elist{:?}, alist {:?}, blist {:?}",
                //           rlist, elist, alist, blist);
            }


            //          call subroutine dqpsrt to maintain the descending ordering
            //          in the list of error estimates and select the subinterval
            //           with the largest error estimate (to be bisected next).

            //  println!("Entering qpsrt with : limit={limit},last={last},maxerr={maxerr},errmax={errmax},\
            //  elist={:?},iord={:?},nrmax={nrmax}",elist,iord);
            qpsrt2(self.limit, last, &mut maxerr, &mut errmax, &elist, &mut iord);
            //  println!("Exiting qpsrt with : limit={limit},last={last},maxerr={maxerr},errmax={errmax},\
            //  elist={:?},iord={:?},nrmax={nrmax}",elist,iord);
            if errsum <= errbnd {
                break
            }
        }
        //           compute final result.

        result = 0.0;
        for k in 1..lastt+1 {
            result += rlist[k-1];
        }
        abserr = errsum;


        if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
        if keyf == 1 { neval = 30 * neval + 15; }


        QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt)
    }
}






