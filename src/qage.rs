use crate::qk::*;
use crate::qk15::*;
use crate::qk21::*;
use crate::qk31::*;
use crate::qk41::*;
use crate::qk51::*;
use crate::qk61::*;
use crate::qpsrt::*;
use crate::qag_integration_result::*;
use crate::qag_integrator_result::QagIntegratorResult;


#[derive(Clone)]
pub struct Qna{}

impl Qna{
    pub fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64, key : i32,
                 limit : usize) -> QagIntegratorResult {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagIntegratorResult::new_error(ResultState::Invalid)
        }
        //            first approximation to the integral

        let mut neval = 0;
        let mut lastt ;
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

        let qk15 = Qk15 {};
        let qk21 = Qk21 {};
        let qk31 = Qk31 {};
        let qk41 = Qk41 {};
        let qk51 = Qk51 {};
        let qk61 = Qk61 {};

        let mut keyf = key;
        if key <= 0 { keyf = 1; }
        if key >= 7 { keyf = 6; }
        match key {
            1 => (result, abserr, defabs, resabs) = qk15.integrate(f, a, b),
            2 => (result, abserr, defabs, resabs) = qk21.integrate(f, a, b),
            3 => (result, abserr, defabs, resabs) = qk31.integrate(f, a, b),
            4 => (result, abserr, defabs, resabs) = qk41.integrate(f, a, b),
            5 => (result, abserr, defabs, resabs) = qk51.integrate(f, a, b),
            6 => (result, abserr, defabs, resabs) = qk61.integrate(f, a, b),
            _ => (),
        }
        //  println!("First integral:result:{result},abserr:{abserr}");
        lastt = 1;
        rlist.push(result);
        elist.push(abserr);
        iord.push(1);

        //           test on accuracy.

        let mut errbnd = epsabs.max(epsrel * result.abs());
        if abserr <= 50.0 * EPMACH * defabs && abserr > errbnd {
            return QagIntegratorResult::new_error(ResultState::BadTolerance)
        }
        if limit == 1 {
            return QagIntegratorResult::new_error(ResultState::MaxIteration)
        }
        if (abserr <= errbnd && abserr != resabs) || abserr == 0.0 {
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
            if keyf == 1 { neval = 30 * neval + 15; }
            return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt)
        }
        //          initialization
        let mut errmax = abserr;
        let mut maxerr = 1;
        let mut area = result;
        let mut errsum = abserr;
        let mut nrmax = 1;
        let mut iroff1 = 0;
        let mut iroff2 = 0;

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.

        for last  in 2..limit + 1 {
            //  println!("{last},alist : {:?}, blist:{:?}",alist,blist);
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
                1 => {
                    (area1, error1, _, defab1) = qk15.integrate(f, a1, b1);
                    (area2, error2, _, defab2) = qk15.integrate(f, a2, b2);
                },
                2 => {
                    (area1, error1, _, defab1) = qk21.integrate(f, a1, b1);
                    (area2, error2, _, defab2) = qk21.integrate(f, a2, b2);
                },
                3 => {
                    (area1, error1, _, defab1) = qk31.integrate(f, a1, b1);
                    (area2, error2, _, defab2) = qk31.integrate(f, a2, b2);
                },
                4 => {
                    (area1, error1, _, defab1) = qk41.integrate(f, a1, b1);
                    (area2, error2, _, defab2) = qk41.integrate(f, a2, b2);
                },
                5 => {
                    (area1, error1, _, defab1) = qk51.integrate(f, a1, b1);
                    (area2, error2, _, defab2) = qk51.integrate(f, a2, b2);
                },
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


            if errsum <= errbnd {} else {

                //           test for roundoff error and eventually set error flag.

                if iroff1 >= 6 || iroff2 >= 20 {
                    return QagIntegratorResult::new_error(ResultState::BadTolerance)
                }

                //           set error flag in the case that the number of subintervals
                //           equals limit.

                if last == limit {
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
            qpsrt(limit, last, &mut maxerr, &mut errmax, &elist, &mut iord, &mut nrmax);
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


