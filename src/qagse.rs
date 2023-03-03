use crate::qag_integrator_result::QagIntegratorResult;
use crate::qelg::qelg;
use crate::qk::{EPMACH, OFLOW, Qk, UFLOW};
use crate::result_state::ResultState;
use crate::qk21::*;
use crate::qpsrt::qpsrt;
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;

#[derive(Clone)]
pub struct Qags {
    pub limit : usize,
}

impl Qags {
    pub fn qintegrate(&self, f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> QagIntegratorResult{
        let mut neval: i32;
        let mut ierro = 0;
        let mut result : f64;
        let mut abserr : f64;
        let mut last : usize;
        let mut lastt : usize;
        let mut defabs : f64;
        let mut resabs : f64;
        let mut dres : f64;
        let mut errbnd : f64;
        let mut errmax : f64;
        let mut erlarg : f64;
        let mut erlast : f64;
        let mut ertest : f64;
        let mut correc : f64;
        let mut maxerr : usize;
        let mut area : f64;
        let mut small : f64;
        let mut errsum : f64;
        let mut nrmax : usize;
        let mut nres : usize;
        let mut numrl2 : usize;
        let mut ktmin : i32;
        let mut extrap : bool;
        let mut noext : bool;
        let mut iroff1 : i32;
        let mut iroff2 : i32;
        let mut iroff3 : i32;
        let mut ksgn : i32;
        let mut alist = vec![];
        let mut blist = vec![];
        let mut rlist = vec![];
        let mut rlist2 = vec![];
        let mut elist = vec![];
        let mut iord = vec![];
        let mut res3la = vec![0.0;3];
        let qk21 = Qk21{};

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagIntegratorResult::new_error(ResultState::Invalid)
        }

        (result,abserr,defabs,resabs) = qk21.integrate(&f,a,b);
        dres = result.abs();
        errbnd = epsabs.max(epsrel * dres);
        last = 1;
        lastt = 1;

        alist.push(a);
        blist.push(b);
        rlist.push(result);
        elist.push(abserr);
        iord.push(1);


        if abserr <= 100.0 * EPMACH * defabs && abserr > errbnd {
            return QagIntegratorResult::new_error(ResultState::BadTolerance)
        }
        if (abserr <= errbnd && abserr != resabs) || abserr == 0.0 {
            neval = 42 * last as i32 - 21 ;
            return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, last)
        }
        if self.limit == 1 {
            return QagIntegratorResult::new_error(ResultState::MaxIteration)
        }


        //  inizialization

        rlist2.push(result);
        errmax = abserr;
        maxerr = 1;
        area = result;
        errsum = abserr;
        erlarg = 0.0;
        erlast = 0.0;
        small = 0.0;
        correc = 0.0;
        abserr = OFLOW;
        nrmax = 1;
        nres = 0;
        numrl2 = 2;
        ktmin = 0;
        extrap = false;
        noext = false;
        iroff1 = 0;
        iroff2 = 0;
        iroff3 = 0;
        ksgn = -1;
        ertest = 0.0;

        if dres >= (1.0 - 50.0 * EPMACH) * defabs { ksgn = 1; }

        //           main do-loop

        for last in 2..self.limit+1{
            //  bisect the subinterval with the nrmax-th largest error estimate.

            let a1 = alist[maxerr - 1];
            let b1 = 0.5 * (alist[maxerr - 1] + blist[maxerr - 1]);
            let a2 = b1;
            let b2 = blist[maxerr - 1];
            erlast = errmax;
            lastt = last;
            let area1: f64;
            let error1: f64;
            let area2: f64;
            let error2: f64;
            let defab1: f64;
            let defab2: f64;
            let reseps : f64;
            let abseps : f64;

            (area1, error1, _, defab1) = qk21.integrate(f, a1, b1);
            (area2, error2, _, defab2) = qk21.integrate(f, a2, b2);

            //  improve previous approximations to integral and error and test for accuracy.
            let area12 = area1 + area2;
            let erro12 = error1 + error2;
            errsum += erro12 - errmax;
            area += area12 - rlist[maxerr - 1];

            if defab1 == error1 || defab2 == error2 {}
            else {
                if (rlist[maxerr - 1] - area12).abs() >= 0.00001 * area12.abs() || erro12 < 0.99 * errmax {}
                else{
                    if extrap { iroff2 += 1; }
                    else { iroff1 += 1; }
                }
                if last > 10 && erro12 > errmax { iroff3 += 1; }
            }

            errbnd = epsabs.max(epsrel * area.abs());

            //  test for roundoff error.

            if iroff1 + iroff2 >= 10 || iroff3 >= 20 {
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
            if iroff2 >= 5 {
                ierro = 3;
            }

            //           append the newly-created intervals to the list.


            if error2 > error1
            {
                alist[maxerr - 1] = a2;
                alist.push(a1);
                blist.push(b1);
                rlist[maxerr - 1] = area2;
                rlist.push(area1);
                elist[maxerr - 1] = error2;
                elist.push(error1);
            } else {
                rlist[maxerr - 1] = area1;
                rlist.push(area2);
                alist.push(a2);
                blist[maxerr - 1] = b1;
                blist.push(b2);
                elist[maxerr - 1] = error1;
                elist.push(error2);
            }

            //           call subroutine dqpsrt to maintain the descending ordering
            //           in the list of error estimates and select the subinterval
            //           with nrmax-th largest error estimate (to be bisected next).

            qpsrt(self.limit, last, &mut maxerr, &mut errmax, &elist, &mut iord, &mut nrmax);


            if errsum <= errbnd { break;}
            if last == 2{
                small = (b-a).abs() * 0.375;
                erlarg = errsum;
                ertest = errbnd;
                rlist2.push(area);
                continue;
            }

                if noext { continue; }

                erlarg -= erlast;

                if (b1 - a1).abs() > small {
                    erlarg += erro12;
                }
                //  test whether the interval to be bisected next is the smallest interval.
                if !extrap {
                    if (blist[maxerr-1] - alist[maxerr-1]).abs() > small { continue; }
                    extrap = true;
                    nrmax = 2;
                }
                if !(ierro==3 || erlarg <= ertest){
                    //           the smallest interval has the largest error.
                    //           before bisecting decrease the sum of the errors over the
                    //           larger intervals (erlarg) and perform extrapolation.
                    let id = nrmax;
                    let mut jupbnd = last;
                    if last > (2+self.limit/2)  { jupbnd = self.limit+3-last; }
                    //  do 50 k = id,jupbnd
                    for _k in id..jupbnd+1{
                        maxerr = iord[nrmax-1];
                        errmax = elist[maxerr-1];
                        if (blist[maxerr-1]-alist[maxerr-1]).abs() > small { break; }
                        nrmax += 1;
                    }
                }
                //          perform extrapolation.
                numrl2 += 1;
                //  rlist2(numrl2) = area !!!
                rlist2.push(area);
                let mut epstab = rlist2.clone();
                while epstab.len() < 52{
                    epstab.push(0.0);
                }
                (reseps,abseps) = qelg(&mut numrl2,&mut epstab, &mut res3la , &mut nres );
                ktmin += 1;
                if ktmin > 5 && abserr < 0.001 * errsum {
                    return QagIntegratorResult::new_error(ResultState::Diverge)
                }
                if abseps >= abserr {// go to 70
                    ktmin = 0;
                    abserr = abseps;
                    result = reseps;
                    correc = erlarg;
                    ertest = epsabs.max(epsrel * reseps.abs());
                    if abserr < ertest { break; }
                }

                //  prepare bisection of the smallest interval.
                if numrl2 == 1 { noext = true; }
                maxerr = iord[0];
                errmax = elist[maxerr-1];
                nrmax = 1;
                extrap = false;
                small = small * 0.5;
                erlarg = errsum;
        }

        //  set final result and error estimate.
        if abserr != OFLOW {
            if ierro == 0 {
                if ksgn == -1 && result.abs().max(area.abs()) <= defabs * 0.01 {
                    neval = 42 * lastt as i32 - 21;
                    return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);
                }
                if 0.01 > result/area || result/area > 100.0 || errsum > area.abs() {
                    return QagIntegratorResult::new_error(ResultState::Diverge)
                }

            }
            if ierro == 3 { abserr += correc; }
            if result != 0.0 && area != 0.0 {
                if !( abserr/result.abs() > errsum/area.abs()) {
                    if ksgn == -1 && result.abs().max(area.abs()) <= defabs * 0.01 {
                        neval = 42 * lastt as i32 - 21;
                        return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);
                    }
                    if 0.01 > result/area || result/area > 100.0 || errsum > area.abs() {
                        return QagIntegratorResult::new_error(ResultState::Diverge)
                    }
                }
                result = 0.0;
                for k in 1..lastt+1{
                    result += rlist[k-1];
                }
                abserr = errsum;
                neval = 42 * lastt as i32 - 21;
                return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);

            }
            if abserr > errsum {
                result = 0.0;
                for k in 1..lastt+1{
                    result += rlist[k-1];
                }
                abserr = errsum;
                neval = 42 * lastt as i32 - 21;
                return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);
            }
            if area == 0.0 {
                neval = 42 * lastt as i32 - 21;
                return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);
            }
            if ksgn == -1 && result.abs().max(area.abs()) <= defabs * 0.01 {
                neval = 42 * lastt as i32 - 21;
                return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);
            }
            if 0.01 > result/area || result/area > 100.0 || errsum > area.abs() {
                return QagIntegratorResult::new_error(ResultState::Diverge)
            }
        }
        result = 0.0;
        for k in 1..lastt+1{
            result += rlist[k-1];
        }
        abserr = errsum;
        neval = 42 * lastt as i32 - 21;
        QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt)

    }
    }

impl QuadIntegralMethod for Qags{
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64) -> QuadIntegratorResult{
        QuadIntegratorResult::new_qags( self.qintegrate(f,a,b,epsabs,epsrel))
    }
}