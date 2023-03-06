use crate::qag_integrator_result::QagIntegratorResult;
use crate::qelg::qelg;
use crate::qk21::Qk21;
use crate::qk::{EPMACH, OFLOW, Qk, UFLOW};
use crate::qpsrt::qpsrt;
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;
use crate::result_state::ResultState;

#[derive(Clone)]
pub struct Qagp {
    pub limit : usize,
    pub npts2 : usize,
    pub points : Vec<f64>,
}

impl Qagp {
    pub fn qintegrate(&self, f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64)
                      -> QagIntegratorResult {
        let mut neval: i32;
        let mut ierro : i32;
        let mut iroff1 : i32;
        let mut iroff2 : i32;
        let mut iroff3 : i32;
        let mut ksgn : i32;
        let mut npts : usize = self.npts2 - 2;
        let mut last: usize;
        let mut lastt: usize;
        let mut maxerr : usize;
        let mut nrmax : usize;
        let mut nres : usize;
        let mut numrl2 : usize;
        let mut ktmin : usize;
        let mut levcur : usize;
        let mut levmax : usize;
        let mut result : f64 = 0.0;
        let mut abserr : f64 = 0.0;
        let mut resabs : f64;
        let mut sign : f64;
        let mut errmax : f64;
        let mut area : f64;
        let mut erlarg : f64;
        let mut correc : f64;
        let mut ertest : f64;
        let mut extrap : bool;
        let mut noext : bool;
        let mut alist = vec![];
        let mut blist = vec![];
        let mut rlist = vec![];
        let mut elist = vec![];
        let mut iord = vec![];
        let mut res3la = vec![0.0;3];
        let mut level = vec![];
        let mut pts = vec![] ;
        let mut ndin = vec![];
        let mut rlist2 = vec![];
        let qk21 = Qk21{};

        if self.npts2 < 2 || self.limit <= npts || ( epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH)) {
            return QagIntegratorResult::new_error(ResultState::Invalid)
        }

        //            If any break points are provided, sort them into an ascending sequence.

        alist.push(a);
        blist.push(b);
        sign = 1.0;
        lastt = 0;
        correc = 0.0;
        if a > b { sign = -1.0;}
        pts.push(a.min(b));

        if npts != 0 {
            for point in self.points.clone(){
                pts.push(point);
            }
        }
        pts.push(a.max(b));
        let nint = npts + 1;
        let mut a1 = pts[0];

        if npts != 0  {
            let nintp1 = nint + 1;
            for i in 1..nint+1{
                let ip1 = i + 1;
                for j in ip1..nintp1+1{
                    if pts[i-1] <= pts[j-1] { continue; }
                    let temp = pts[i-1];
                    pts[i-1] = pts[j-1];
                    pts[j-1] = temp;
                }
                if pts[0] != a.min(b) || pts[nintp1-1] != a.max(b) {
                    return QagIntegratorResult::new_error(ResultState::Invalid)
                }
            }
        }

        //  Compute first integral and error approximations.

        resabs = 0.0;
        for i in 1..nint+1{
            let b1 = pts[i];
            let area1 : f64;
            let error1 : f64;
            let defabs : f64;
            let resa : f64;

            (area1, error1, defabs, resa) = qk21.integrate(&f,a1,b1);
            abserr += error1;
            result += area1;
            if error1 == resa && error1 != 0.0 { ndin.push(1); }
            else { ndin.push(0); }
            resabs += defabs;
            level.push(0);
            elist.push(error1);
            alist.push(a1);
            blist.push(b1);
            rlist.push(area1);
            iord.push(i);
            a1 = b1;
        }

        let mut errsum = 0.0;
        for i in 1..nint+1{
            if ndin[i-1] == 1 { elist[i-1] = abserr; }
            errsum += elist[i-1];
        }

        // Test on accuracy

        last = nint;
        neval = 21 * nint as i32;
        let mut dres = result.abs();
        let mut errbnd = epsabs.max(epsrel * dres);
        if abserr <= 0.0001 * EPMACH * resabs && abserr > errbnd {
            return QagIntegratorResult::new_error(ResultState::BadTolerance)
        }
        if nint != 1 {
            for i in 1..npts+1{
                let jlow = i + 1;
                let mut k : usize = 0; // safe to initialize it at 0?
                let mut ind1 = iord[i-1];
                for j in jlow..nint+1{
                    let ind2 = iord[j-1];
                    if elist[ind1-1] > elist[ind2-1] { continue; }
                    ind1 = ind2;
                    k = j;
                }
                if ind1 == iord[i-1] { continue; }
                iord[k-1] = iord[i-1];
                iord[i-1] = ind1;
            }

        if self.limit < self.npts2 {
            return QagIntegratorResult::new_error(ResultState::MaxIteration)
        }
                }
        if abserr <= errbnd {
            result = result * sign;
            return QagIntegratorResult::new(result,abserr,neval,alist,blist,rlist,elist,iord,last)
        }

        //  Initialization.

        rlist2.push(result);
        maxerr = iord[0];
        errmax = elist[maxerr-1];
        area = result;
        nrmax = 1;
        nres = 0;
        numrl2 = 1;
        ktmin = 0;
        extrap = false;
        noext = false;
        erlarg = errsum;
        ertest = errbnd;
        levmax = 1;
        iroff1 = 0;
        iroff2 = 0;
        iroff3 = 0;
        ierro = 0;
        abserr = OFLOW;
        ksgn = -1;
        if dres >= (1.0 - 50.0 * EPMACH) * resabs { ksgn = 1; }

        //  main do-loop

        for last in self.npts2..self.limit+1{
            //  bisect the subinterval with the nrmax-th largest error estimate.
            lastt = last;
            levcur = level[maxerr-1] + 1;
            let a1 = alist[maxerr-1];
            let b1 = 0.5 * ( alist[maxerr-1] + blist[maxerr-1] );
            let a2 = b1;
            let b2 = blist[maxerr-1];
            let erlast = errmax;
            let area1: f64;
            let error1: f64;
            let area2: f64;
            let error2: f64;
            let defab1: f64;
            let defab2: f64;
            let reseps : f64;
            let abseps : f64;

            (area1,error1,_,defab1) = qk21.integrate(f,a1,b1);
            (area2,error2,_,defab2) = qk21.integrate(f,a2,b2);

            //  Improve previous approximations to integral and error and test for accuracy.

            neval += 42;
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
            level[maxerr-1] = levcur;
            level.push(levcur);

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
            //           at a point of the integration range.

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
            if noext { continue; }
            erlarg -= erlast;

            if levcur+1 <= levmax  { erlarg+= erro12; }

            //  test whether the interval to be bisected next is the smallest interval.
            if !extrap {
                if level[maxerr-1]+1 <= levmax { continue; }
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
                    if level[maxerr-1]+1 <= levmax { break; }
                    nrmax += 1;
                }
            }
            //          perform extrapolation.
            numrl2 += 1;
            // if numrl2 == last { rlist2.push(area); } !!!!!!!!!!!!
            //  rlist2(numrl2) = area !!!!!!!!!!1
            rlist2.push(area);
            let mut epstab = rlist2.clone();
            while epstab.len() < 52{
                epstab.push(0.0);
            }
            if numrl2 > 2 {
                (reseps, abseps) = qelg(&mut numrl2, &mut epstab, &mut res3la, &mut nres);
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
            }
            maxerr = iord[0];
            errmax = elist[maxerr-1];
            nrmax = 1;
            extrap = false;
            levmax += 1;
            erlarg = errsum;
        }

        //  set final result and error estimate.
        if abserr != OFLOW {
            if ierro == 0 {
                if ksgn == -1 && result.abs().max(area.abs()) <= resabs * 0.01 {
                    result = result * sign;
                    return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);
                }
                if 0.01 > result/area || result/area > 100.0 || errsum > area.abs() {
                    return QagIntegratorResult::new_error(ResultState::Diverge)
                }

            }
            if ierro == 3 { abserr += correc; }
            if result != 0.0 && area != 0.0 {
                if !( abserr/result.abs() > errsum/area.abs()) {
                    if ksgn == -1 && result.abs().max(area.abs()) <= resabs * 0.01 {
                        result = result * sign;
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
                result = result * sign;
                return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);

            }
            if abserr > errsum {
                result = 0.0;
                for k in 1..lastt+1{
                    result += rlist[k-1];
                }
                abserr = errsum;
                result = result * sign;
                return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);
            }
            if area == 0.0 {
                result = result * sign;
                return QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt);
            }
            if ksgn == -1 && result.abs().max(area.abs()) <= resabs * 0.01 {
                result = result * sign;
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
        result = result * sign;
        QagIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, iord, lastt)

    }
}


impl QuadIntegralMethod for Qagp{
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs : f64, epsrel : f64) -> QuadIntegratorResult{
        QuadIntegratorResult::new_qags( self.qintegrate(f,a,b,epsabs,epsrel))
    }
}

