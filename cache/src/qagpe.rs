///     Purpose:
///         The routine calculates an approximation result to a given
///         definite integral i = integral of f over (a,b), hopefully
///         satisfying following claim for accuracy abs(i-result) <=
///         max(epsabs,epsrel*abs(i)). Break points of the integration
///         interval, where local difficulties of the integrand may
///         occur(e.g. singularities,discontinuities),provided by user.
///
///     Parameters:
///
///     On entry:
///         f       :   f64
///                     Integrand function f(x).
///
///         a       :   f64
///                     Lower limit of integration.
///
///         b       :   f64
///                     Upper limit of integration.
///
///         npts2   :   i32
///                     Number equal to two more than the number of
///                     user-supplied break points within the integration
///                     range, npts2 >= 2.
///                     If npts2 < 2, the function will end with error_state = Invalid.
///
///         points  :   Vec<f64>
///                     Vector of dimension npts2, the first (npts2-2)
///                     elements of which are the user provided break
///                     points. If these points do not constitute an
///                     ascending sequence there will be an automatic.
///                     sorting.
///
///         epsabs  :   f64
///                     Absolute accuracy requested.
///         epsrel  :   f64
///                     Relative accuracy requested.
///                     If  epsabs <= 0 && epsrel < max(50*rel.mach.acc.,0.5d-28),
///                     the function will end with error_state = Invalid.
///
///         limit   :   i32
///                     Gives an upper bound on the number of subintervals
///                     in the partition of (a,b), limit >= npts2
///                     if limit.lt.npts2, the routine will end with
///                     error_state = Invalid.
///
///     On return:
///         result  :   f64
///                     Approximation to the integral.
///
///         abserr  :   f64
///                     Estimate of the modulus of the absolute error,
///                     which should equal or exceed abs(i-result).
///
///         neval   :   i32
///                     Number of integrand evaluations.
///
///         error_state :
///
///             Success :   normal and reliable termination of the routine. it is assumed that the
///                         requested accuracy has been achieved.
///
///             MaxIteration    :   Maximum number of subdivisions allowed has been achieved.
///                                 One can allow more subdivisions by increasing the value of
///                                 limit (and taking the according dimension adjustments into
///                                 account). However, if this yields no improvement it is advised
///                                 to analyze the integrand in order to determine the integration
///                                 difficulties. If the position of a local difficulty can be
///                                 determined (i.e. singularity,discontinuity within the interval),
///                                 it should be supplied to the function as an element of the
///                                 vector points. If necessary an appropriate special-purpose
///                                 integrator must be used, which is designed for handling the
///                                 type of difficulty involved.
///             BadTolerance    :   The occurrence of roundoff error is detected, which prevents
///                                 the requested tolerance from being achieved.
///                                 The error may be under-estimated.
///             BadFunction     :   Extremely bad integrand behaviour occurs at some points of the
///                                 integration interval.
///                 = 4 the algorithm does not converge.
///                     roundoff error is detected in the
///                     extrapolation table. it is presumed that
///                     the requested tolerance cannot be
///                     achieved, and that the returned result is
///                     the best which can be obtained.
///                 = 5 the integral is probably divergent, or
///                     slowly convergent. it must be noted that
///                     divergence can occur with any other value
///                     of ier.gt.0.
///                 = 6 the input is invalid because
///                     npts2.lt.2 or
///                     break points are specified outside
///                     the integration range or
///                     (epsabs.le.0 and
///                      epsrel.lt.max(50*rel.mach.acc.,0.5d-28))
///                     or limit.lt.npts2.
///                     result, abserr, neval, last, rlist(1),
///                     and elist(1) are set to zero. alist(1) and
///                     blist(1) are set to a and b respectively.
///
///            alist  - double precision
///                     vector of dimension at least limit, the first
///                      last  elements of which are the left end points
///                     of the subintervals in the partition of the given
///                     integration range (a,b)
///
///            blist  - double precision
///                     vector of dimension at least limit, the first
///                      last  elements of which are the right end points
///                     of the subintervals in the partition of the given
///                     integration range (a,b)
///
///            rlist  - double precision
///                     vector of dimension at least limit, the first
///                      last  elements of which are the integral
///                     approximations on the subintervals
///
///            elist  - double precision
///                     vector of dimension at least limit, the first
///                      last  elements of which are the moduli of the
///                     absolute error estimates on the subintervals
///
///            pts    - double precision
///                     vector of dimension at least npts2, containing the
///                     integration limits and the break points of the
///                     interval in ascending sequence.
///
///            level  - integer
///                     vector of dimension at least limit, containing the
///                     subdivision levels of the subinterval, i.e. if
///                     (aa,bb) is a subinterval of (p1,p2) where p1 as
///                     well as p2 is a user-provided break point or
///                     integration limit, then (aa,bb) has level l if
///                     abs(bb-aa) = abs(p2-p1)*2**(-l).
///
///            ndin   - integer
///                     vector of dimension at least npts2, after first
///                     integration over the intervals (pts(i)),pts(i+1),
///                     i = 0,1, ..., npts2-2, the error estimates over
///                     some of the intervals may have been increased
///                     artificially, in order to put their subdivision
///                     forward. if this happens for the subinterval
///                     numbered k, ndin(k) is put to 1, otherwise
///                     ndin(k) = 0.
///
///            iord   - integer
///                     vector of dimension at least limit, the first k
///                     elements of which are pointers to the
///                     error estimates over the subintervals,
///                     such that elist(iord(1)), ..., elist(iord(k))
///                     form a decreasing sequence, with k = last
///                     if last.le.(limit/2+2), and k = limit+1-last
///                     otherwise
///
///            last   - integer
///                     number of subintervals actually produced in the
///                     subdivisions process


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

            let (area1, error1, defabs, resa) = qk21.integrate(&f,a1,b1);
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
            println!("{i} abserr: {abserr}");
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
        println!("abserr:{abserr} errbnd:{errbnd}");
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
            let reseps: f64;
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

            //if a1.abs().max(b2.abs()) <= (1.0 + 100.0 * EPMACH) * (a2.abs() + 1000.0 * UFLOW) {
            //    return QagIntegratorResult::new_error(ResultState::BadFunction)
            //}
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
            if numrl2 - 1 == rlist2.len()  { rlist2.push(area); }
            else { rlist2[numrl2-1] = area ;}
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
                if abseps >= abserr {
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

