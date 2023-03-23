///     Purpose:
///         The function determines the limit of a given sequence of
///         approximations, by means of the epsilon algorithm of
///         p.wynn. An estimate of the absolute error is also given.
///         The condensed epsilon table is computed. Only those
///         elements needed for the computation of the next diagonal
///         are preserved.
///
///     Parameters:
///         n   :   i32
///                 epstab(n) contains the new element in the
///                 first column of the epsilon table.
///
///         epstab  :   Vec<f64>
///                     Vector of dimension 52 containing the elements
///                     of the two lower diagonals of the triangular
///                     epsilon table. The elements are numbered
///                     starting at the right-hand corner of the
///                     triangle.
///
///         result  :   f64
///                     resulting approximation to the integral
///
///         abserr  :   f64
///                     estimate of the absolute error computed from
///                     result and the 3 previous results
///
///         res3la  :   Vec<f64>
///                     vector of dimension 3 containing the last 3
///                     results
///
///         nres    :   i32
///                     number of calls to the routine
///                     (should be zero at first call)
///
///     List of major variables:
///
///         e0,e1,e2,e3 :   the 4 elements on which the computation of a new
///                         element in the epsilon table is based
///                                      e0
///                                e3    e1    new
///                                      e2
///         newelm  :   number of elements to be computed in the new
///                     diagonal
///         error   :   error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
///         result  :   the element in the new diagonal with least value
///                     of error
///

use crate::qk::{EPMACH, OFLOW};

pub fn qelg(n : &mut usize, epstab : &mut Vec<f64>, res3la : &mut Vec<f64>, nres : &mut usize) -> (f64, f64){
    println!("{n}");
    *nres += 1;
    let mut abserr = OFLOW;
    let mut result = epstab[*n-1];
    if *n < 3 {
        abserr = abserr.max(5.0 * EPMACH * result.abs());
        return (result,abserr);
    }
    let limexp = 50;
    epstab[*n+1] = epstab[*n-1];
    let newelm = (*n-1)/2;
    epstab[*n-1] = OFLOW;
    let num = *n;
    let mut k1 = *n;
    for i in 1..newelm+1 {
        let k2 = k1 - 1;
        let k3 = k1 - 2;
        let mut res = epstab[k1 + 1];
        let e0 = epstab[k3 - 1];
        let e1 = epstab[k2 - 1];
        let e2 = res;
        let e1abs = e1.abs();
        let delta2 = e2 - e1;
        let err2 = delta2.abs();
        let tol2 = e1abs.max(e2.abs()) * EPMACH;
        let delta3 = e1 - e0;
        let err3 = delta3.abs();
        let tol3 = e1abs.max(e0.abs()) * EPMACH;

        if !(err2 > tol2 || err3 > tol3) {
            //           if e0, e1 and e2 are equal to within machine
            //           accuracy, convergence is assumed.
            //           result = e2
            //           abserr = abs(e1-e0)+abs(e2-e1)
            result = res;
            abserr = err2 + err3;
            abserr = abserr.max(5.0 * EPMACH * result.abs());
            if *n == limexp-1 {
                *n = 2 * (limexp / 2) - 2;
            }
            print!("convergence is assumed");
            return (result, abserr);
        }
        let e3 = epstab[k1 - 1];
        epstab[k1 - 1] = e1;
        let delta1 = e1 - e3;
        let err1 = delta1.abs();
        let tol1 = e1abs.max(e3.abs()) * EPMACH;

        //  if two elements are very close to each other, omit a part of the table by adjusting
        //  the value of n

        if err1 <= tol1 || err2 <= tol2 || err3 <= tol3 {
            *n = 2 * i - 1;
            break;
        }
        let ss = 1.0 / delta1 + 1.0 / delta2 - 1.0 / delta3;
        let epsinf = (ss * e1).abs();

        //           compute a new element and eventually adjust the value of result.

        if epsinf <= 0.0001 {
            *n = 2 * i - 1;
            break;
        }
        res = e1 + 1.0 / ss;
        epstab[k1 - 1] = res;
        k1 = k1 - 2;
        let error = err2 + (res - e2).abs() + err3;
        if error > abserr { continue; }
        abserr = error;
        result = res;
    }


            //           shift the table.

            if *n == limexp-1 {
                *n = 2 * (limexp / 2) - 2;
            }
            print!(",{n}");
            let mut ib = 1;
            if (num / 2) * 2 == num {
                ib = 2;
            }
            let ie = newelm + 1;
            for _i in 1..ie+1{
                let ib2 = ib + 2;
                epstab[ib-1] = epstab[ib2-1];
                ib = ib2;
            }
            if num != *n {
                let mut indx = num - *n + 1;
                //  do 70 i = 1,n
                for i in 1..*n+1{
                    epstab[i-1]= epstab[indx];
                    indx += 1;
                }

            }

            if *nres >= 4 {
                //  compute error estimate.
                abserr = (result - res3la[2]).abs() + (result - res3la[1]).abs() +
                    (result - res3la[0]).abs();
                res3la[0] = res3la[1];
                res3la[1] = res3la[2];
                res3la[2] = result;
                abserr = abserr.max(5.0 * EPMACH * result.abs());
                return (result,abserr);
            }

            res3la[*nres-1] = result;
            abserr = OFLOW;
            abserr = abserr.max(5.0 * EPMACH * result.abs());

            (result,abserr)
        }

