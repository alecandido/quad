use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::hash;
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
use crate::qk61_vec_norm::Qk61VecNorm;
use crate::quad_integral_method::QuadIntegralMethod;
use crate::quad_integrator_result::QuadIntegratorResult;
use crate::result_state::*;




#[derive(Clone)]
pub struct QagVecNorm {
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





impl QagVecNorm {
    pub fn qintegrate(&self, f : &dyn Fn(f64)->Vec<f64>, a : f64, b : f64, epsabs : f64, epsrel : f64)
                      ->  (Vec<f64>,f64) {

        //if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
        //    return QagVecIntegratorResult::new_error(ResultState::Invalid)
        //}


        //            first approximation to the integral

        let mut neval = 0;
        let mut last= 1 ;
        let mut result = vec![];
        let mut abserr = 0.0;
        let mut rounderr  = 0.0;

        let qk61 = Qk61VecNorm {};

        let mut keyf = self.key;
        if self.key <= 0 { keyf = 1; }
        if self.key >= 7 { keyf = 6; }
        match keyf {
            6 => (result, abserr, rounderr) = qk61.integrate(f, a, b),
            _ => (),
        }

        //           test on accuracy.

        let mut errbnd = epsabs.max(epsrel * norm_vec(&result));

        let mut interval_cache = HashMap::from([((Myf64{x:a},Myf64{x:b}),result.clone())]);
        let mut heap = BinaryHeap::new();
        heap.push(HeapItem::new((a,b),abserr));

        //  DA VEDERE !!!!!!!!!!!!!
        //if abserr <= 50.0 * EPMACH * defabs[k] && abserr[k] > errbnd[k] {
        //    return QagVecIntegratorResult::new_error(ResultState::BadTolerance)
        //}


        if abserr <= errbnd{
            if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
            if keyf == 1 { neval = 30 * neval + 15; }
            println!("first iter is enough ");
            return (result,abserr)

            //return QagVecIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist, last)
        }

        if self.limit == 1 {
            //return QagVecIntegratorResult::new_error(ResultState::MaxIteration)
        }

        //          initialization
        let mut iroff1 = 0;
        let mut iroff2 = 0;

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.




        while last < self.limit{
            last += 1;
            let old_interval = heap.pop().unwrap();
            let ((x,y),old_err) = (old_interval.interval,old_interval.err);
            let old_res = interval_cache.remove(&(Myf64{x}, Myf64{x:y})).unwrap();

            let mut result1 = vec![];
            let mut abserr1 = 0.0;
            let mut rounderr1  = 0.0;

            let mut result2 = vec![];
            let mut abserr2 = 0.0;
            let mut rounderr2  = 0.0;

            let a1 = x;
            let b1 = 0.5 * (x + y);
            let a2 = b1;
            let b2 = y;


            match keyf {
                6 => {
                    (result1, abserr1, rounderr1) = qk61.integrate(f, a1, b1);
                    (result2, abserr2, rounderr2) = qk61.integrate(f, a2, b2);
                },
                _ => (),
            }

            //println!("{:?}",result);
            res_update(&mut result, & result1, & result2, &old_res);

            interval_cache.insert((Myf64{x:a1},Myf64{x:b1}),result1);
            interval_cache.insert((Myf64{x:a2},Myf64{x:b2}),result2);

            heap.push(HeapItem::new((a1,b1),abserr1));
            heap.push(HeapItem::new((a2,b2),abserr2));


            abserr += -old_err + abserr1 + abserr2;
            rounderr += rounderr1 + rounderr2;

            errbnd = epsabs.max(epsrel * norm_vec(&result));

            //println!("{:?}",abserr);

            if abserr <= errbnd{ break;}
        }









        //println!("{:?}",interval_cache);
        //println!("{:?}",heap);
        //println!("{:?}",result);

        return (result,abserr)

/*

        while last  < self.limit + 1 {
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
                    errsum[k] += erro12[k] - elist[i][k];
                    area[k] += area12[k] - rlist[i][k];
                }


                if defab1[comp] == error1[comp] || defab2[comp] == error2[comp] {} else {
                    if (rlist[i][comp] - area12[comp]).abs() <= 0.00001 * area12[comp].abs() && erro12[comp] >= 0.99 * errmax[comp] {
                        iroff1 += 1;
                    }
                    if last > 10 && erro12[comp] > errmax[comp] { iroff2 += 1; }
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





            }
           //           compute final result.

        result = [0.0;4];
        for k in 1..last+1 {
            for j in 0..4{
                result[j] += rlist[k-1][j];
            }
        }
        abserr = errsum;


        if keyf != 1 { neval = (10 * keyf + 1) * (2 * neval + 1); }
        if keyf == 1 { neval = 30 * neval + 15; }


        QagVecIntegratorResult::new(result, abserr, neval, alist, blist, rlist, elist,  last)


 */

    }
}


pub fn norm_vec(v : &[f64]) -> f64{
    let mut norm = 0.0;
    for comp in v{
        norm += comp.powi(2);
    }
    norm = norm.sqrt();
    norm
}

pub fn res_update(v : &mut[f64], w: &[f64], z : &[f64], y : &[f64]){
    for k  in 0..v.len(){
        v[k] += w[k] + z[k] - y[k];
    }
}

//#[derive(PartialEq, PartialOrd, Eq, Ord)]
#[derive(Debug)]
pub struct HeapItem {
    interval : (f64,f64),
    err : f64,
}

impl HeapItem {
    pub fn new( interval : (f64,f64) , err : f64) -> Self{
        Self{ interval,err}
    }
}

impl Eq for HeapItem{}

impl PartialEq for HeapItem{
    fn eq(&self, other: &Self) -> bool {
        self.err == other.err
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.err).partial_cmp(&other.err).unwrap()
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}


/*
#[derive(Hash)]
pub struct Interval{
    start : f64,
    end : f64,
}

impl Interval{
    pub fn new( start : f64 , end : f64) -> Self{
        Self{ start, end}
    }
    fn key(&self) -> u64 {
        self.0.to_bits()
    }
}

impl PartialEq for Interval{
    fn eq(&self, other: &Self) -> bool {
        (self.start,self.end) == (other.start,other.end)
    }
}

impl Eq for Interval{}
 */

#[derive(Debug)]
pub struct Myf64{
    x : f64,
}
impl Myf64 {
    fn key(&self) -> u64 {
        self.x.to_bits()
    }
}

impl hash::Hash for Myf64 {
    fn hash<H>(&self, state: &mut H)
        where
            H: hash::Hasher,
    {
        self.key().hash(state)
    }
}

impl PartialEq for Myf64 {
    fn eq(&self, other: &Myf64) -> bool {
        self.key() == other.key()
    }
}

impl Eq for Myf64{}


#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qage::Qag;
    use crate::qage_vec_norm::QagVecNorm;
    use crate::qage_vec_nosort::QagVecNosort;
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

        let f = |x:f64| vec![f1(x),f2(x),f3(x),f4(x)];


        let fv = |x:f64| [f1(x),f2(x),f3(x),f4(x)];

        let a = 0.0;
        let b = 10000000.0;
        let key = 6;
        let epsabs = 1.0e-2;
        let epsrel = 0.0;
        let limit = 1000000;
        let max = 3;
        let qag_vec_norm = QagVecNorm{key,limit};
        let qag_vec_nosort = QagVecNosort{key,limit};


        let mut res;
        let mut res_vec;



        for k in 0..max {
            let start = Instant::now();
            res = qag_vec_norm.qintegrate(&f,a,b,epsabs,epsrel);
            println!("scipy time : {:?}",start.elapsed());
            let start = Instant::now();
            res_vec = qag_vec_nosort.qintegrate(&fv, a, b,epsabs,epsrel);
            println!("my time : {:?}",start.elapsed());

            if k == max-1{
                println!("{:?}",res);
                println!("{:?}",res_vec.integration_result.result);
                println!("{:?}",res_vec.integration_result.abserr);
                println!("{:?}",res_vec.integration_result.last);
            }

        }

    }
}


