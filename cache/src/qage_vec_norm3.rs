use std::cmp::Ordering;
use std::collections::{BinaryHeap, HashMap};
use std::hash;
use crate::qk::*;
use crate::qag_vec_norm_integrator_result::QagVecNormIntegratorResult;
use crate::qk61_vec_norm2::Qk61VecNorm2;
use crate::result_state::*;




#[derive(Clone)]
pub struct QagVecNorm3 {
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





impl QagVecNorm3 {
    pub fn qintegrate<const n :usize>(&self, f : &dyn Fn(f64)->[f64;n], a : f64, b : f64, epsabs : f64, epsrel : f64)
                                      ->  QagVecNormIntegratorResult<n> {

        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return QagVecNormIntegratorResult::new_error(ResultState::Invalid)
        }


        //            first approximation to the integral

        let mut neval = 0;
        let mut last= 1 ;
        let mut result = [0.0;n];
        let mut abserr = 0.0;
        let mut rounderr  = 0.0;

        let qk61 = Qk61VecNorm2 {};

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
            return QagVecNormIntegratorResult::new(result,abserr,neval,last)
        }

        if self.limit == 1 {
            return QagVecNormIntegratorResult::new_error(ResultState::MaxIteration)
        }

        //          initialization

        //          main do-loop
        //           bisect the subinterval with the largest error estimate.




        while last < self.limit{
            let mut to_process = vec![];
            let mut err_sum = 0.0;


            while to_process.len() < 128 && heap.len() != 0{
                let old_interval = heap.pop().unwrap();
                let ((x,y),old_err) = (old_interval.interval,old_interval.err);
                let old_res = interval_cache.remove(&(Myf64{x}, Myf64{x:y})).unwrap();
                err_sum += old_err;
                to_process.push((x,y,old_err,old_res));
                if err_sum > abserr - errbnd / 8.0 { break}
            }


            for comp in to_process{
                last += 1;

                let mut result1 = [0.0;n];
                let mut abserr1 = 0.0;
                let mut rounderr1  = 0.0;

                let mut result2 = [0.0;n];
                let mut abserr2 = 0.0;
                let mut rounderr2  = 0.0;

                let a1 = comp.0;
                let b1 = 0.5 * (comp.0 + comp.1);
                let a2 = b1;
                let b2 = comp.1;


                match keyf {
                    6 => {
                        (result1, abserr1, rounderr1) = qk61.integrate(f, a1, b1);
                        (result2, abserr2, rounderr2) = qk61.integrate(f, a2, b2);
                    },
                    _ => (),
                }

                res_update(&mut result, & result1, & result2, &comp.3);

                interval_cache.insert((Myf64{x:a1},Myf64{x:b1}),result1);
                interval_cache.insert((Myf64{x:a2},Myf64{x:b2}),result2);

                heap.push(HeapItem::new((a1,b1),abserr1));
                heap.push(HeapItem::new((a2,b2),abserr2));


                abserr += -comp.2 + abserr1 + abserr2;
                rounderr += rounderr1 + rounderr2;

            }
            errbnd = epsabs.max(epsrel * norm_vec(&result));


            if abserr <= errbnd / 8.0{ break;}
            if abserr < rounderr {
                return QagVecNormIntegratorResult::new_error(ResultState::BadTolerance)
             }
        }


        if last >= self.limit {
            return QagVecNormIntegratorResult::new_error(ResultState::MaxIteration)
        }

        return QagVecNormIntegratorResult::new(result,abserr,neval,last)

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
    use crate::qage_vec_norm2::QagVecNorm2;
    use crate::qage_vec_norm3::QagVecNorm3;
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
        let ff = |x:f64| [f1(x),f2(x),f3(x),f4(x)];


        let fv = |x:f64| [f1(x),f2(x),f3(x),f4(x)];

        let a = 0.0;
        let b = 1000000000.0;
        let key = 6;
        let epsabs = 1.0e-3;
        let epsrel = 0.0;
        let limit = 10000000;
        let max = 2;
        let qag_vec_norm = QagVecNorm{key,limit};
        let qag_vec_norm2 = QagVecNorm2{key,limit};
        let qag_vec_norm3 = QagVecNorm3{key,limit};
        let qag_vec_nosort = QagVecNosort{key,limit};


        let mut res;
        let mut res_vec;
        let mut res2;
        let mut res3;



        for k in 0..max {
            let start = Instant::now();
            res = qag_vec_norm.qintegrate(&f,a,b,epsabs,epsrel);
            println!("scipy time : {:?}",start.elapsed());
            let start = Instant::now();
            res_vec = qag_vec_nosort.qintegrate(&fv, a, b,epsabs,epsrel);
            println!("my time : {:?}",start.elapsed());
            let start = Instant::now();
            res2= qag_vec_norm2.qintegrate(&ff,a,b,epsabs,epsrel);
            println!("scipy time2 : {:?}",start.elapsed());
            let start = Instant::now();
            res3= qag_vec_norm3.qintegrate(&ff,a,b,epsabs,epsrel);
            println!("scipy time3 : {:?}",start.elapsed());

            if k == max-1{
                println!("{:?}",res);
                println!("{:?}",res2.unwrap());
                println!("{:?}",res3.unwrap());
                println!("{:?}",res_vec.integration_result.result);
                println!("{:?}",res_vec.integration_result.abserr);
                println!("{:?}",res_vec.integration_result.last);
            }

        }

    }
}


