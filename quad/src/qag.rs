use ::rayon::prelude::*;

use crate::constants::*;
use crate::errors::QagError;
use crate::qag_integration_result::QagIntegrationResult;
use crate::qk15::qk15_quadrature;
use crate::qk21::qk21_quadrature;
use crate::qk31::qk31_quadrature;
use crate::qk41::qk41_quadrature;
use crate::qk51::qk51_quadrature;
use crate::qk61::qk61_quadrature;
use crate::semi_infinite_function::{double_infinite_function, semi_infinite_function};
use ndarray::Array1;
use std::collections::{BinaryHeap, HashMap};
use std::sync::Arc;

#[derive(Clone)]
pub struct Qag {
    pub key: i32,
    pub limit: usize,
    pub points: Vec<f64>,
    pub number_of_thread: usize,
    pub more_info: bool,
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

impl Qag {
    pub fn integrate(
        &self,
        fun: &FnVec,
        a: f64,
        b: f64,
        epsabs: f64,
        epsrel: f64,
    ) -> Result<QagIntegrationResult, QagError> {
        let f = &fun.components;
        if b == f64::INFINITY && a.is_finite()
            || a == f64::NEG_INFINITY && b.is_finite()
            || a == f64::NEG_INFINITY && b == f64::INFINITY
        {
            let points = points_transformed(self.points.clone(), a, b);
            let qag = Qag {
                key: self.key,
                limit: self.limit,
                points,
                number_of_thread: self.number_of_thread,
                more_info: self.more_info,
            };

            if b == f64::INFINITY && a.is_finite() {
                let f2 = FnVec {
                    components: Arc::new(|x: f64| semi_infinite_function(&**f, x, a, b)),
                };
                return qag.qintegrate(&f2, 0.0, 1.0, epsabs, epsrel);
            } else if a == f64::NEG_INFINITY && b.is_finite() {
                let f2 = FnVec {
                    components: Arc::new(|x: f64| semi_infinite_function(&**f, x, b, a)),
                };
                return qag.qintegrate(&f2, 0.0, 1.0, epsabs, epsrel);
            } else if a == f64::NEG_INFINITY && b == f64::INFINITY {
                let f2 = FnVec {
                    components: Arc::new(|x: f64| double_infinite_function(&**f, x)),
                };
                return qag.qintegrate(&f2, -1.0, 1.0, epsabs, epsrel);
            };
        }

        self.qintegrate(&fun, a, b, epsabs, epsrel)
    }

    pub fn qintegrate(
        &self,
        fun: &FnVec,
        a: f64,
        b: f64,
        epsabs: f64,
        epsrel: f64,
    ) -> Result<QagIntegrationResult, QagError> {
        if epsabs <= 0.0 && epsrel < 0.5e-28_f64.max(50.0 * EPMACH) {
            return Err(QagError::Invalid);
        }

        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.number_of_thread)
            .build()
            .unwrap();

        let mut initial_intervals = vec![];
        let mut points = self.points.clone();
        points.sort_by(|a, b| a.partial_cmp(b).unwrap());

        if points.is_empty() {
            initial_intervals.push((a, b));
        } else {
            let mut prev = a;
            for p in points {
                if p > a && p < b {
                    initial_intervals.push((prev, p));
                    prev = p;
                }
            }
            initial_intervals.push((prev, b));
        }

        let f = &fun.components;
        let n: usize = f(0.0).len();
        let mut neval = 0;
        let mut last = 1;
        let mut interval_cache = HashMap::new();
        let mut heap = BinaryHeap::new();
        let mut result = Array1::<f64>::zeros(n);
        let mut abserr = 0.0;
        let mut rounderr = 0.0;
        let mut iroff1 = 0;
        let mut iroff2 = 0;
        let mut keyf = self.key;
        if self.key <= 0 {
            keyf = 1;
        }
        if self.key >= 7 {
            keyf = 6;
        }

        for comp in initial_intervals {
            let (result_temp, abserr_temp, rounderr_temp) = match keyf {
                1 => qk15_quadrature(&**f, comp.0, comp.1),
                2 => qk21_quadrature(&**f, comp.0, comp.1),
                3 => qk31_quadrature(&**f, comp.0, comp.1),
                4 => qk41_quadrature(&**f, comp.0, comp.1),
                5 => qk51_quadrature(&**f, comp.0, comp.1),
                6 => qk61_quadrature(&**f, comp.0, comp.1),
                _ => (Array1::<f64>::from_vec(vec![0.0; f(0.0).len()]), 0.0, 0.0),
            };
            result += &(Array1::<f64>::from(result_temp.clone()));
            abserr += abserr_temp;
            rounderr += rounderr_temp;
            heap.push(HeapItem::new((comp.0, comp.1), abserr_temp));
            interval_cache.insert((Myf64 { x: comp.0 }, Myf64 { x: comp.1 }), result_temp);
        }

        let mut errbnd = epsabs.max(epsrel * norm_ar(&result));

        if abserr + rounderr <= errbnd {
            if keyf != 1 {
                neval = (10 * keyf + 1) * (2 * last as i32 - 1);
            }
            if keyf == 1 {
                neval = 30 * last as i32 + 15;
            }
            abserr = abserr + rounderr;
            if self.more_info {
                return Ok(QagIntegrationResult::new_more_info(
                    result,
                    abserr,
                    neval,
                    last,
                    interval_cache,
                    heap,
                ));
            } else {
                return Ok(QagIntegrationResult::new(result, abserr));
            }
        }

        if self.limit == 1 {
            return Err(QagError::MaxIteration);
        }

        if abserr < rounderr {
            return Err(QagError::BadTolerance);
        }

        while last < self.limit {
            let mut to_process = vec![];
            let mut err_sum = 0.0;
            let mut old_result = Array1::<f64>::zeros(n);
            let max_new_divison = self.limit - last;

            while to_process.len() < 128.min(max_new_divison) && heap.len() != 0 {
                let old_interval = heap.pop().unwrap();
                let ((x, y), old_err) = (old_interval.interval, old_interval.err);
                if bad_function_flag(x, y) {
                    return Err(QagError::BadFunction);
                }
                let old_res = interval_cache
                    .remove(&(Myf64 { x }, Myf64 { x: y }))
                    .unwrap();
                err_sum += old_err;
                old_result += &Array1::<f64>::from(old_res);
                to_process.push((x, y));
                if err_sum > abserr - errbnd / 8.0 {
                    break;
                }
            }

            last += to_process.len();

            let new_result: (Vec<_>, Vec<_>) = pool.install(|| {
                to_process
                    .par_iter()
                    .map(|comp| {
                        let mut result1 = Array1::<f64>::from_elem(1, 0.0);
                        let mut abserr1 = 0.0;
                        let mut rounderr1 = 0.0;

                        let mut result2 = Array1::<f64>::from_elem(1, 0.0);
                        let mut abserr2 = 0.0;
                        let mut rounderr2 = 0.0;

                        let a1 = comp.0;
                        let b1 = 0.5 * (comp.0 + comp.1);
                        let a2 = b1;
                        let b2 = comp.1;

                        match keyf {
                            1 => {
                                (result1, abserr1, rounderr1) = qk15_quadrature(&**f, a1, b1);
                                (result2, abserr2, rounderr2) = qk15_quadrature(&**f, a2, b2);
                            }
                            2 => {
                                (result1, abserr1, rounderr1) = qk21_quadrature(&**f, a1, b1);
                                (result2, abserr2, rounderr2) = qk21_quadrature(&**f, a2, b2);
                            }
                            3 => {
                                (result1, abserr1, rounderr1) = qk31_quadrature(&**f, a1, b1);
                                (result2, abserr2, rounderr2) = qk31_quadrature(&**f, a2, b2);
                            }
                            4 => {
                                (result1, abserr1, rounderr1) = qk41_quadrature(&**f, a1, b1);
                                (result2, abserr2, rounderr2) = qk41_quadrature(&**f, a2, b2);
                            }
                            5 => {
                                (result1, abserr1, rounderr1) = qk51_quadrature(&**f, a1, b1);
                                (result2, abserr2, rounderr2) = qk51_quadrature(&**f, a2, b2);
                            }
                            6 => {
                                (result1, abserr1, rounderr1) = qk61_quadrature(&**f, a1, b1);
                                (result2, abserr2, rounderr2) = qk61_quadrature(&**f, a2, b2);
                            }
                            _ => (),
                        }
                        (
                            (a1, b1, result1, abserr1, rounderr1),
                            (a2, b2, result2, abserr2, rounderr2),
                        )
                    })
                    .collect()
            });

            let mut new_res = Array1::<f64>::zeros(n);
            let mut new_abserr = 0.0;

            for k in 0..new_result.0.len() {
                new_res += &(Array1::<f64>::from(new_result.0[k].2.clone()));
                new_res += &(Array1::<f64>::from(new_result.1[k].2.clone()));
                new_abserr += new_result.0[k].3 + new_result.1[k].3;
                rounderr += new_result.0[k].4 + new_result.1[k].4;
                interval_cache.insert(
                    (
                        Myf64 {
                            x: new_result.0[k].0,
                        },
                        Myf64 {
                            x: new_result.0[k].1,
                        },
                    ),
                    new_result.0[k].2.clone(),
                );
                interval_cache.insert(
                    (
                        Myf64 {
                            x: new_result.1[k].0,
                        },
                        Myf64 {
                            x: new_result.1[k].1,
                        },
                    ),
                    new_result.1[k].2.clone(),
                );
                heap.push(HeapItem::new(
                    (new_result.0[k].0, new_result.0[k].1),
                    new_result.0[k].3,
                ));
                heap.push(HeapItem::new(
                    (new_result.1[k].0, new_result.1[k].1),
                    new_result.1[k].3,
                ));
            }
            if iroff1_flag(&old_result, &new_res, new_abserr, err_sum) {
                iroff1 += 1;
            }
            if last > 10 && new_abserr > err_sum {
                iroff2 += 1;
            }
            result += &new_res;
            result -= &old_result;
            abserr += new_abserr - err_sum;

            errbnd = epsabs.max(epsrel * norm_ar(&result));

            if abserr <= errbnd / 8.0 {
                break;
            }
            if abserr < rounderr || iroff1 >= IROFF1_THRESHOLD || iroff2 >= IROFF2_THRESHOLD {
                return Err(QagError::BadTolerance);
            }
        }

        if abserr > errbnd / 8.0 && last >= self.limit {
            return Err(QagError::MaxIteration);
        }

        if keyf != 1 {
            neval = (10 * keyf + 1) * (2 * last as i32 - 1);
        }
        if keyf == 1 {
            neval = 30 * last as i32 + 15;
        }

        abserr = abserr + rounderr;

        if self.more_info {
            return Ok(QagIntegrationResult::new_more_info(
                result,
                abserr,
                neval,
                last,
                interval_cache,
                heap,
            ));
        } else {
            return Ok(QagIntegrationResult::new(result, abserr));
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::constants::{FnVec, Myf64};
    use crate::errors::QagError;
    use crate::qag::Qag;
    use ndarray::array;
    use std::sync::Arc;

    #[test]
    fn max_iteration1() {
        let a = 0.0;
        let b = 10000.0;
        let epsrel = 0.0;
        let epsabs = 1.0e-2;
        let limit = 1;
        let key = 6;

        let qag = Qag {
            key,
            limit,
            points: vec![0.0; 0],
            number_of_thread: 8,
            more_info: true,
        };

        let f = FnVec {
            components: Arc::new(|x: f64| array![x.sin(), x.cos()]),
        };
        let res = qag.integrate(&f, a, b, epsabs, epsrel);
        let error = res.unwrap_err();

        assert_eq!(error, QagError::MaxIteration);
    }
    #[test]
    fn max_iteration2() {
        let a = 0.0;
        let b = 1000000.0;
        let epsrel = 0.0;
        let epsabs = 1.0e-2;
        let limit = 30;
        let key = 6;

        let qag = Qag {
            key,
            limit,
            points: vec![0.0; 0],
            number_of_thread: 8,
            more_info: true,
        };

        let f = FnVec {
            components: Arc::new(|x: f64| array![x.sin(), x.cos()]),
        };
        let res = qag.integrate(&f, a, b, epsabs, epsrel);
        let error = res.unwrap_err();

        assert_eq!(error, QagError::MaxIteration);
    }

    #[test]
    fn invalid() {
        let a = 0.0;
        let b = 1000000.0;
        let epsrel = 1.0e-30;
        let epsabs = 0.0;
        let limit = 30;
        let key = 6;

        let qag = Qag {
            key,
            limit,
            points: vec![0.0; 0],
            number_of_thread: 8,
            more_info: true,
        };

        let f = FnVec {
            components: Arc::new(|x: f64| array![x.sin(), x.cos()]),
        };
        let res = qag.integrate(&f, a, b, epsabs, epsrel);
        let error = res.unwrap_err();

        assert_eq!(error, QagError::Invalid);
    }

    #[test]
    fn key() {
        let a = 0.0;
        let b = 10000.0;
        let epsrel = 0.0;
        let epsabs = 1.0e-3;
        let limit = 10000;
        let correct_result = [1.0 - 10000.0_f64.cos(), 10000.0_f64.sin()];

        for key in 1..7 {
            let qag = Qag {
                key,
                limit,
                points: vec![0.0; 0],
                number_of_thread: 8,
                more_info: true,
            };

            let f = FnVec {
                components: Arc::new(|x: f64| array![x.sin(), x.cos()]),
            };
            let res = qag.integrate(&f, a, b, epsabs, epsrel).unwrap();

            assert!(
                res.result[0] - correct_result[0] < epsabs
                    && res.result[1] - correct_result[1] < epsabs
            );
        }
    }
    #[test]
    fn semi_infinite() {
        let a = 0.0;
        let b = f64::INFINITY;
        let c = f64::NEG_INFINITY;
        let epsrel = 0.0;
        let epsabs = 1.0e-12;
        let limit = 10000;
        let key = 6;
        let correct_result = [0.4, 0.6];

        let qag = Qag {
            key,
            limit,
            points: vec![0.0; 0],
            number_of_thread: 8,
            more_info: true,
        };

        let f = FnVec {
            components: Arc::new(|x: f64| {
                array![
                    x.sin().powi(2) / x.abs().exp(),
                    x.cos().powi(2) / x.abs().exp(),
                ]
            }),
        };

        let res1 = qag.integrate(&f, a, b, epsabs, epsrel).unwrap();
        let res2 = qag.integrate(&f, c, a, epsabs, epsrel).unwrap();

        assert!(
            res1.result[0] - correct_result[0] < epsabs
                && res1.result[1] - correct_result[1] < epsabs
        );
        assert!(
            res2.result[0] - correct_result[0] < epsabs
                && res2.result[1] - correct_result[1] < epsabs
        );
    }
    #[test]
    fn double_infinite() {
        let a = f64::NEG_INFINITY;
        let b = f64::INFINITY;
        let epsrel = 0.0;
        let epsabs = 1.0e-10;
        let limit = 10000;
        let key = 6;
        let correct_result = [1.2879903316984565533522585284072106913, 1.5974];

        let qag = Qag {
            key,
            limit,
            points: vec![0.0; 0],
            number_of_thread: 8,
            more_info: true,
        };

        let f = FnVec {
            components: Arc::new(|x: f64| {
                array![
                    x.sin().powi(2) / x.abs().exp2(),
                    x.cos().powi(2) / x.abs().exp2(),
                ]
            }),
        };

        let res = qag.integrate(&f, a, b, epsabs, epsrel).unwrap();
        assert!(
            res.result[0] - correct_result[0] < epsabs
                && res.result[1] - correct_result[1] < epsabs
        );
    }
    #[test]
    fn additional_points() {
        let a = 0.0;
        let b = 1.0;
        let epsrel = 0.0;
        let epsabs = 1.0;
        let limit = 10000;
        let key = 6;
        let points = vec![0.0, 0.2, 0.4, 0.6, 0.8, 1.0];

        let qag = Qag {
            key,
            limit,
            points: points.clone(),
            number_of_thread: 8,
            more_info: true,
        };
        let f = FnVec {
            components: Arc::new(|x: f64| array![x.cos(), x.sin()]),
        };
        let res = qag.integrate(&f, a, b, epsabs, epsrel).unwrap();
        let mut res_hash = res.more_info.unwrap().hash.clone();
        assert_eq!(res_hash.len(), qag.points.len() - 1);
        for k in 0..points.len() - 1 {
            res_hash.remove(&((Myf64 { x: points[k] }, Myf64 { x: points[k + 1] })));
        }
        assert_eq!(res_hash.len(), 0);
    }
}
