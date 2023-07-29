#[cfg(doc)]
use crate::errors::QagError;
use crate::qag::Qag;

use ndarray::Array1;
use std::cmp::Ordering;
use std::hash;
use std::sync::Arc;
/// Vector of function.
#[derive(Clone)]
pub struct FnVec<'a> {
    pub components: Arc<dyn Fn(f64) -> Array1<f64> + Send + Sync + 'a>,
}
/// [Machine epsilon] value for `f64`.
///
/// This is the difference between `1.0` and the next larger representable number.
///
/// [Machine epsilon]: https://en.wikipedia.org/wiki/Machine_epsilon
pub const EPMACH: f64 = f64::EPSILON;
/// Smallest positive normal `f64` value.
pub const UFLOW: f64 = f64::MIN_POSITIVE;
/// Parameter of [iroff1_flag].
pub const IROFF_PARAMETER1: f64 = 0.00001;
/// Parameter of [iroff1_flag].
pub const IROFF_PARAMETER2: f64 = 0.99;
/// Threshold for 'iroff1' beyond which the [BadTolerance](QagError::BadTolerance) is set.
pub const IROFF1_THRESHOLD: i32 = 6;
/// Threshold for 'iroff1' beyond which the [BadTolerance](QagError::BadTolerance) is set.
pub const IROFF2_THRESHOLD: i32 = 20;
/// Parameter of [bad_function_flag].
pub const BAD_FUNCTION_PARAMETER1: f64 = 100.0;
/// Parameter of [bad_function_flag].
pub const BAD_FUNCTION_PARAMETER2: f64 = 1000.0;
/// Norm of an [Array1].
pub fn norm_ar(ar: &Array1<f64>) -> f64 {
    ar.iter().map(|x| x.powi(2)).sum::<f64>().sqrt()
}
/// Transform the list of additional points in case of semi-infinite or infinite interval.
pub fn points_transformed(mut points: Vec<f64>, a: f64, b: f64) -> Vec<f64> {
    points.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut points_transformed = vec![0.0; 0];
    for point in &points {
        points_transformed.push(if b == f64::INFINITY && a.is_finite() {
            1.0 / (*point - a + 1.0)
        } else if a == f64::NEG_INFINITY && b.is_finite() {
            1.0 / (b - *point + 1.0)
        } else {
            point.signum() / (point.abs() + 1.0)
        });
    }
    points_transformed
}
/// Condition to increase iroff1.
pub fn iroff1_flag(
    old_res: &Array1<f64>,
    new_res: &Array1<f64>,
    new_abserr: f64,
    old_abserr: f64,
) -> bool {
    for k in 0..old_res.len() {
        if !((old_res[k] - new_res[k]).abs() <= IROFF_PARAMETER1 * new_res[k].abs()
            && new_abserr >= IROFF_PARAMETER2 * old_abserr)
        {
            return false;
        }
    }
    return true;
}
/// Condition to return a [BadFunction](QagError::BadFunction) .
pub fn bad_function_flag(x: f64, y: f64) -> bool {
    if x.abs().max(y.abs())
        <= (1.0 + BAD_FUNCTION_PARAMETER1 * EPMACH)
            * (((x + y) / 2.0).abs() + BAD_FUNCTION_PARAMETER2 * UFLOW)
    {
        return true;
    }
    false
}
/// Heap used in [qintegrate](Qag::qintegrate) to store the sub-intervals and their errors.
#[derive(Debug, Clone)]
pub struct HeapItem {
    pub interval: (f64, f64),
    pub err: f64,
}

impl HeapItem {
    pub fn new(interval: (f64, f64), err: f64) -> Self {
        Self { interval, err }
    }
}

impl Eq for HeapItem {}

impl PartialEq for HeapItem {
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
/// 'f64' implementing Hash.
///
/// Needed to used an interval as key in a [HashMap](std::collections::HashMap).
#[derive(Debug, Clone)]
pub struct Myf64 {
    pub x: f64,
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

impl Eq for Myf64 {}
