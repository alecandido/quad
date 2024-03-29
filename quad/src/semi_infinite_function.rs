use crate::constants::*;
use ndarray::Array1;
/// Transform the function in case of semi-infinite interval.
///
/// For an interval (start,+∞) integrand is transformed using the transformation x = start + (1-t)/t.
/// For an interval (-∞,start) integrand is transformed using the transformation x = start - (1-t)/t.
pub fn semi_infinite_function<F>(f: &F, x: f64, start: f64, infty: f64) -> Array1<f64>
where
    F: Fn(f64) -> Array1<f64> + ?Sized,
{
    if x < UFLOW.sqrt() {
        return Array1::<f64>::zeros(f(0.0).len());
    }
    let sgn = if infty.is_sign_positive() { 1.0 } else { -1.0 };
    let z = start + sgn * (1.0 - x) / x;
    let res: Array1<f64> = f(z);
    res / (sgn * x * x)
}
/// Transform the function in case of infinite interval.
///
/// For an interval (-∞,+∞) integrand is transformed using the transformation x = (1-t)/t.
pub fn double_infinite_function<F>(f: &F, x: f64) -> Array1<f64>
where
    F: Fn(f64) -> Array1<f64> + ?Sized,
{
    if x.abs() < UFLOW.sqrt() {
        return Array1::<f64>::zeros(f(0.0).len());
    }
    let z = (1.0 - x.abs()) / x;
    let res: Array1<f64> = f(z);
    res / (x * x)
}
