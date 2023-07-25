use crate::constants::*;
use ndarray::Array1;

pub fn semi_infinite_function<F>(f: &F, x: f64, start: f64, infty: f64) -> Array1<f64>
where
    F: Fn(f64) -> Array1<f64> + ?Sized,
{
    if x < UFLOW.sqrt() {
        return Array1::<f64>::from_vec(vec![0.0; f(0.0).len()]);
    }
    let sgn = if infty.is_sign_positive() { 1.0 } else { -1.0 };
    let z = start + sgn * (1.0 - x) / x;
    let mut res: Array1<f64> = f(z);
    div(&mut res, sgn * x * x);
    res
}

pub fn double_infinite_function<F>(f: &F, x: f64) -> Array1<f64>
where
    F: Fn(f64) -> Array1<f64> + ?Sized,
{
    if x.abs() < UFLOW.sqrt() {
        return Array1::<f64>::from_vec(vec![0.0; f(0.0).len()]);
    }
    let z = (1.0 - x.abs()) / x;
    let mut res: Array1<f64> = f(z);
    div(&mut res, x * x);
    res
}

pub fn div(arr: &mut Array1<f64>, scalar: f64) {
    for k in 0..arr.len() {
        arr[k] /= scalar;
    }
}
