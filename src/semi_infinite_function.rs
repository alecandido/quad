use crate::constants::*;

pub fn SemiInfiniteFunction<const N : usize,F>(f : &F,x : f64, start : f64, infty : f64) -> [f64;N]
where F : Fn(f64) -> [f64;N] {
    if x < UFLOW.sqrt() { return [0.0;N] }
    let sgn = if infty.is_sign_positive() { 1.0 } else { - 1.0};
    let z = start + sgn * (1.0 - x) / x;
    let mut res: [f64;N] = f(z);
    div( &mut res, sgn * x * x);
    res
}

pub fn DoubleInfiniteFunction<const N : usize,F>(f : &F,x : f64) -> [f64;N]
    where F : Fn(f64) -> [f64;N] {
    if x.abs() < UFLOW.sqrt() { return [0.0;N] }
    let z = (1.0 - x.abs()) / x;
    let mut res: [f64;N] = f(z);
    div( &mut res,  x * x);
    res



}




pub fn div( arr : &mut [f64], scalar : f64){
    for k in 0..arr.len(){
        arr[k] /= scalar ;
    }
}