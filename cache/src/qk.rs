use crate::funct_vector::FnVec;

pub trait Qk{
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64,) -> (f64, f64, f64, f64);
}

pub trait QkVec{
    fn integrate(&self,f : &FnVec, a : f64, b : f64,) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>);
}

pub const EPMACH : f64 = f64::EPSILON;          // the largest relative spacing.
pub const UFLOW : f64 = f64::MIN_POSITIVE;      // the smallest positive magnitude.
pub const OFLOW : f64 = f64::MAX;               // oflow is the largest positive magnitude.

