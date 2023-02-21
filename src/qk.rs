pub trait Qk{
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64,) -> (f64, f64, f64, f64);
}

pub const EPMACH : f64 = f64::EPSILON;
pub const UFLOW : f64 = f64::MIN_POSITIVE;