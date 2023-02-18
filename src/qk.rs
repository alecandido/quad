pub trait Qk{
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64,) -> (f64, f64, f64, f64);
}