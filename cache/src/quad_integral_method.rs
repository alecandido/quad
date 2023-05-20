use crate::quad_integrator_result::QuadIntegratorResult;

pub trait QuadIntegralMethod{
    fn integrate(&self, fun : &dyn Fn(f64)->f64, a : f64, b : f64, epsabs: f64, epsrel : f64) -> QuadIntegratorResult;
}