use std::f64::consts::PI;

pub fn rastrigin(x : f64 ) -> f64 {
    10.0 + x.powi(2) - 10.0 * ( 2.0 * PI * x ).cos()
}