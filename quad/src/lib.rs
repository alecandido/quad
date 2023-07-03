pub mod constants;
mod errors;
pub mod qag;
mod qag_integration_result;
mod qag_par;
mod qk;
pub mod qk15;
pub mod qk21;
pub mod qk31;
pub mod qk41;
pub mod qk51;
pub mod qk61;
pub mod semi_infinite_function;

use crate::constants::FnVec;
use crate::errors::QagError;
use crate::qag::Qag;
use crate::qag_integration_result::QagIntegrationResult;
use crate::qag_par::QagPar;

pub fn integrate<F>(
    f: &F,
    a: f64,
    b: f64,
    epsabs: f64,
    epsrel: f64,
    key: i32,
    limit: usize,
    points: Vec<f64>,
    more_info: bool,
) -> Result<QagIntegrationResult, QagError>
where
    F: Fn(f64) -> Vec<f64>,
{
    let qag = Qag {
        key,
        limit,
        points,
        more_info,
    };
    qag.integrate(&f, a, b, epsabs, epsrel)
}

pub fn integrate_par(
    f: &FnVec,
    a: f64,
    b: f64,
    epsabs: f64,
    epsrel: f64,
    key: i32,
    limit: usize,
    points: Vec<f64>,
    number_of_thread: usize,
    more_info: bool,
) -> Result<QagIntegrationResult, QagError> {
    let qag = QagPar {
        key,
        limit,
        points,
        number_of_thread,
        more_info,
    };
    qag.qintegrate(&f, a, b, epsabs, epsrel)
}
