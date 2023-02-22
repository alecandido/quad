use crate::qag_integration_result::*;

#[derive(Clone,Debug)]
pub struct QagIntegratorResult{
    result_state : ResultState,
    integration_result : QagIntegrationResult,
}

const MAX_ITERATION_ERROR_MESSAGE : &str = "Maximum number of subdivisions allowed has been achieved.
            One can allow more subdivisions by increasing the value of limit. However, if this
            yields no improvement it is rather advised to analyze the integrand in order to
            determine the integration difficulties. If the position of a local difficulty can be
            determined(e.g. singularity, discontinuity within the interval) one will probably gain
            from splitting up the interval at this point and calling the integrator on the
            subranges. If possible, an appropriate special-purpose integrator should be used
            which is designed for handling the type of difficulty involved.";
const BAD_TOLERANCE_ERROR_MESSAGE : &str = "The occurrence of roundoff error is detected, which \
            prevents the requested tolerance from being achieved.";
const INVALID_ERROR_MESSAGE : &str = "The input is invalid, because epsabs <= 0 and\
            epsrel < max(50 * rel.mach.acc.,0.5d-28)";
const BAD_FUNCTION_ERROR_MESSAGE : &str = "Extremely bad integrand behaviour occurs at some points\
            of the integration interval.";


impl QagIntegratorResult{
    pub fn new(result : f64, abserr : f64, neval : i32, alist : Vec<f64>,
               blist : Vec<f64>, rlist : Vec<f64>, elist : Vec<f64>, iord : Vec<usize>,
               last : usize) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagIntegrationResult::new(result , abserr , neval , alist ,
                                                           blist , rlist , elist , iord , last ),
        }
    }
    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QagIntegrationResult::new_error(),
        }
    }


    pub fn unwrap(&self) -> QagIntegrationResult {
        match self.result_state{
            ResultState::Success => self.integration_result.clone(),
            ResultState::Failure => panic!("Generic Fail"),
            ResultState::MaxIteration => panic!("{}", MAX_ITERATION_ERROR_MESSAGE),
            ResultState::BadTolerance => panic!("{}", BAD_TOLERANCE_ERROR_MESSAGE),
            ResultState::Invalid => panic!("{}", INVALID_ERROR_MESSAGE),
            ResultState::BadFunction => panic!("{}", BAD_FUNCTION_ERROR_MESSAGE),

        }
    }
}
