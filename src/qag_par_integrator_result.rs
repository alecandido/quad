use crate::qag_par_integration_result::QagParIntegrationResult;
use crate::result_state::*;

#[derive(Clone,Debug)]
pub struct QagParIntegratorResult<const N:usize> {
    pub result_state : ResultState,
    pub integration_result : QagParIntegrationResult<N>,
}



impl<const N:usize> QagParIntegratorResult<N> {
    pub fn new(result : [f64; N], abserr : f64, neval : i32, last : usize) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagParIntegrationResult::new(result, abserr, neval, last),
        }
    }
    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QagParIntegrationResult::new_error(),
        }
    }


    pub fn unwrap(&self) -> QagParIntegrationResult<N> {
        match self.result_state{
            ResultState::Success => self.integration_result.clone(),
            ResultState::Failure => panic!("Generic Fail"),
            ResultState::MaxIteration => panic!("{}", MAX_ITERATION_ERROR_MESSAGE),
            ResultState::BadTolerance => panic!("{}", BAD_TOLERANCE_ERROR_MESSAGE),
            ResultState::Invalid => panic!("{}", INVALID_ERROR_MESSAGE),
            ResultState::BadFunction => panic!("{}", BAD_FUNCTION_ERROR_MESSAGE),
            ResultState::Diverge => panic!("{}", DIVERGE_ERROR_MESSAGE),

        }
    }
}
