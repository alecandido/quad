use crate::qag_vec_norm_integration_result::QagVecNormIntegrationResult;
use crate::result_state::*;

#[derive(Clone,Debug)]
pub struct QagVecNormIntegratorResult<const N:usize> {
    pub result_state : ResultState,
    pub integration_result : QagVecNormIntegrationResult<N>,
}



impl<const N:usize> QagVecNormIntegratorResult<N> {
    pub fn new(result : [f64; N], abserr : f64, neval : i32, last : usize) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagVecNormIntegrationResult::new(result, abserr, neval,last),
        }
    }
    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QagVecNormIntegrationResult::new_error(),
        }
    }


    pub fn unwrap(&self) -> QagVecNormIntegrationResult<N> {
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
