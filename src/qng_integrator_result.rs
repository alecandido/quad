use crate::qng_integration_result::*;
use crate::result_state::*;

#[derive(Clone,Debug)]
pub struct QngIntegratorResult{
    pub result_state : ResultState,
    pub integration_result : QngIntegrationResult,
}

impl QngIntegratorResult{
    pub fn new(result : f64, abserr : f64, neval : i32) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QngIntegrationResult::new(result , abserr , neval),
        }
    }
    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QngIntegrationResult::new_error(),
        }
    }


    pub fn unwrap(&self) -> QngIntegrationResult {
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