use std::simd::Simd;
use crate::qag_vec8_integration_result::{QagVec8IntegrationResult, ResultVec8};
use crate::qag_vec_integration_result::*;
use crate::result_state::*;
use crate::qage_vec::*;

#[derive(Clone,Debug)]
pub struct QagVec8IntegratorResult {
    pub result_state : ResultState,
    pub integration_result : QagVec8IntegrationResult,
}



impl QagVec8IntegratorResult {
    pub fn new(result : Simd<f64,8>, abserr : Simd<f64,8>, neval : i32, list : Vec<ResultVec8>,
               last : usize) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagVec8IntegrationResult::new(result, abserr, neval, list, last ),
        }
    }
    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QagVec8IntegrationResult::new_error(),
        }
    }


    pub fn unwrap(&self) -> QagVec8IntegrationResult {
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
