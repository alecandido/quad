use std::simd::Simd;
use crate::qag_vec4_integration_result::{QagVec4IntegrationResult, ResultVec4};
use crate::qag_vec_integration_result_pre::*;
use crate::result_state::*;
use crate::qage_vec_pre::*;

#[derive(Clone,Debug)]
pub struct QagVec4IntegratorResult {
    pub result_state : ResultState,
    pub integration_result : QagVec4IntegrationResult,
}



impl QagVec4IntegratorResult {
    pub fn new(result : Simd<f64,4>, abserr : Simd<f64,4>, neval : i32, list : Vec<ResultVec4>,
               last : usize) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagVec4IntegrationResult::new(result, abserr, neval, list, last ),
        }
    }
    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QagVec4IntegrationResult::new_error(),
        }
    }


    pub fn unwrap(&self) -> QagVec4IntegrationResult {
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
