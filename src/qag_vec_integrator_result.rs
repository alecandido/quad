use std::collections::{BinaryHeap, HashMap};
use crate::constants::{HeapItem, Myf64};
use crate::qag_vec_integration_result::QagVecIntegrationResult;
use crate::result_state::*;

#[derive(Clone,Debug)]
pub struct QagVecIntegratorResult {
    pub result_state : ResultState,
    pub integration_result : QagVecIntegrationResult,
}



impl QagVecIntegratorResult{
    pub fn new(result : Vec<f64>, abserr : f64) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagVecIntegrationResult::new(result, abserr),
        }
    }

    pub fn new_more_info(result : Vec<f64>, abserr : f64, neval : i32, last : usize,
                         hash : HashMap<(Myf64,Myf64),Vec<f64>>, heap : BinaryHeap<HeapItem>) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagVecIntegrationResult::new_more_info(result, abserr, neval, last, hash, heap)
        }
    }

    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QagVecIntegrationResult::new_error(),
        }
    }


    pub fn unwrap(&self) -> QagVecIntegrationResult {
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
