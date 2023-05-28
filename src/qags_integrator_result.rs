use std::collections::{BinaryHeap, HashMap};
use crate::constants::{HeapItem, Myf64};
use crate::qags_integration_result::QagsIntegrationResult;
use crate::result_state::*;

#[derive(Clone,Debug)]
pub struct QagsIntegratorResult<const N:usize> {
    pub result_state : ResultState,
    pub integration_result : QagsIntegrationResult<N>,
}



impl<const N:usize> QagsIntegratorResult<N> {
    pub fn new(result : [f64; N], abserr : f64) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagsIntegrationResult::new(result, abserr),
        }
    }

    pub fn new_more_info(result : [f64; N], abserr : f64, neval : i32, last : usize,
                         hash : HashMap<(Myf64,Myf64),[f64;N]>, heap : BinaryHeap<HeapItem>) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagsIntegrationResult::new_more_info(result,abserr,neval,last,hash,heap)
        }
    }

    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QagsIntegrationResult::new_error(),
        }
    }


    pub fn unwrap(&self) -> QagsIntegrationResult<N> {
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
