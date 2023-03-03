use crate::qag_integration_result::*;
use crate::result_state::*;

#[derive(Clone,Debug)]
pub struct QagIntegratorResult{
    pub result_state : ResultState,
    pub integration_result : QagIntegrationResult,
}



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
            ResultState::Diverge => panic!("{}", DIVERGE_ERROR_MESSAGE),

        }
    }
}
