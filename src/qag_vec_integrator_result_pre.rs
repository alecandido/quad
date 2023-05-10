use crate::qag_vec_integration_result_pre::*;
use crate::result_state::*;
use crate::qage_vec_pre::*;

#[derive(Clone,Debug)]
pub struct QagVecIntegratorResultPre {
    pub result_state : ResultState,
    pub integration_result : QagVecIntegrationResultPre,
}



impl QagVecIntegratorResultPre {
    pub fn new(result : Vec<f64>, abserr : Vec<f64>, neval : i32, list : Vec<ResultVecPre>,
               last : usize) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QagVecIntegrationResultPre::new(result, abserr, neval, list, last ),
        }
    }
    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QagVecIntegrationResultPre::new_error(),
        }
    }


    pub fn unwrap(&self) -> QagVecIntegrationResultPre {
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
