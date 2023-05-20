use crate::qag_integrator_result::QagIntegratorResult;
use crate::qag_1dvec_integrator_result::Qag1DVecIntegratorResult;
use crate::qng_integrator_result::QngIntegratorResult;
use crate::quad_integration_result::*;
use crate::result_state::*;

#[derive(Clone,Debug)]
pub struct QuadIntegratorResult{
    result_state : ResultState,
    integration_result : QuadIntegrationResult,
}



impl QuadIntegratorResult{
    pub fn new(result : f64, abserr : f64, neval : i32) -> Self {
        Self{
            result_state : ResultState::Success,
            integration_result : QuadIntegrationResult::new(result , abserr , neval),
        }
    }
    pub fn new_error(result_state : ResultState) -> Self{
        Self{
            result_state,
            integration_result : QuadIntegrationResult::new_error(),
        }
    }

    pub fn new_qng(integrator_result : QngIntegratorResult) -> Self{
        Self{
            result_state : integrator_result.result_state,
            integration_result : QuadIntegrationResult::new_qng(integrator_result.integration_result),
        }
    }
    pub fn new_qag(integrator_result : QagIntegratorResult) -> Self{
        Self{
            result_state : integrator_result.result_state,
            integration_result : QuadIntegrationResult::new_qag(integrator_result.integration_result),
        }
    }

    pub fn new_qag_vec(integrator_result : Qag1DVecIntegratorResult) -> Self{
        Self{
            result_state : integrator_result.result_state,
            integration_result : QuadIntegrationResult::new_qag_vec(integrator_result.integration_result),
        }
    }

    pub fn new_qags(integrator_result : QagIntegratorResult) -> Self{
        Self{
            result_state : integrator_result.result_state,
            integration_result : QuadIntegrationResult::new_qags(integrator_result.integration_result),
        }
    }


    pub fn unwrap(&self) -> QuadIntegrationResult {
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
