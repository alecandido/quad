use crate::qag_integration_result::QagIntegrationResult;
use crate::qng_integration_result::QngIntegrationResult;

#[derive(Debug,Clone)]
pub struct QuadIntegrationResult {
    pub result : f64,
    pub abserr : f64,
    pub neval : i32,
}

impl QuadIntegrationResult {
    pub fn new(result : f64, abserr : f64, neval : i32) -> Self {
        Self{
            result, abserr, neval,
        }
    }

    pub fn new_error() -> Self {
        let zerof = 0.0 ;
        let zeroi = 0;
        Self{
            result : zerof, abserr : zerof, neval : zeroi,
        }
    }

    pub fn new_qng(integration_result : QngIntegrationResult) -> Self{
        Self{
            result : integration_result.result,
            abserr : integration_result.abserr,
            neval : integration_result.neval
        }
    }

    pub fn new_qag(integration_result : QagIntegrationResult) -> Self{
        Self{
            result : integration_result.result,
            abserr : integration_result.abserr,
            neval : integration_result.neval,
        }
    }

    pub fn new_qags(integration_result : QagIntegrationResult) -> Self{
        Self{
            result : integration_result.result,
            abserr : integration_result.abserr,
            neval : integration_result.neval,
        }
    }
}
