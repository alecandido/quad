#[derive(Debug,Clone)]
pub struct QngIntegrationResult {
    pub result : f64,
    pub abserr : f64,
    pub neval : i32,
}

impl QngIntegrationResult {
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
}
