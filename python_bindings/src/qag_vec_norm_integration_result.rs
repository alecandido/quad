
#[derive(Debug,Clone)]
pub struct QagVecNormIntegrationResult<const n :usize> {
    pub result : [f64;n],
    pub abserr : f64,
    pub neval : i32,
    pub last : usize,
}

impl<const n :usize> QagVecNormIntegrationResult<n>{
    pub fn new(result : [f64;n], abserr : f64, neval : i32, last : usize) -> Self {
        Self{
            result, abserr, neval, last
        }
    }

    pub fn new_error() -> Self {
        Self{
            result : [0.0;n], abserr : 0.0, neval : 0, last : 0,
        }
    }
}

