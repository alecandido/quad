
#[derive(Debug,Clone)]
pub struct QagParIntegrationResult<const N:usize> {
    pub result : [f64; N],
    pub abserr : f64,
    pub neval : i32,
    pub last : usize,
}

impl<const N:usize> QagParIntegrationResult<N>{
    pub fn new(result : [f64; N], abserr : f64, neval : i32, last : usize) -> Self {
        Self{
            result, abserr, neval, last
        }
    }

    pub fn new_error() -> Self {
        Self{
            result : [0.0; N], abserr : 0.0, neval : 0, last : 0,
        }
    }
}

