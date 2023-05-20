
#[derive(Debug,Clone)]
pub struct QagVecIntegrationResult {
    pub result : [f64;4],
    pub abserr : [f64;4],
    pub neval : i32,
    pub alist : Vec<f64>,
    pub blist : Vec<f64>,
    pub rlist : Vec<[f64;4]>,
    pub elist : Vec<[f64;4]>,
    pub last : usize,
}

impl QagVecIntegrationResult {
    pub fn new(result : [f64;4], abserr : [f64;4], neval : i32, alist : Vec<f64>,
               blist : Vec<f64>, rlist : Vec<[f64;4]>, elist : Vec<[f64;4]>,
               last : usize) -> Self {
        Self{
            result, abserr, neval, alist, blist, rlist, elist, last
        }
    }

    pub fn new_error() -> Self {
        Self{
            result : [0.0;4], abserr : [0.0;4], neval : 0, alist : vec![0.0],
            blist : vec![0.0], rlist : vec![[0.0;4]], elist : vec![[0.0;4]], last : 0,
        }
    }
}

