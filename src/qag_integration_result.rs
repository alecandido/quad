#[derive(Debug,Clone)]
pub struct QagIntegrationResult {
    result : f64,
    abserr : f64,
    neval : i32,
    alist : Vec<f64>,
    blist : Vec<f64>,
    rlist : Vec<f64>,
    elist : Vec<f64>,
    iord : Vec<usize>,
    last : usize,
}
#[derive(Clone,Debug)]
pub enum ResultState{
    Success,
    Failure,
    Invalid,
    MaxIteration,
    BadTolerance,
    BadFunction,

}

impl QagIntegrationResult {
    pub fn new(result : f64, abserr : f64, neval : i32, alist : Vec<f64>,
               blist : Vec<f64>, rlist : Vec<f64>, elist : Vec<f64>, iord : Vec<usize>,
               last : usize) -> Self {
        Self{
            result, abserr, neval, alist, blist, rlist, elist, iord, last
        }
    }

    pub fn new_error() -> Self {
        let zerof = 0.0 ;
        let zerovecf = vec![zerof];
        let zeroi = 0;
        let zerou : usize = 0 ;
        let zerovecu = vec![zerou];
        Self{
            result : zerof, abserr : zerof, neval : zeroi, alist : zerovecf.clone(),
            blist : zerovecf.clone(), rlist : zerovecf.clone(), elist : zerovecf.clone(),
            iord : zerovecu, last : zerou,
        }
    }
}




