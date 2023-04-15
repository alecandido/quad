use std::sync::{Arc, Mutex};
use crate::qage_1dvec2::*;

#[derive(Clone,Debug)]
pub struct ArcResult {
    pub inner : Arc<Mutex<Result>>,
}
impl ArcResult{
    pub fn new(a : f64, b : f64, result : f64, error : f64) -> Self{
        Self{ inner : Arc::new(Mutex::new(Result::new(a,b,result,error)))}
    }
}
unsafe impl Send for ArcResult{}
unsafe impl Sync for ArcResult{}

#[derive(Debug,Clone)]
pub struct Qag1DVecParIntegrationResult {
    pub result : f64,
    pub abserr : f64,
    pub neval : i32,
    pub list : Vec<ArcResult>,
    pub last : usize,
}

impl Qag1DVecParIntegrationResult {
    pub fn new(result : f64, abserr : f64, neval : i32, list : Vec<ArcResult>,
               last : usize) -> Self {
        Self{
            result, abserr, neval, list, last
        }
    }

    pub fn new_error() -> Self {
        let zerof = 0.0 ;
        let zeroi = 0;
        let zerou : usize = 0 ;
        Self{
            result : zerof, abserr : zerof, neval : zeroi, list :
            vec![ArcResult::new(0.0,0.0, 0.0,0.0)].clone(), last : zerou,
        }
    }
}




