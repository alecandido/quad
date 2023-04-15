use crate::qage_1dvec2::*;

#[derive(Debug,Clone)]
pub struct Qag1DVecIntegrationResult {
    pub result : f64,
    pub abserr : f64,
    pub neval : i32,
    pub list : Vec<Result>,
    pub last : usize,
}

impl Qag1DVecIntegrationResult {
    pub fn new(result : f64, abserr : f64, neval : i32, list : Vec<Result>,
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
            vec![Result::new(0.0,0.0, 0.0,0.0)].clone(), last : zerou,
        }
    }
}




