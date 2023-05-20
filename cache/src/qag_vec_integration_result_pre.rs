use crate::qage_vec_pre::*;

#[derive(Debug,Clone)]
pub struct QagVecIntegrationResultPre {
    pub result : Vec<f64>,
    pub abserr : Vec<f64>,
    pub neval : i32,
    pub list : Vec<ResultVecPre>,
    pub last : usize,
}

impl QagVecIntegrationResultPre {
    pub fn new(result : Vec<f64>, abserr : Vec<f64>, neval : i32, list : Vec<ResultVecPre>,
               last : usize) -> Self {
        Self{
            result, abserr, neval, list, last
        }
    }

    pub fn new_error() -> Self {
        let zerof = 0.0 ;
        let zeroi = 0;
        let zerou : usize = 0 ;
        let zerovecf = vec![zerof];
        Self{
            result : zerovecf.clone(), abserr : zerovecf.clone(), neval : zeroi, list :
            vec![ResultVecPre::new(0.0, 0.0, zerovecf.clone(), zerovecf)].clone(), last : zerou,
        }
    }
}




