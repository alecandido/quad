use std::simd::{f64x8, Simd};
use crate::qage_vec_pre::*;

#[derive(Clone,Debug)]
pub struct ResultVec8 {
    pub a : f64,
    pub b : f64,
    pub result : Simd<f64,8>,
    pub error : Simd<f64,8>,
}

impl ResultVec8 {
    pub fn new(a : f64, b : f64, result : Simd<f64,8>, error : Simd<f64,8>) -> Self{
        Self{
            a,b,result,error
        }
    }
}

#[derive(Debug,Clone)]
pub struct QagVec8IntegrationResult {
    pub result : Simd<f64,8>,
    pub abserr : Simd<f64,8>,
    pub neval : i32,
    pub list : Vec<ResultVec8>,
    pub last : usize,
}

impl QagVec8IntegrationResult {
    pub fn new(result : Simd<f64,8>, abserr : Simd<f64,8>, neval : i32, list : Vec<ResultVec8>,
               last : usize) -> Self {
        Self{
            result, abserr, neval, list, last
        }
    }

    pub fn new_error() -> Self {
        Self{
            result : f64x8::splat(0.0), abserr : f64x8::splat(0.0), neval : 0, list :
            vec![ResultVec8::new(0.0, 0.0, f64x8::splat(0.0), f64x8::splat(0.0))].clone(), last : 0,
        }
    }
}




