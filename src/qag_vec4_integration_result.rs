use std::simd::{f64x4, Simd};
use crate::qage_vec::*;

#[derive(Clone,Debug)]
pub struct ResultVec4 {
    pub a : f64,
    pub b : f64,
    pub result : Simd<f64,4>,
    pub error : Simd<f64,4>,
}

impl ResultVec4{
    pub fn new(a : f64, b : f64, result : Simd<f64,4>, error : Simd<f64,4>) -> Self{
        Self{
            a,b,result,error
        }
    }
}

#[derive(Debug,Clone)]
pub struct QagVec4IntegrationResult {
    pub result : Simd<f64,4>,
    pub abserr : Simd<f64,4>,
    pub neval : i32,
    pub list : Vec<ResultVec4>,
    pub last : usize,
}

impl QagVec4IntegrationResult {
    pub fn new(result : Simd<f64,4>, abserr : Simd<f64,4>, neval : i32, list : Vec<ResultVec4>,
               last : usize) -> Self {
        Self{
            result, abserr, neval, list, last
        }
    }

    pub fn new_error() -> Self {
        Self{
            result : f64x4::splat(0.0), abserr : f64x4::splat(0.0), neval : 0, list :
            vec![ResultVec4::new(0.0,0.0, f64x4::splat(0.0),f64x4::splat(0.0))].clone(), last : 0,
        }
    }
}




