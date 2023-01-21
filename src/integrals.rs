use crate::functions::*;

// all the integrals are computed on the range [0,1]

pub trait Integrator {
    fn integrate(&self,funct : &impl Function, number_of_point : i32) -> f32 ;
    // fn error(&self, funct: &impl Function, number_of_point : i32) -> f32 ;
}

pub struct Rectangular {}
pub struct Trapezoiadal {}
pub struct Simpson1 {}
pub struct Simpson2 {}

impl Integrator for Rectangular {
    fn integrate(&self,funct: &impl Function, number_of_point: i32) -> f32 {
        let mut tot : f32 = funct.evaluate(&0.0) ;
        let m = number_of_point as f32 ;
        for i in 1..number_of_point {
            let j =  i as f32 ;
            tot += funct.evaluate( &(j/m) ) ;
        }
        tot / m
    }
}

impl Integrator for Trapezoiadal {
    fn integrate(&self,funct: &impl Function, number_of_point: i32) -> f32 {
        let mut tot : f32 = 0.5 * (funct.evaluate(&0.0) + funct.evaluate(&1.0)) ;
        let m = number_of_point as f32 ;
        for i in 1..number_of_point {
            let j =  i as f32 ;
            tot += funct.evaluate( &(j/m) ) ;
        }
        tot / m
    }
}

impl Integrator for Simpson1 {
    fn integrate(&self,funct: &impl Function, number_of_point: i32) -> f32 {
        let mut tot : f32 =  funct.evaluate(&0.0) - funct.evaluate(&1.0) ;
        let m = number_of_point as f32 ;
        for i in 1..(number_of_point/2 + 1 ) {
            let j =  i as f32 ;
            tot += 4.0 * funct.evaluate( &((2.0 * j - 1.0 )/m) ) + 2.0 * funct.evaluate(&((2.0 * j)/m) ) ;
        }
        tot / ( 3.0 * m )
    }
}

impl Integrator for Simpson2 {
    fn integrate(&self,funct: &impl Function, number_of_point: i32) -> f32 {
        let mut tot : f32 = funct.evaluate(&0.0) + funct.evaluate(&1.0) ;
        let m = number_of_point as f32 ;
        for i in 1..(number_of_point ) {
            let j =  i as f32 ;
            if i % 3 == 0 {
                tot += 2.0 * funct.evaluate(&(j/m)) ;
            }
            else {
                tot += 3.0 * funct.evaluate(&(j/m)) ;
            }
        }
        3.0 * tot / ( 8.0 * m )
    }
}

