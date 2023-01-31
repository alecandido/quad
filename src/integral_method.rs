use crate::*;

// all the integrals are computed on the range [0,1]
// !? What is the optimal number_of_points for the abs_max function !?

pub trait IntegralMethod {
    fn integrate(&self,funct : &impl Function, number_of_point : i32) -> f32 ;
    fn error(&self, funct: &impl Function, number_of_ponts : i32) -> f32 ;
    fn relative_error(&self, funct : &impl Function , number_of_points : i32) -> f32 {
        self.error(funct,number_of_points) / self.integrate(funct, number_of_points).abs()
    }
    fn even_interval(&self) -> bool ; // necessary for Integrator::FixedPrecision.integrate_uniform
}

pub struct Rectangular {}
pub struct Trapezoiadal {}
pub struct Simpson1 {}
pub struct Simpson2 {}

impl IntegralMethod for Rectangular {
    fn integrate(&self,funct: &impl Function, number_of_point: i32) -> f32 {
        let mut tot : f32 = funct.evaluate(&0.0) ;
        let m = number_of_point as f32 ;
        for i in 1..number_of_point {
            let j =  i as f32 ;
            tot += funct.evaluate( &(j/m) ) ;
        }
        tot / m
    }

    fn error(&self, funct: &impl Function, number_of_points : i32) -> f32 {
        let number_of_points_float = number_of_points as f32;
        let max: f32 = funct.first_derivative_abs_max();
        let error: f32 = max / (2.0 * number_of_points_float);
        error
    }

    fn even_interval(&self) -> bool {
        true
    }
}

impl IntegralMethod for Trapezoiadal {
    fn integrate(&self,funct: &impl Function, number_of_point: i32) -> f32 {
        let mut tot : f32 = 0.5 * (funct.evaluate(&0.0) + funct.evaluate(&1.0)) ;
        let m = number_of_point as f32 ;
        for i in 1..number_of_point {
            let j =  i as f32 ;
            tot += funct.evaluate( &(j/m) ) ;
        }
        tot / m
    }

    fn error(&self, funct: &impl Function, number_of_points : i32) -> f32 {
        let number_of_points_float = number_of_points as f32;
        let max: f32 = funct.second_derivative_abs_max();
        let error: f32 = max / (12.0 * number_of_points_float.powi(2));
        error
    }
    fn even_interval(&self) -> bool {
        true
    }
}



// number_of_point == number of sub intervals

// Simpson1 ( 1/3 ) : number_of_point has to be even
impl IntegralMethod for Simpson1 {
    fn integrate(&self,funct: &impl Function, number_of_point: i32) -> f32 {
        let mut tot : f32 =  funct.evaluate(&0.0) - funct.evaluate(&1.0) ;
        let m = number_of_point as f32 ;
        for i in 1..(number_of_point/2 + 1 ) {
            let j =  i as f32 ;
            tot += 4.0 * funct.evaluate( &((2.0 * j - 1.0 )/m) ) + 2.0 * funct.evaluate(&((2.0 * j)/m) ) ;
        }
        tot / ( 3.0 * m )
    }

    fn error(&self, funct: &impl Function , number_of_points : i32) -> f32 {
        let number_of_points_float = number_of_points as f32;
        let max: f32 = funct.fourth_derivative_abs_max();
        let error: f32 = max / (180.0 * number_of_points_float.powi(4));
        error
    }
    fn even_interval(&self) -> bool {
        true
    }
}

// Simpson2 ( 3/8 ) : number_of_points has to be multiple of 3 !!
impl IntegralMethod for Simpson2 {
    fn integrate(&self,funct: &impl Function, number_of_points: i32) -> f32 {
        let mut tot : f32 = funct.evaluate(&0.0) + funct.evaluate(&1.0) ;
        let m = number_of_points as f32 ;
        for i in 1..(number_of_points) {
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

    fn error(&self, funct: &impl Function, number_of_points : i32) -> f32 {
        let number_of_points_float = number_of_points as f32;
        let max: f32 = funct.fourth_derivative_abs_max();
        let error: f32 =  max / (80.0 * number_of_points_float.powi(4));
        error
    }

    fn even_interval(&self) -> bool {
        false
    }
}

