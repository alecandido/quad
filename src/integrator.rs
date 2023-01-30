use crate::*;

pub trait Integrator{
    fn integrate(&self, integral_method : &impl IntegralMethod, funct : &impl Function) -> f32 ;
    fn number_of_points(&self, integral_method : &impl IntegralMethod, funct : &impl Function) -> i32;
    fn relitive_error(&self, integral_method : &impl IntegralMethod, funct : &impl Function) -> f32 {
        integral_method.relative_error(funct,self.number_of_points(integral_method,funct))
    }
}

pub struct FixedPrecision {
    pub precision : f32 ,
}

impl Integrator for FixedPrecision {
    fn integrate(&self, integral_method: &impl IntegralMethod, funct: &impl Function) -> f32 {
        let mut number_of_points: i32 ;

        if  integral_method.even_interval() {
            number_of_points = 2;
            while integral_method.relative_error(funct, number_of_points) > self.precision {
                number_of_points += 2;
            }
        }
        else {
            number_of_points = 3;
            while integral_method.relative_error(funct, number_of_points) > self.precision {
                number_of_points += 3;
            }
        }
        integral_method.integrate(funct,number_of_points)
    }
    fn number_of_points(&self, integral_method : &impl IntegralMethod, funct : &impl Function) -> i32 {
        let mut number_of_points: i32 ;

        if  integral_method.even_interval() {
            number_of_points = 2;
            while integral_method.relative_error(funct, number_of_points) > self.precision {
                number_of_points += 2;
            }
        }
        else {
            number_of_points = 3;
            while integral_method.relative_error(funct, number_of_points) > self.precision {
                number_of_points += 3;
            }
        }
        number_of_points
    }
}