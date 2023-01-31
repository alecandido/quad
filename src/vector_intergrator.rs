use crate::*;

pub trait VecIntegrator{
    fn integrate_uniform(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> Vec<f32> ;
    fn integrate_non_uniform(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> Vec<f32> ;
    fn number_of_points_uniform(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> i32;
    fn number_of_points_non_uniform(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> Vec<i32>;
    fn relative_errors(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> Vec<f32> {
        integral_method.relative_error(funct,self.number_of_points_uniform(integral_method,funct))
    }
}

pub struct VecFixedPrecision {
    pub precision : f32 ,
}

impl VecIntegrator for VecFixedPrecision {
    fn integrate_uniform(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> Vec<f32>  {
        let mut number_of_points: i32 ;
        let dimension = funct.components.len();
        let mut condition : Vec<bool> = vec![false;dimension];
        let check : Vec<bool> = vec![true;dimension];

        if  integral_method.even_interval() {
            number_of_points = 0;
            while condition != check {
                number_of_points += 2;
                for i in 0..dimension {
                    if integral_method.relative_error(funct, number_of_points)[i] > self.precision {
                        condition[i] = false;
                    }
                    else {condition[i] = true;}
                }
                }
            }

        else {
            number_of_points = 0;
            while condition != check {
                number_of_points += 3;
                for i in 0..dimension {
                    if integral_method.relative_error(funct, number_of_points)[i] > self.precision {
                        condition[i] = false;
                    }
                    else {condition[i] = true;}
                }
            }
        }
        integral_method.integrate_uniform(funct, number_of_points)
    }

    fn number_of_points_uniform(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> i32  {
        let mut number_of_points: i32 ;
        let dimension = funct.components.len();
        let mut condition : Vec<bool> = vec![false;dimension];
        let check : Vec<bool> = vec![true;dimension];

        if  integral_method.even_interval() {
            number_of_points = 0;
            while condition != check {
                number_of_points += 2;
                for i in 0..dimension {
                    if integral_method.relative_error(funct, number_of_points)[i] > self.precision {
                        condition[i] = false;
                    }
                    else {condition[i] = true;}
                }
            }
        }

        else {
            number_of_points = 0;
            while condition != check {
                number_of_points += 3;
                for i in 0..dimension {
                    if integral_method.relative_error(funct, number_of_points)[i] > self.precision {
                        condition[i] = false;
                    }
                    else {condition[i] = true;}
                }
            }
        }
        number_of_points
    }

    fn integrate_non_uniform(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> Vec<f32>  {
        let mut number_of_points: i32 ;
        let dimension = funct.components.len();
        let mut number_of_points_vec : Vec<i32> = vec![1;dimension];
        let mut condition : Vec<bool> = vec![false;dimension];
        let check : Vec<bool> = vec![true;dimension];

        if  integral_method.even_interval() {
            number_of_points = 0;
            while condition != check {
                number_of_points += 2;
                for i in 0..dimension {
                    if integral_method.relative_error(funct, number_of_points)[i] > self.precision {
                        condition[i] = false;
                    }
                    else {
                        condition[i] = true;
                        number_of_points_vec[i] = number_of_points;
                    }
                }
            }

        }

        else {
            number_of_points = 0;
            while condition != check {
                number_of_points += 3;
                for i in 0..dimension {
                    if integral_method.relative_error(funct, number_of_points)[i] > self.precision {
                        condition[i] = false;
                    }
                    else {
                        condition[i] = true;
                        number_of_points_vec[i] = number_of_points;
                    }
                }
            }
        }
        integral_method.integrate_non_uniform(funct, number_of_points_vec)
    }

    fn number_of_points_non_uniform(&self, integral_method : &impl VecIntegralMethod, funct : & FunctVector) -> Vec<i32>  {
        let mut number_of_points: i32 ;
        let dimension = funct.components.len();
        let mut condition : Vec<bool> = vec![false;dimension];
        let mut number_of_points_vec : Vec<i32> = vec![1;dimension];
        let check : Vec<bool> = vec![true;dimension];

        if  integral_method.even_interval() {
            number_of_points = 0;
            while condition != check {
                number_of_points += 2;
                for i in 0..dimension {
                    if integral_method.relative_error(funct, number_of_points)[i] > self.precision {
                        condition[i] = false;
                    }
                    else {
                        condition[i] = true;
                        number_of_points_vec[i] = number_of_points;
                    }
                }
            }
        }

        else {
            number_of_points = 0;
            while condition != check {
                number_of_points += 3;
                for i in 0..dimension {
                    if integral_method.relative_error(funct, number_of_points)[i] > self.precision {
                        condition[i] = false;
                    }
                    else {
                        condition[i] = true;
                        number_of_points_vec[i] = number_of_points;
                    }
                }
            }
        }
        number_of_points_vec
    }

}