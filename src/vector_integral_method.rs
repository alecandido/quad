use crate::*;
use std::thread;

// all the integrals are computed on the range [0,1]
// !? What is the optimal number_of_points for the abs_max function !?

pub trait VecIntegralMethod {
    fn integrate_non_uniform(&self, funct : & FunctVector, number_of_points : Vec<i32>) -> Vec<f32>;
    fn integrate_non_uniform_parallel(&self, funct: & FunVector , number_of_points_vec: Vec<i32>) -> Vec<f32>;
    fn errors_non_uniform(&self, funct: & FunctVector, number_of_points : Vec<i32>) -> Vec<f32> ;
    fn errors_non_uniform_parallel(&self, funct: & FunVector, number_of_points : Vec<i32>) -> Vec<f32> ;
    fn integrate_uniform(&self, funct : & FunctVector, number_of_points : i32) -> Vec<f32> {
        let number_of_points_vec = vec![number_of_points;funct.components.len()];
        self.integrate_non_uniform(funct,number_of_points_vec)
    }
    fn integrate_uniform_parallel(&self, funct: & FunVector , number_of_points: i32) -> Vec<f32>{
        let number_of_points_vec = vec![number_of_points;funct.components.len()];
        self.integrate_non_uniform_parallel(funct,number_of_points_vec)
    }
    fn errors_uniform(&self, funct: & FunctVector, number_of_points : i32) -> Vec<f32> {
        let number_of_points_vec = vec![number_of_points;funct.components.len()];
        self.errors_non_uniform(funct,number_of_points_vec)
    }
    fn errors_uniform_parallel(&self, funct: & FunVector, number_of_points : i32) -> Vec<f32> {
        let number_of_points_vec = vec![number_of_points;funct.components.len()];
        self.errors_non_uniform_parallel(funct,number_of_points_vec)
    }
    fn relative_error_uniform(&self, funct : & FunctVector, number_of_points : i32) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut relative_errors :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            relative_errors[k] = (self.errors_uniform(funct, number_of_points)[k] / self.integrate_uniform(funct, number_of_points)[k]).abs();
        }
        relative_errors
    }
    fn relative_error_non_uniform(&self, funct : & FunctVector, number_of_points : Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut relative_errors :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            relative_errors[k] = (self.errors_non_uniform(funct, number_of_points.to_vec())[k] / self.integrate_non_uniform(funct, number_of_points.to_vec())[k]).abs();
        }
        relative_errors
    }
    fn relative_error_uniform_parallel(&self, funct : & FunVector, number_of_points : i32) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut relative_errors: Vec<f32> = vec![0.0; dimension];
        for k in 0..dimension {
            relative_errors[k] = (self.errors_uniform_parallel(funct, number_of_points)[k] / self.integrate_uniform_parallel(funct, number_of_points)[k]).abs();
        }
        relative_errors
    }
    fn relative_error_non_uniform_parallel(&self, funct : & FunVector, number_of_points : Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut relative_errors :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            relative_errors[k] = (self.errors_non_uniform_parallel(funct, number_of_points.to_vec())[k] / self.integrate_non_uniform_parallel(funct, number_of_points.to_vec())[k]).abs();
        }
        relative_errors
    }
    fn even_interval(&self) -> bool ; // necessary for Integrator::FixedPrecision.integrate_uniform

}

pub struct VecRectangular {}
pub struct VecTrapezoiadal {}
pub struct VecSimpson1 {}
pub struct VecSimpson2 {}

impl VecIntegralMethod for VecRectangular {
    fn integrate_non_uniform(&self, funct: & FunctVector, number_of_points: Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut results :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            let mut tot : f32 = funct.components[k].evaluate(&0.0);
            let m = number_of_points[k] as f32 ;
            for i in 1..number_of_points[k] {
                let j =  i as f32 ;
                tot += funct.components[k].evaluate( &(j/m) ) ;
            }
            results[k] = tot / m;
        }
        results
    }

    fn integrate_non_uniform_parallel(&self, funct: & FunVector , number_of_points_vec: Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut results: Vec<f32> = vec![];
        let mut handles = vec![];
        for k in 0..dimension {
            let pointer = funct.components[k].clone();
            let number_of_points= number_of_points_vec[k];
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let mut tot: f32 = function.evaluate(&0.0);
                let m = number_of_points as f32;
                for i in 1..number_of_points {
                    let j = i as f32;
                    tot += function.evaluate(&(j / m));
                }
                tot / m
            });
            handles.push(handle);
        }
        for handle in handles {
            results.push(handle.join().unwrap());
        }
        results
    }

    fn errors_non_uniform(&self, funct: & FunctVector, number_of_points : Vec<i32>) -> Vec<f32>{
        let dimension = funct.components.len();
        let mut errors :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            let number_of_points_float = number_of_points[k] as f32;
            let max: f32 = funct.components[k].first_derivative_abs_max();
            errors[k] = max / (2.0 * number_of_points_float);
        }
        errors
    }

    fn errors_non_uniform_parallel(&self, funct: & FunVector, number_of_points : Vec<i32>) -> Vec<f32>{
        let dimension = funct.components.len();
        let mut errors :Vec<f32> = vec![];
        let mut handles = vec![];
        for k in 0..dimension {
            let pointer = funct.components[k].clone();
            let number_of_points_float = number_of_points[k] as f32;
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let max: f32 = function.first_derivative_abs_max();
                max / (2.0 * number_of_points_float)
            });
            handles.push(handle);
        }
        for handle in handles {
            errors.push(handle.join().unwrap());
        }
        errors
    }

    fn even_interval(&self) -> bool {
        true
    }

}

impl VecIntegralMethod for VecTrapezoiadal {
    fn integrate_non_uniform(&self, funct: & FunctVector, number_of_points: Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut results :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            let mut tot: f32 = 0.5 * (funct.components[k].evaluate(&0.0) + funct.components[k].evaluate(&1.0));
            let m = number_of_points[k] as f32;
            for i in 1..number_of_points[k] {
                let j = i as f32;
                tot += funct.components[k].evaluate(&(j / m));
            }
            results[k] = tot / m ;
        }
        results
    }

    fn integrate_non_uniform_parallel(&self, funct: & FunVector , number_of_points_vec: Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut results: Vec<f32> = vec![];
        let mut handles = vec![];
        for k in 0..dimension {
            let pointer = funct.components[k].clone();
            let number_of_points= number_of_points_vec[k];
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let mut tot: f32 = 0.5 * ( function.evaluate(&0.0) + function.evaluate(&1.0) );
                let m = number_of_points as f32;
                for i in 1..number_of_points {
                    let j = i as f32;
                    tot += function.evaluate(&(j / m));
                }
                tot / m
            });
            handles.push(handle);
        }
        for handle in handles {
            results.push(handle.join().unwrap());
        }
        results
    }

    fn errors_non_uniform(&self, funct: & FunctVector, number_of_points : Vec<i32>) -> Vec<f32>{
        let dimension = funct.components.len();
        let mut errors :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            let number_of_points_float = number_of_points[k] as f32;
            let max: f32 = funct.components[k].second_derivative_abs_max();
            errors[k] = max / (12.0 * number_of_points_float.powi(2));
        }
        errors
    }

    fn errors_non_uniform_parallel(&self, funct: & FunVector, number_of_points : Vec<i32>) -> Vec<f32>{
        let dimension = funct.components.len();
        let mut errors :Vec<f32> = vec![];
        let mut handles = vec![];
        for k in 0..dimension {
            let pointer = funct.components[k].clone();
            let number_of_points_float = number_of_points[k] as f32;
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let max: f32 = function.second_derivative_abs_max();
                max / (12.0 * number_of_points_float.powi(2))
            });
            handles.push(handle);
        }
        for handle in handles {
            errors.push(handle.join().unwrap());
        }
        errors
    }

    fn even_interval(&self) -> bool {
        true
    }

}



// number_of_point == number of sub intervals

// Simpson1 ( 1/3 ) : number_of_point has to be even
impl VecIntegralMethod for VecSimpson1 {
    fn integrate_non_uniform(&self, funct: & FunctVector, number_of_points: Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut results :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            let mut tot: f32 = funct.components[k].evaluate(&0.0) - funct.components[k].evaluate(&1.0);
            let m = number_of_points[k] as f32;
            for i in 1..(number_of_points[k] / 2 + 1) {
                let j = i as f32;
                tot += 4.0 * funct.components[k].evaluate(&((2.0 * j - 1.0) / m)) + 2.0 * funct.components[k].evaluate(&((2.0 * j) / m));
            }
            results[k] = tot / (3.0 * m);
        }
        results
    }

    fn integrate_non_uniform_parallel(&self, funct: & FunVector , number_of_points_vec: Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut results: Vec<f32> = vec![];
        let mut handles = vec![];
        for k in 0..dimension {
            let pointer = funct.components[k].clone();
            let number_of_points = number_of_points_vec[k] ;
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let mut tot: f32 = function.evaluate(&0.0) - function.evaluate(&1.0);
                let m = number_of_points as f32;
                for i in 1..(number_of_points/2 + 1) {
                    let j = i as f32;
                    tot += 4.0 * function.evaluate(&((2.0 * j - 1.0) / m)) + 2.0 * function.evaluate(&((2.0 * j) / m));
                }
                tot / (3.0 * m)
            });
            handles.push(handle);
        }
        for handle in handles {
            results.push(handle.join().unwrap());
        }
        results
    }

    fn errors_non_uniform(&self, funct: & FunctVector, number_of_points : Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut errors :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            let number_of_points_float = number_of_points[k] as f32;
            let max: f32 = funct.components[k].fourth_derivative_abs_max();
            errors[k] = max / (180.0 * number_of_points_float.powi(4));
        }
        errors
    }

    fn errors_non_uniform_parallel(&self, funct: & FunVector, number_of_points : Vec<i32>) -> Vec<f32>{
        let dimension = funct.components.len();
        let mut errors :Vec<f32> = vec![];
        let mut handles = vec![];
        for k in 0..dimension {
            let pointer = funct.components[k].clone();
            let number_of_points_float = number_of_points[k] as f32;
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let max: f32 = function.fourth_derivative_abs_max();
                max / (180.0 * number_of_points_float.powi(4))
            });
            handles.push(handle);
        }
        for handle in handles {
            errors.push(handle.join().unwrap());
        }
        errors
    }

    fn even_interval(&self) -> bool {
        true
    }

}

// Simpson2 ( 3/8 ) : number_of_points has to be multiple of 3 !!
impl VecIntegralMethod for VecSimpson2 {
       fn integrate_non_uniform(&self, funct: & FunctVector, number_of_points: Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut results :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            let mut tot: f32 = funct.components[k].evaluate(&0.0) + funct.components[k].evaluate(&1.0);
            let m = number_of_points[k] as f32;
            for i in 1..(number_of_points[k]) {
                let j = i as f32;
                if i % 3 == 0 {
                    tot += 2.0 * funct.components[k].evaluate(&(j / m));
                } else {
                    tot += 3.0 * funct.components[k].evaluate(&(j / m));
                }
            }
            results[k] = 3.0 * tot / (8.0 * m);
        }
        results
    }

    fn integrate_non_uniform_parallel(&self, funct: & FunVector , number_of_points_vec: Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut results: Vec<f32> = vec![];
        let mut handles = vec![];
        for k in 0..dimension {
            let pointer = funct.components[k].clone();
            let number_of_points= number_of_points_vec[k];
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let mut tot: f32 = function.evaluate(&0.0) + function.evaluate(&1.0);
                let m = number_of_points as f32;
                for i in 1..(number_of_points) {
                    let j = i as f32;
                    if i % 3 == 0 {
                        tot += 2.0 * function.evaluate(&(j / m));
                    } else {
                        tot += 3.0 * function.evaluate(&(j / m));
                    }
                }
                3.0 * tot / (8.0 * m)
            });
            handles.push(handle);
        }
        for handle in handles {
            results.push(handle.join().unwrap());
        }
        results
    }

    fn errors_non_uniform(&self, funct: & FunctVector, number_of_points : Vec<i32>) -> Vec<f32> {
        let dimension = funct.components.len();
        let mut errors :Vec<f32> = vec![0.0;dimension];
        for k in 0..dimension {
            let number_of_points_float = number_of_points[k] as f32;
            let max: f32 = funct.components[k].fourth_derivative_abs_max();
            errors[k] =  max / (80.0 * number_of_points_float.powi(4));
        }
        errors
    }

    fn errors_non_uniform_parallel(&self, funct: & FunVector, number_of_points : Vec<i32>) -> Vec<f32>{
        let dimension = funct.components.len();
        let mut errors :Vec<f32> = vec![];
        let mut handles = vec![];
        for k in 0..dimension {
            let pointer = funct.components[k].clone();
            let number_of_points_float = number_of_points[k] as f32;
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let max: f32 = function.fourth_derivative_abs_max();
                max / (80.0 * number_of_points_float.powi(4))
            });
            handles.push(handle);
        }
        for handle in handles {
            errors.push(handle.join().unwrap());
        }
        errors
    }

    fn even_interval(&self) -> bool {
        false
    }

}

