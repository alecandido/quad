extern crate nalgebra as na;

pub mod functions;
pub mod integral_method;
pub mod integrator;
pub mod funct_vector;
pub mod vector_integral_method;
pub mod vector_intergrator;
pub mod gauss_legendre;
pub mod qng;
pub mod quad_integrator;

use functions::*;
use integral_method::*;
use integrator::*;
use funct_vector::*;
use vector_integral_method::*;
use vector_intergrator::*;
use std::time::Instant;
use std::sync::{Arc, Mutex};
use qng::*;
use quad_integrator::*;


pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use crate::gauss_legendre::{Gauss, GaussLegendre};
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }

    #[test]
    fn qng_test(){
        let (a,b) = (0.0,1.0);
        let f = |x:f64|  (x*x.sin()).powf(1.5);
        let integral_method = Qng{};
        let res = integral_method.integrate(&f, a, b, 2.2436292642488107e-15, 0.0);
        println!("{:?}",res);
        let f1 = |x:f64|  (x*x.sin()).powf(1.5);
        let f2 = |x:f64|  (x.powi(3)*x.cos()).powf(2.5);
        let f3 = |x:f64|  (x.powi(5)*x.tan()).powf(2.5);
        let f4 = |x:f64|  (x.powi(7)*x.atan()).powf(2.5);
        let fun = FnVec  {
            components : vec![Box::new(f1),Box::new(f2),Box::new(f3),Box::new(f4)],
        };
        let fun2 = FnVecP{
            components : vec![Arc::new(Mutex::new(f1)),Arc::new(Mutex::new(f2)),
                            Arc::new(Mutex::new(f3)),Arc::new(Mutex::new(f4))],
        };

        let int = QngIntegrator{
            a : vec![0.0,0.0,0.0,0.0],
            b : vec![1.0,1.0,1.0,1.0],
            epsabs : vec![0.0,0.0,0.0,0.0],
            epsrel : vec![1.0,1.0,1.0,1.0],
        };
        let qng = Qng{};
        let qngg1 = Qng{};
        let qngg2 = Qng{};
        let qngg3 = Qng{};
        let qngg4 = Qng{};
        let qng2 = IntVec{
            components : vec![Arc::new(Mutex::new(qngg1)),Arc::new(Mutex::new(qngg2)),
                              Arc::new(Mutex::new(qngg3)),Arc::new(Mutex::new(qngg4))],
        };

        //let time1 = Instant::now();
        let res1 =int.integrate(&qng,&fun);
        //let time2 = Instant::now();
        let res2 = int.integrate_p(&qng2,&fun2);
        //let time3 = Instant::now();

        //println!("{:?}\n{:?}\n{:?},{:?}", res1,res2,time2-time1,time3-time2 );


    }

    #[test]
    fn gauss_legendre(){
        let precision : f64 = 0.001 ; // set the relative errors required
        let parameters = random_vector();
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        //let polynomial = Sin{i : 0} ;
        let simpson2 = Simpson2{};
        let integral_precision = FixedPrecision{precision};
        let exact_integral = exact_polynom_integral(&parameters);
        //let exact_integral = - 1.0_f64.cos() + 1.0 ;

        let gauss = GaussLegendre::new( 8);
        println!("{}\n{:?}\n{:?}\n{:?}", parameters.len(), integral_precision.integrate(&simpson2,&polynomial),
                exact_integral, gauss.integrate(&polynomial));

    }



    //-------------------   RELATIVE ERRORS   -----------------

    #[test]
    fn parallel() {
        let test : bool = false; // set to true if you want to have the results printed out
        let integration_method = VecSimpson2{};
        let number_of_points = vec![10000;9];
        let number_of_points2 = vec![10000;9];

        let parameters1 = [random_vector(),random_vector(),random_vector()];
        let parameters2 = [random_vector(),random_vector(),random_vector()];
        let parameters3 = [random_vector(),random_vector(),random_vector()];
        let (polynomial1, polynomial2,polynomial3) =
            (PolynomialFunction::new(parameters1[0].to_vec()),
             PolynomialFunction::new(parameters1[1].to_vec()),
             PolynomialFunction::new(parameters1[2].to_vec()) );
        let (polynomial4, polynomial5,polynomial6) =
            (PolynomialFunction::new(parameters2[0].to_vec()),
             PolynomialFunction::new(parameters2[1].to_vec()),
             PolynomialFunction::new(parameters2[2].to_vec()) );
        let (polynomial7, polynomial8,polynomial9) =
            (PolynomialFunction::new(parameters3[0].to_vec()),
             PolynomialFunction::new(parameters3[1].to_vec()),
             PolynomialFunction::new(parameters3[2].to_vec()) );
        let collection : FunVector = FunVector{ components : vec![Arc::new(Mutex::new(polynomial1.clone())),
                                                                  Arc::new(Mutex::new(polynomial2.clone())),
                                                                  Arc::new(Mutex::new(polynomial3.clone())),
                                                                  Arc::new(Mutex::new(polynomial4.clone())),
                                                                  Arc::new(Mutex::new(polynomial5.clone())),
                                                                  Arc::new(Mutex::new(polynomial6.clone())),
                                                                  Arc::new(Mutex::new(polynomial7.clone())),
                                                                  Arc::new(Mutex::new(polynomial8.clone())),
                                                                  Arc::new(Mutex::new(polynomial9.clone()))]};
        let collection2 : FunctVector = FunctVector{ components : vec![Box::new(polynomial1),Box::new(polynomial2),Box::new(polynomial3),
                                                                       Box::new(polynomial4),Box::new(polynomial5),Box::new(polynomial6),
                                                                       Box::new(polynomial7),Box::new(polynomial8),Box::new(polynomial9)]};

        let time3 = Instant::now();
        let res_parall = integration_method.integrate_non_uniform_parallel(&collection, number_of_points);
        let time4 = Instant::now();
        let res = integration_method.integrate_non_uniform(&collection2, number_of_points2);
        let time5 = Instant::now();
        let parallel_time = time4 - time3;
        let non_parallel_time = time5 - time4;
        let proportion = parallel_time.as_secs_f64()/non_parallel_time.as_secs_f64();
        assert!(res == res_parall && parallel_time < non_parallel_time );

        if test {
            println!("{:?},{:?}", res, res_parall);
            println!("Parallel function duration :{:?},\nNon parallel function duration :{:?}\n{:?}", parallel_time, non_parallel_time,proportion);
        }

    }


    #[test]
    fn error_comparison() {
        let test : bool = false; // set to true if you want to have the results printed out
        let precision : f64 = 0.5 ; // set the relative errors required

        let parameters = random_vector();
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        //let polynomial = Sin{i : 0} ;
        let integral_precision = FixedPrecision{precision};
        let rectangular = Rectangular{};
        let trapezoidal = Trapezoiadal{};
        let simpson1 = Simpson1{};
        let simpson2 = Simpson2{};
        let exact_integral = exact_polynom_integral(&parameters) ;
        //let exact_integral = - 1.0_f64.cos() + 1.0 ;


        let rel_error_real_rectangular = ((integral_precision.integrate(&rectangular,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_trapezoidal = ((integral_precision.integrate(&trapezoidal,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_simpson1 = ((integral_precision.integrate(&simpson1,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_simpson2 = ((integral_precision.integrate(&simpson2,&polynomial) - exact_integral).abs() / exact_integral).abs() ;

        if test {
            println!("Exact integral : {}", exact_integral);
            println!("The random polynomial function has order : {}", parameters.len());
            println!("Rectangular: {} , with {} points ( relative errors bound : {}; real errors : {})",
                     integral_precision.integrate(&rectangular, &polynomial),
                     integral_precision.number_of_points(&rectangular, &polynomial),
                     integral_precision.relative_error(&rectangular, &polynomial),
                     rel_error_real_rectangular);
            println!("Trapezoidal: {}, with {} points ( relative errors bound : {}; real errors : {})",
                     integral_precision.integrate(&trapezoidal, &polynomial),
                     integral_precision.number_of_points(&trapezoidal, &polynomial),
                     integral_precision.relative_error(&trapezoidal, &polynomial),
                     rel_error_real_trapezoidal);
            println!("Simpson1: {}, with {} points ( relative errors bound : {}; real errors : {})",
                     integral_precision.integrate(&simpson1, &polynomial),
                     integral_precision.number_of_points(&simpson1, &polynomial),
                     integral_precision.relative_error(&simpson1, &polynomial),
                     rel_error_real_simpson1);
            println!("Simpson2: {}, with {} points ( relative errors bound : {}; real errors : {})",
                     integral_precision.integrate(&simpson2, &polynomial),
                     integral_precision.number_of_points(&simpson2, &polynomial),
                     integral_precision.relative_error(&simpson2, &polynomial),
                     rel_error_real_simpson2);
        }

        assert!(rel_error_real_rectangular < precision && rel_error_real_trapezoidal < precision
                && rel_error_real_simpson1 < precision && rel_error_real_simpson2 < precision)

    }

//----------------------- TEST OF FUNCTVECTORS


    #[test]
    fn funct_vector() {
        let test : bool = false; // set to true if you want to have the results printed out
        let precision :f64 = 0.5; // set the relative errors required

        let parameters = [random_vector(),random_vector(),random_vector()];
        let functions_order = [ parameters[0].len(),parameters[1].len(),parameters[2].len()];
        let (polynomial1, polynomial2,polynomial3) =
            (PolynomialFunction::new(parameters[0].to_vec()),
            PolynomialFunction::new(parameters[1].to_vec()),
            PolynomialFunction::new(parameters[2].to_vec()) );
        let collection : FunctVector = FunctVector{ components : vec![Box::new(polynomial1),Box::new(polynomial2),Box::new(polynomial3)]};
        let rectangular = VecRectangular{};
        let trapezoidal = VecTrapezoiadal{};
        let simpson1 = VecSimpson1{};
        let simpson2 = VecSimpson2{};
        let vec_integral_method = VecFixedPrecision{precision};
        let exact_integral = [exact_polynom_integral(&parameters[0]),
            exact_polynom_integral(&parameters[1]),exact_polynom_integral(&parameters[2])];

        // UNIFORM
        let result = [vec_integral_method.integrate_uniform(&rectangular,&collection),
            vec_integral_method.integrate_uniform(&trapezoidal,&collection),
            vec_integral_method.integrate_uniform(&simpson1,&collection),
            vec_integral_method.integrate_uniform(&simpson2,&collection)];
        let number_of_points = [vec_integral_method.number_of_points_uniform(&rectangular,&collection),
            vec_integral_method.number_of_points_uniform(&trapezoidal,&collection),
            vec_integral_method.number_of_points_uniform(&simpson1,&collection),
            vec_integral_method.number_of_points_uniform(&simpson2,&collection)];
        let relative_errors = [vec_integral_method.relative_errors_uniform(&rectangular, &collection),
            vec_integral_method.relative_errors_uniform(&trapezoidal, &collection),
            vec_integral_method.relative_errors_uniform(&simpson1, &collection),
            vec_integral_method.relative_errors_uniform(&simpson2, &collection)];
        let rel_error_real = [[((result[0][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
            ((result[0][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
            ((result[0][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result[1][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result[1][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result[1][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result[2][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result[2][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result[2][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result[3][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result[3][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result[3][2] - exact_integral[2]).abs() / exact_integral[2]).abs()]];


        // NON UNIFORM (nu)
        let result_nu = [vec_integral_method.integrate_non_uniform(&rectangular,&collection),
            vec_integral_method.integrate_non_uniform(&trapezoidal,&collection),
            vec_integral_method.integrate_non_uniform(&simpson1,&collection),
            vec_integral_method.integrate_non_uniform(&simpson2,&collection)];
        let number_of_points_nu = [vec_integral_method.number_of_points_non_uniform(&rectangular,&collection),
            vec_integral_method.number_of_points_non_uniform(&trapezoidal,&collection),
            vec_integral_method.number_of_points_non_uniform(&simpson1,&collection),
            vec_integral_method.number_of_points_non_uniform(&simpson2,&collection)];
        let relative_errors_nu = [vec_integral_method.relative_errors_non_uniform(&rectangular, &collection),
            vec_integral_method.relative_errors_non_uniform(&trapezoidal, &collection),
            vec_integral_method.relative_errors_non_uniform(&simpson1, &collection),
            vec_integral_method.relative_errors_non_uniform(&simpson2, &collection)];
        let rel_error_real_nu = [[((result_nu[0][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
            ((result_nu[0][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
            ((result_nu[0][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result_nu[1][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result_nu[1][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result_nu[1][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result_nu[2][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result_nu[2][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result_nu[2][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result_nu[3][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result_nu[3][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result_nu[3][2] - exact_integral[2]).abs() / exact_integral[2]).abs()]];

        if test {
            println!("UNIFORM:\nThe random polynomial functions have order : {:?}", functions_order);
            println!("Notation : [Rectangular, Trapezoidal, Simpson1, Simpson2]\nNumber of point used : {:?} ", number_of_points);
            println!("Integral result: [first function, second function, third function]\n{:?}", result);
            //println!("{:?}", rectangular.integrate_uniform_parallel(&collection,number_of_points[0]));
            println!("Relative errors bound:\n{:?}", relative_errors);
            println!("NON-UNIFORM:\nThe random polynomial functions have order : {:?}", functions_order);
            println!("Notation : [Rectangular, Trapezoidal, Simpson1, Simpson2]\nNumber of point used : {:?} ", number_of_points_nu);
            println!("Integral result: [first function, second function, third function]\n{:?}", result_nu);
            println!("Relative errors bound:\n{:?}", relative_errors_nu);
        }

        assert!(rel_error_real[0][0] < precision && rel_error_real[0][1] < precision &&
            rel_error_real[0][2] < precision &&
            rel_error_real[1][0] < precision && rel_error_real[1][1] < precision &&
            rel_error_real[1][2] < precision &&
            rel_error_real[2][0] < precision && rel_error_real[2][1] < precision &&
            rel_error_real[2][2] < precision &&
            rel_error_real[3][0] < precision && rel_error_real[3][1] < precision &&
            rel_error_real[3][2] < precision &&
            rel_error_real_nu[0][0] < precision && rel_error_real_nu[0][1] < precision &&
            rel_error_real_nu[0][2] < precision &&
            rel_error_real_nu[1][0] < precision && rel_error_real_nu[1][1] < precision &&
            rel_error_real_nu[1][2] < precision &&
            rel_error_real_nu[2][0] < precision && rel_error_real_nu[2][1] < precision &&
            rel_error_real_nu[2][2] < precision &&
            rel_error_real_nu[3][0] < precision && rel_error_real_nu[3][1] < precision &&
            rel_error_real_nu[3][2] < precision)
        }





}




/* test for integrate_parallel ( in integral_method.rs )
        let parameters = random_vector();
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        let rectangular = Rectangular{};
        let integration_method = VecRectangular{};
        let now = Instant::now();


        let res_parall = rectangular.integrate_parallel(&polynomial,1000000);
        let time1 = Instant::now();
        let res = rectangular.integrate(&polynomial,1000000);
        let time2 = Instant::now();
        let exact = exact_polynom_integral(&parameters);

        println!("{res},{res_parall},{exact}");
        println!("{:?},{:?},{:?} , {:?}, {:?}",now,time1,time2, (time1 - now), time2-time1);
         */

/*
fn gauss() {
    let n = 3;
    let f = |x:f32|   x.powi(n) ;
    fn square<F>(f : F, x : f32) -> f32
        where F: Fn(f32)->f32{
        f(x).powi(2)
    }
    println!("{},{}", f(3.0), square(f,3.0));
}

 */