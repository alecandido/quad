pub mod functions;
pub mod integral_method;
pub mod integrator;
pub mod funct_vector;
pub mod vector_integral_method;
pub mod vector_intergrator;

use functions::*;
use integral_method::*;
use integrator::*;
use funct_vector::*;
use vector_integral_method::*;
use vector_intergrator::*;

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use crate::integrator::FixedPrecision;
    use crate::vector_integral_method::VecRectangular;
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }

    //------------------- INTEGRAL OF POLYNOMIAL FUNCTION --------------------

    #[test] // test if rectangular_integral is in [ ( 1 - precision) exact_polynom_integral), (1 + precision) exact_polynom_integral]
    fn rectangular_on_polynomial() {
        let number_of_points = 100;
        let parameters = random_vector() ;
        let precision: f32 = 0.1;
        let polynomial = PolynomialFunction::new(parameters.to_vec());

        let integral = Rectangular{} ;

        assert!(!((1.0 - precision) * exact_polynom_integral(&parameters) >= integral.integrate(&polynomial, number_of_points)
            || integral.integrate(&polynomial, number_of_points) >= (1.0 + precision) * exact_polynom_integral(&parameters) ))
    }

    #[test]
    fn trapezoidal_on_polynomial() {
        let parameters = random_vector();
        let precision: f32 = 0.1;
        let number_of_point = 100;
        let polynomial = PolynomialFunction::new(parameters.to_vec());

        let integral = Trapezoiadal{} ;

        assert!(!((1.0 - precision) * exact_polynom_integral(&parameters) >= integral.integrate(&polynomial, number_of_point)
            || integral.integrate(&polynomial, number_of_point) >= (1.0 + precision) * exact_polynom_integral(&parameters) ))
    }

    #[test]
    fn simpson1_on_polynomial() {
        let parameters = random_vector();
        let precision: f32 = 0.1;
        let number_of_point = 100;
        let polynomial = PolynomialFunction::new(parameters.to_vec());

        let integral = Simpson1{} ;

        assert!(!((1.0 - precision) * exact_polynom_integral(&parameters) >= integral.integrate(&polynomial, number_of_point)
            || integral.integrate(&polynomial, number_of_point) >= (1.0 + precision) * exact_polynom_integral(&parameters) ))
    }


    #[test]
    fn simpson2_on_polynomial() {
        let parameters = random_vector();
        let precision: f32 = 0.1;
        let number_of_point = 300; // number_of_points has to be multiple of 3 !!
        let polynomial = PolynomialFunction::new(parameters.to_vec());

        let integral = Simpson2{} ;

        assert!(!((1.0 - precision) * exact_polynom_integral(&parameters) >= integral.integrate(&polynomial, number_of_point)
            || integral.integrate(&polynomial, number_of_point) >= (1.0 + precision) * exact_polynom_integral(&parameters) ))
    }


    //-------------------   RELATIVE ERRORS   -----------------

    #[test]
    fn error_comparison() {
        let parameters = random_vector();
        let precision : f32 = 0.01 ;
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        //let polynomial = Sin{i : 0} ;

        let integral_precision = FixedPrecision{precision};

        let rectangular = Rectangular{};
        let trapezoidal = Trapezoiadal{};
        let simpson1 = Simpson1{};
        let simpson2 = Simpson2{};
        let exact_integral = exact_polynom_integral(&parameters) ;
        //let exact_integral = - 1.0_f32.cos() + 1.0 ;


        //let vector: Vec<Box<dyn IntegralMethod>> = vec![Box::new(rectangular), Box::new(trapezoidal),
        //                                                Box::new(simpson1), Box::new(simpson2)];

        let rel_error_real_rectangular = ((integral_precision.integrate(&rectangular,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_trapezoidal = ((integral_precision.integrate(&trapezoidal,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_simpson1 = ((integral_precision.integrate(&simpson1,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_simpson2 = ((integral_precision.integrate(&simpson2,&polynomial) - exact_integral).abs() / exact_integral).abs() ;


        println!("Exact integral : {}", exact_integral);
        println!("The random polynomial function has order : {}", parameters.len());
        println!("Rectangular: {} , with {} points ( relative error bound : {}; real error : {})",
                 integral_precision.integrate(&rectangular,&polynomial),
                 integral_precision.number_of_points(&rectangular,&polynomial),
                 integral_precision.relative_error(&rectangular, &polynomial),
                 rel_error_real_rectangular);
        println!("Trapezoidal: {}, with {} points ( relative error bound : {}; real error : {})",
                 integral_precision.integrate(&trapezoidal,&polynomial),
                 integral_precision.number_of_points(&trapezoidal,&polynomial),
                 integral_precision.relative_error(&trapezoidal, &polynomial),
                 rel_error_real_trapezoidal);
        println!("Simpson1: {}, with {} points ( relative error bound : {}; real error : {})",
                 integral_precision.integrate(&simpson1,&polynomial),
                 integral_precision.number_of_points(&simpson1,&polynomial),
                 integral_precision.relative_error(&simpson1, &polynomial),
                 rel_error_real_simpson1);
        println!("Simpson2: {}, with {} points ( relative error bound : {}; real error : {})",
                 integral_precision.integrate(&simpson2,&polynomial),
                 integral_precision.number_of_points(&simpson2,&polynomial),
                 integral_precision.relative_error(&simpson2, &polynomial),
                 rel_error_real_simpson2);



    }

//----------------------- TEST OF FUNCTVECTORS


    #[test]
    fn funct_vector() {
        let parameters1 = random_vector();
        let parameters2 = random_vector();
        let parameters3 = random_vector();

        let functions_order = [ parameters1.len(),parameters2.len(),parameters3.len()];

        let polynomial1 = PolynomialFunction::new(parameters1.to_vec());
        let polynomial2 = PolynomialFunction::new(parameters2.to_vec());
        let polynomial3 = PolynomialFunction::new(parameters3.to_vec());


        let collection : FunctVector = FunctVector{ components : vec![Box::new(polynomial1),Box::new(polynomial2),Box::new(polynomial3)]};
        let rectangular = VecRectangular{};
        let trapezoidal = VecTrapezoiadal{};
        let simpson1 = VecSimpson1{};
        let simpson2 = VecSimpson2{};


        let precision :f32 = 0.01;
        let vec_integral_method = VecFixedPrecision{precision};

        let result = [vec_integral_method.integrate_uniform(&rectangular,&collection),
            vec_integral_method.integrate_uniform(&trapezoidal,&collection),
            vec_integral_method.integrate_uniform(&simpson1,&collection),
            vec_integral_method.integrate_uniform(&simpson2,&collection)];
        let number_of_points = [vec_integral_method.number_of_points_uniform(&rectangular,&collection),
            vec_integral_method.number_of_points_uniform(&trapezoidal,&collection),
            vec_integral_method.number_of_points_uniform(&simpson1,&collection),
            vec_integral_method.number_of_points_uniform(&simpson2,&collection)];
        let relative_errors = [vec_integral_method.relative_errors(&rectangular, &collection),
            vec_integral_method.relative_errors(&trapezoidal, &collection),
            vec_integral_method.relative_errors(&simpson1, &collection),
            vec_integral_method.relative_errors(&simpson2, &collection)];

        println!("The random polynomial functions have order : {:?}", functions_order);
        println!("Rectangular ; Trapezoidal ; Simpson1 ; Simpson2");
        println!("Number of point used : {:?} ",number_of_points);
        println!("Integral result: [first function, second function, third function]");
        println!("{:?}", result);
        println!("Relative errors bound:");
        println!("{:?}", relative_errors );

    }





}