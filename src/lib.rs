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
    use super::*;

    #[test]
    fn it_works() {
        let result = add(2, 2);
        assert_eq!(result, 4);
    }


    //-------------------   RELATIVE ERRORS   -----------------

    #[test]
    fn error_comparison() {
        let test : bool = false; // set to true if you want to have the results printed out
        let precision : f32 = 0.1 ; // set the relative errors required

        let parameters = random_vector();
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        //let polynomial = Sin{i : 0} ;
        let integral_precision = FixedPrecision{precision};
        let rectangular = Rectangular{};
        let trapezoidal = Trapezoiadal{};
        let simpson1 = Simpson1{};
        let simpson2 = Simpson2{};
        let exact_integral = exact_polynom_integral(&parameters) ;
        //let exact_integral = - 1.0_f32.cos() + 1.0 ;


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

        assert!(rel_error_real_rectangular < precision || rel_error_real_trapezoidal < precision
                || rel_error_real_simpson1 < precision || rel_error_real_simpson2 < precision)

    }

//----------------------- TEST OF FUNCTVECTORS


    #[test]
    fn funct_vector() {
        let test : bool = false; // set to true if you want to have the results printed out
        let precision :f32 = 0.1; // set the relative errors required

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
            println!("UNIFORM:");
            println!("The random polynomial functions have order : {:?}", functions_order);
            println!("Notation : [Rectangular, Trapezoidal, Simpson1, Simpson2");
            println!("Number of point used : {:?} ", number_of_points);
            println!("Integral result: [first function, second function, third function]");
            println!("{:?}", result);
            println!("Relative errors bound:");
            println!("{:?}", relative_errors);

            println!("NON-UNIFORM:");
            println!("The random polynomial functions have order : {:?}", functions_order);
            println!("Notation : [Rectangular, Trapezoidal, Simpson1, Simpson2");
            println!("Number of point used : {:?} ", number_of_points_nu);
            println!("Integral result: [first function, second function, third function]");
            println!("{:?}", result_nu);
            println!("Relative errors bound:");
            println!("{:?}", relative_errors_nu);
        }

        assert!(rel_error_real[0][0] < precision || rel_error_real[0][1] < precision ||
            rel_error_real[0][2] < precision ||
            rel_error_real[1][0] < precision || rel_error_real[1][1] < precision ||
            rel_error_real[1][2] < precision ||
            rel_error_real[2][0] < precision || rel_error_real[2][1] < precision ||
            rel_error_real[2][2] < precision ||
            rel_error_real[3][0] < precision || rel_error_real[3][1] < precision ||
            rel_error_real[3][2] < precision ||
            rel_error_real_nu[0][0] < precision || rel_error_real_nu[0][1] < precision ||
            rel_error_real_nu[0][2] < precision ||
            rel_error_real_nu[1][0] < precision || rel_error_real_nu[1][1] < precision ||
            rel_error_real_nu[1][2] < precision ||
            rel_error_real_nu[2][0] < precision || rel_error_real_nu[2][1] < precision ||
            rel_error_real_nu[2][2] < precision ||
            rel_error_real_nu[3][0] < precision || rel_error_real_nu[3][1] < precision ||
            rel_error_real_nu[3][2] < precision)
        }





}