pub mod functions;
pub mod integrals;

use functions::*;
use integrals::*;

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

    #[test] // test if rectangular_integral is in [ ( 1 - precision) exact_polynom_integral), (1 + precision) exact_polynom_integral]
    fn rectangular_on_polynomial() {
        let number_of_points = 100;
        let parameters = random_vector() ;
        let precision: f32 = 0.1;

        let integral = Rectangular{} ;
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        assert!(!((1.0 - precision) * exact_polynom_integral(&parameters) >= integral.integrate(&polynomial, number_of_points)
            || integral.integrate(&polynomial, number_of_points) >= (1.0 + precision) * exact_polynom_integral(&parameters) ))
    }


    #[test]
    fn trapezoidal_on_polynomial() {
        let parameters = random_vector();
        let precision: f32 = 0.1;
        let number_of_point = 100;
        let integral = Trapezoiadal{} ;
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        assert!(!((1.0 - precision) * exact_polynom_integral(&parameters) >= integral.integrate(&polynomial, number_of_point)
            || integral.integrate(&polynomial, number_of_point) >= (1.0 + precision) * exact_polynom_integral(&parameters) ))
    }

    #[test]
    fn simpson1_on_polynomial() {
        let parameters = random_vector();
        let precision: f32 = 0.1;
        let number_of_point = 100;
        let integral = Simpson1{} ;
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        assert!(!((1.0 - precision) * exact_polynom_integral(&parameters) >= integral.integrate(&polynomial, number_of_point)
            || integral.integrate(&polynomial, number_of_point) >= (1.0 + precision) * exact_polynom_integral(&parameters) ))
    }


    #[test]
    fn simpson2_on_polynomial() {
        let parameters = random_vector();
        let precision: f32 = 0.1;
        let number_of_point = 300; // number_of_points has to be multiple of 3 !!
        let integral = Simpson2{} ;
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        assert!(!((1.0 - precision) * exact_polynom_integral(&parameters) >= integral.integrate(&polynomial, number_of_point)
            || integral.integrate(&polynomial, number_of_point) >= (1.0 + precision) * exact_polynom_integral(&parameters) ))
    }

    #[test]
    fn error_bound_rectangular() {
        let number_of_points = 100;
        let parameters1 = &vec![1.1, 2.2];
        let number_of_points_float = number_of_points as f32;
        let line = PolynomialFunction::new(parameters1.to_vec());
        let derivative_line = line.constr_derivative();
        let mut max: f32 = derivative_line.evaluate(&0.0).abs();
        for i in 1..number_of_points {
            let j: f32 = i as f32;
            if derivative_line.evaluate(&(j / number_of_points_float)) > max {
                max = derivative_line.evaluate(&(j / number_of_points_float));
            }
        }
        let error: f32 = max / (2.0 * number_of_points_float);
        println!("The error_bound is {error}");
    }
/*
    #[test]
    fn random_vec() {
        let random_vector : Vec<f32> = random_vector() ;
        println!("{}", random_vector.len()) ;
        for i in random_vector {
            println!("{i}")
        }


    }


 */

}