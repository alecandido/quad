pub mod functions;
pub mod integral_method;
pub mod integrator;

use functions::*;
use integral_method::*;
use integrator::*;

pub fn add(left: usize, right: usize) -> usize {
    left + right
}

#[cfg(test)]
mod tests {
    use crate::integrator::FixedPrecision;
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
        let mut number_of_points_rec = 2;
        let mut number_of_points_trap = 2;
        let mut number_of_points_simp1 = 2;
        let mut number_of_points_simp2 = 3;
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        //let polynomial = Sin{i : 0} ;

        let integral_precision = FixedPrecision{ precision};


        let rectangular = Rectangular{};
        let trapezoidal = Trapezoiadal{};
        let simpson1 = Simpson1{};
        let simpson2 = Simpson2{};
        let exact_integral = exact_polynom_integral(&parameters) ;
        //let exact_integral = - 1.0_f32.cos() + 1.0 ;


        let mut rel_error_bound_rec: f32 = rectangular.relative_error(&polynomial, number_of_points_rec);
        let mut rel_error_bound_trap: f32 = trapezoidal.relative_error(&polynomial, number_of_points_trap);
        let mut rel_error_bound_simp1: f32 = simpson1.relative_error(&polynomial, number_of_points_simp1);
        let mut rel_error_bound_simp2: f32 = simpson2.relative_error(&polynomial, number_of_points_simp2);


        println!("Exact integral : {}", exact_integral);
        while rel_error_bound_rec > precision {
            number_of_points_rec += 2 ;
            rel_error_bound_rec = rectangular.relative_error(&polynomial,number_of_points_rec);
        }
        while rel_error_bound_trap > precision {
            number_of_points_trap += 2 ;
            rel_error_bound_trap = trapezoidal.relative_error(&polynomial,number_of_points_trap);
        }
        while rel_error_bound_simp1 > precision {
            number_of_points_simp1 += 2 ;
            rel_error_bound_simp1 = simpson1.relative_error(&polynomial,number_of_points_simp1);
        }
        while rel_error_bound_simp2 > precision {
            number_of_points_simp2 += 3 ;
            rel_error_bound_simp2 = simpson2.relative_error(&polynomial,number_of_points_simp2);
        }


        let rectangular_integral : f32 = rectangular.integrate(&polynomial, number_of_points_rec);
        let trapezoidal_integral : f32 = trapezoidal.integrate(&polynomial, number_of_points_trap);
        let simpson1_integral : f32 = simpson1.integrate(&polynomial, number_of_points_simp1);
        let simpson2_integral : f32 = simpson2.integrate(&polynomial, number_of_points_simp2);

        let rel_error_real_rectangular = (rectangular_integral - exact_integral).abs() / exact_integral ;
        let rel_error_real_trapezoidal = (trapezoidal_integral - exact_integral).abs() / exact_integral ;
        let rel_error_real_simpson1 = (simpson1_integral - exact_integral).abs() / exact_integral ;
        let rel_error_real_simpson2 = (simpson2_integral - exact_integral).abs() / exact_integral ;

        println!("The random polynomial function has order : {}", parameters.len());
        println!("Rectangular: {} and {} , with {} points ( relative error bound : {} ; real error : {})",
                 integral_precision.integrate(&rectangular,&polynomial),
                 rectangular_integral, number_of_points_rec, rel_error_bound_rec, rel_error_real_rectangular);
        println!("Trapezoidal: {} and {} , with {} points ( relative error bound : {} ; real error : {})",
                 integral_precision.integrate(&trapezoidal,&polynomial),
                 trapezoidal_integral, number_of_points_trap, rel_error_bound_trap, rel_error_real_trapezoidal);
        println!("Simpson1: {} and {} , with {} points ( relative error bound: {} ; real error : {})",
                 integral_precision.integrate(&simpson1,&polynomial),
                 simpson1_integral, number_of_points_simp1, rel_error_bound_simp1, rel_error_real_simpson1);
        println!("Simpson2: {} and {} , with {} points ( relative error bound: {} ; real error : {})",
                 integral_precision.integrate(&simpson2,&polynomial),
                 simpson2_integral, number_of_points_simp2, rel_error_bound_simp2, rel_error_real_simpson2);



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