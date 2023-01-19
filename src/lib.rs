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
        let parameters1 = &vec![1.1,2.2] ;
        let parameters2 = &vec![3.4,6.5,7.8] ;
        let parameters3 = &vec![7.6,1.9,9.1,11.2] ;
        let precision : f32 = 0.1 ;
        let number_of_point = 100 ;
        let line = PolynomialFunction{ parameters : parameters1.to_vec()} ;
        let parabola = PolynomialFunction{ parameters : parameters2.to_vec()} ;
        let cubic = PolynomialFunction{ parameters : parameters3.to_vec()} ;
        assert!(!( (1.0 - precision) * exact_polynom_integral(parameters1) >= rectangular_integral(&line, number_of_point)
            || rectangular_integral(&line, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters1)
            || (1.0 - precision) * exact_polynom_integral(parameters2) >= rectangular_integral(&parabola, number_of_point)
            || rectangular_integral(&parabola, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters2)
            || (1.0 - precision) * exact_polynom_integral(parameters3) >= rectangular_integral(&cubic, number_of_point)
            || rectangular_integral(&cubic, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters3)))
    }
}

    #[test]
    fn trapezoidal_on_polynomial() {
        let parameters1 = &vec![1.1,2.2] ;
        let parameters2 = &vec![3.4,6.5,7.8] ;
        let parameters3 = &vec![7.6,1.9,9.1,11.2] ;
        let precision : f32 = 0.1 ;
        let number_of_point = 100 ;
        let line = PolynomialFunction{ parameters : parameters1.to_vec()} ;
        let parabola = PolynomialFunction{ parameters : parameters2.to_vec()} ;
        let cubic = PolynomialFunction{ parameters : parameters3.to_vec()} ;
        assert!(!( (1.0 - precision) * exact_polynom_integral(parameters1) >= trapezoidal_integral(&line, number_of_point)
            || trapezoidal_integral(&line, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters1)
            || (1.0 - precision) * exact_polynom_integral(parameters2) >= trapezoidal_integral(&parabola, number_of_point)
            || trapezoidal_integral(&parabola, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters2)
            || (1.0 - precision) * exact_polynom_integral(parameters3) >= trapezoidal_integral(&cubic, number_of_point)
            || trapezoidal_integral(&cubic, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters3)))
    }

    #[test]
    fn simpson1_on_polynomial() {
        let parameters1 = &vec![1.1,2.2] ;
        let parameters2 = &vec![3.4,6.5,7.8] ;
        let parameters3 = &vec![7.6,1.9,9.1,11.2] ;
        let precision : f32 = 0.1 ;
        let number_of_point = 100 ;
        let line = PolynomialFunction{ parameters : parameters1.to_vec()} ;
        let parabola = PolynomialFunction{ parameters : parameters2.to_vec()} ;
        let cubic = PolynomialFunction{ parameters : parameters3.to_vec()} ;
        assert!(!( (1.0 - precision) * exact_polynom_integral(parameters1) >= simpson1(&line, number_of_point)
            || simpson1(&line, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters1)
            || (1.0 - precision) * exact_polynom_integral(parameters2) >= simpson1(&parabola, number_of_point)
            || simpson1(&parabola, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters2)
            || (1.0 - precision) * exact_polynom_integral(parameters3) >= simpson1(&cubic, number_of_point)
            || simpson1(&cubic, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters3)))
    }


    #[test]
    fn simpson2_on_polynomial() {
        let parameters1 = &vec![1.1,2.2] ;
        let parameters2 = &vec![3.4,6.5,7.8] ;
        let parameters3 = &vec![7.6,1.9,9.1,11.2] ;
        let precision : f32 = 0.1 ;
        let number_of_point = 100 ;
        let line = PolynomialFunction{ parameters : parameters1.to_vec()} ;
        let parabola = PolynomialFunction{ parameters : parameters2.to_vec()} ;
        let cubic = PolynomialFunction{ parameters : parameters3.to_vec()} ;
        assert!(!( (1.0 - precision) * exact_polynom_integral(parameters1) >= simpson2(&line, number_of_point)
            || simpson2(&line, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters1)
            || (1.0 - precision) * exact_polynom_integral(parameters2) >= simpson2(&parabola, number_of_point)
            || simpson2(&parabola, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters2)
            || (1.0 - precision) * exact_polynom_integral(parameters3) >= simpson2(&cubic, number_of_point)
            || simpson2(&cubic, number_of_point) >= (1.0 + precision) * exact_polynom_integral(parameters3)))
    }