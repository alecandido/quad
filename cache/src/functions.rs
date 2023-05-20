use rand::{thread_rng, Rng};
use rand::distributions::Standard;

pub trait Function {
    fn evaluate(&self, x : &f64) -> f64;
    fn abs_max(&self, number_of_points : i32) -> f64;
    fn first_derivative_abs_max(&self) -> f64;
    fn second_derivative_abs_max(&self) -> f64;
    fn fourth_derivative_abs_max(&self) -> f64;
}
// PolynomialFunction : a + bx + cx^2 + dx^3 + .... ; parameters = [a,b,c,d,...]
#[derive(Debug,Clone)]
pub struct PolynomialFunction {
    pub parameters : Vec<f64>,
}

impl PolynomialFunction {

    pub fn new(parameters : Vec<f64>) -> Self {
        Self{ parameters}
    }

    pub fn constr_derivative(&self) -> Self {
        let mut derivative_parameters : Vec<f64>;
        if self.parameters.len() == 1 {
            derivative_parameters = vec![0.0];
        }
        else {
            derivative_parameters = vec![0.0;self.parameters.len()-1];
        }
        for i in 0..self.parameters.len()-1 {
            let j : f64 = i as f64 ;
            derivative_parameters[i] = ( j + 1.0 ) * self.parameters[i+1] ;
        }
        Self{
            parameters : derivative_parameters,
        }
    }
}

impl Function for PolynomialFunction {
    fn evaluate(&self, x: &f64) -> f64 {
        let mut total : f64 = self.parameters[0] ;
        for i  in 1..self.parameters.len() {
            let n : i32 = i as i32 ;
            total += self.parameters[i] * x.powi(n) ;
        }
        total
    }


    fn abs_max(&self, number_of_points : i32) -> f64 {
        let mut max: f64 = self.evaluate(&0.0).abs();
        let number_of_points_float = number_of_points as f64;
        for i in 1..number_of_points {
            let j: f64 = i as f64;
            if self.evaluate(&(j / number_of_points_float)).abs() > max {
                max = self.evaluate(&(j / number_of_points_float)).abs();
            }
        }
        max
    }


    fn first_derivative_abs_max(&self) -> f64 {
        let derivative = self.constr_derivative();
        let first_derivative_abs_max: f64 = derivative.abs_max(1000);
        first_derivative_abs_max
    }

    fn second_derivative_abs_max(&self) -> f64 {
        let derivative = self.constr_derivative();
        let second_derivative = derivative.constr_derivative();
        let second_derivative_abs_max: f64 = second_derivative.abs_max(1000);
        second_derivative_abs_max
    }

    fn fourth_derivative_abs_max(&self) -> f64 {
        let derivative = self.constr_derivative();
        let second_derivative = derivative.constr_derivative();
        let third_derivative = second_derivative.constr_derivative();
        let fourth_derivative = third_derivative.constr_derivative();
        let fourth_derivative_abs_max: f64 = fourth_derivative.abs_max(1000);
        fourth_derivative_abs_max
    }


}

pub fn exact_polynom_integral( parameters : &Vec<f64>) -> f64 {
    let mut total : f64 = parameters[0] ;
    for i in 1..parameters.len() {
        let j : f64 = i as f64 ;
        total += parameters[i] / ( j + 1.0 );
    }
    total
}

// constructs a vector of random size (between 1 and 15) with random component
pub fn random_vector() -> Vec<f64> {
    let mut rng = thread_rng();
    let capacity : usize = rng.gen_range(20..25) ;
    let v: Vec<f64> = (&mut rng).sample_iter(Standard).take(capacity).collect();
    v
}


// i = 0 : sin ;
// i = 1 : cos ;
// i must be positive

pub struct Sin {
    pub i: i32,
}

impl Sin {
    pub fn constr_derivative(&self) -> Self {
        Self{ i : self.i + 1}
    }
}

impl Function for Sin{
    fn evaluate(&self, x : &f64) -> f64 {
        match self.i%4{
            0 => x.sin(),
            1 => x.cos(),
            2 => - x.sin(),
            3 => - x.cos(),
            _ => 0.0,
        }
    }


    fn abs_max(&self, _number_of_points: i32) -> f64 {
        if self.i % 2 == 1 { 1.0 }
        else { 1.0_f64.sin() }
    }

    fn first_derivative_abs_max(&self) -> f64 {
        let derivative = self.constr_derivative();
        derivative.abs_max(1000)
    }

    fn second_derivative_abs_max(&self) -> f64 {
        let derivative = self.constr_derivative();
        let second_derivative = derivative.constr_derivative();
        second_derivative.abs_max(1000)
    }

    fn fourth_derivative_abs_max(&self) -> f64 {
        let derivative = self.constr_derivative();
        let second_derivative = derivative.constr_derivative();
        let third_derivative = second_derivative.constr_derivative();
        let fourth_derivative = third_derivative.constr_derivative();
        fourth_derivative.abs_max(1000)
    }

}

pub fn exact_sin_integral() -> f64 { - 1.0_f64.cos() + 1.0 }

// e^(ax)
pub struct Exp {
    pub exponent : f64,
}

impl Function for Exp{
    fn evaluate(&self, x : &f64) -> f64 {
        ( self.exponent * x ).exp()
    }


    fn abs_max(&self, _number_of_points: i32) -> f64 {
        let result: f64;

        if self.exponent > 0.0 { result = self.exponent.exp(); } else { result = 1.0; }
        result
    }

    fn first_derivative_abs_max(&self) -> f64 {
        self.abs_max(1000) * self.exponent
    }

    fn second_derivative_abs_max(&self) -> f64 {
        self.abs_max(1000) * self.exponent.powi(2)
    }

    fn fourth_derivative_abs_max(&self) -> f64 {
        self.abs_max(1000) * self.exponent.powi(4)
    }

}

pub fn exact_exp_integral(exponential : f64) -> f64{
    if exponential != 0.0 { (exponential.exp() - 1.0)/exponential }
    else {0.0}
}


