use rand::{thread_rng, Rng};
use rand::distributions::Standard;

pub trait Function {
    fn evaluate(&self, x : &f32) -> f32 ;
    fn constr_derivative(&self) -> Self ;
    fn abs_max(&self, number_of_points : i32) -> f32 ;
}
// PolynomialFunction : a + bx + cx^2 + dx^3 + ....
pub struct PolynomialFunction {
    pub parameters : Vec<f32>,
}

impl PolynomialFunction {
    pub fn new(parameters : Vec<f32>) -> Self {
        Self{ parameters}
    }
}

impl Function for PolynomialFunction {
    fn evaluate(&self, x: &f32) -> f32 {
        let mut total : f32 = self.parameters[0] ;
        for i  in 1..self.parameters.len() {
            // let mut contr = self.parameters[i] ;
            // for _n in 1..i+1 {
            //    contr = contr * x ;
            // }
            let n : i32 = i as i32 ;
            total += self.parameters[i] * x.powi(n) ;
        }
        total
    }

    fn constr_derivative(&self) -> Self {
        let mut derivative_parameters : Vec<f32>;
        if self.parameters.len() == 1 {
            derivative_parameters = vec![0.0];
        }
        else {
            derivative_parameters = vec![0.0;self.parameters.len()-1];
        }
        for i in 0..self.parameters.len()-1 {
            let j : f32 = i as f32 ;
            derivative_parameters[i] = ( j + 1.0 ) * self.parameters[i+1] ;
        }
        Self{
            parameters : derivative_parameters,
        }
    }

    fn abs_max(&self, number_of_points : i32) -> f32 {
        let mut max: f32 = self.evaluate(&0.0).abs();
        let number_of_points_float = number_of_points as f32;
        for i in 1..number_of_points {
            let j: f32 = i as f32;
            if self.evaluate(&(j / number_of_points_float)).abs() > max {
                max = self.evaluate(&(j / number_of_points_float)).abs();
            }
        }
        max
    }

}

pub fn exact_polynom_integral( parameters : &Vec<f32>) -> f32 {
    let mut total : f32 = parameters[0] ;
    for i in 1..parameters.len() {
        let j : f32 = i as f32 ;
        total += parameters[i] / ( j + 1.0 );
    }
    total
}

pub fn random_vector() -> Vec<f32> {
    let mut rng = thread_rng();
    let capacity : usize = rng.gen_range(1..15) ;
    let v: Vec<f32> = (&mut rng).sample_iter(Standard).take(capacity).collect();
    v
}



pub struct Sin {
    pub i: i32,
}

impl Function for Sin{
    fn evaluate(&self, x : &f32) -> f32 {
        match self.i%4{
            0 => x.sin(),
            1 => x.cos(),
            2 => - x.sin(),
            3 => - x.cos(),
            _ => 0.0,
        }
    }
    fn constr_derivative(&self) -> Self {
        Self{ i : self.i + 1}
    }

    fn abs_max(&self, _number_of_points: i32) -> f32 {
        if self.i % 2 == 1 { 1.0 }
        else { 1.0_f32.sin() }
    }


}

/*
impl Function for Cos{
    fn evaluate(&self, x : &f32) -> f32 {
        x.cos()
    }
    fn constr_derivative(&self) -> Self {
        Sin{}
    }
}

*/
