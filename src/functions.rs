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
        for i in 1..self.parameters.len(){
            let mut contr = self.parameters[i] ;
            for _n in 1..i+1 {
                contr = contr * x ;
            }
            total += contr ;
        }
        total
    }
    fn constr_derivative(&self) -> Self {
        let mut derivative_parameters : Vec<f32> = vec![0.0;self.parameters.len()-1] ;
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

/*
pub struct Sin{}
pub struct Cos{}

impl Function for Sin{
    fn evaluate(&self, x : &f32) -> f32 {
        x.sin()
    }
    fn constr_derivative(&self) -> Self {
        Cos{}
    }
}

impl Function for Cos{
    fn evaluate(&self, x : &f32) -> f32 {
        x.cos()
    }
    fn constr_derivative(&self) -> Self {
        Sin{}
    }
}

*/

/*
pub struct Line{
    pub a : f32,
    pub b : f32,
}



impl Function for Line {
    fn evaluate(&self, x : &f32) -> f32 {
        self.a + self.b * x
    }
}

pub struct Parabola{
    pub a : f32,
    pub b : f32,
    pub c : f32,
}

impl Function for Parabola{
    fn evaluate(&self, x : &f32) -> f32 {
        self.a + self.b * x + self.c * x * x
    }
}

pub struct Cubic{
    pub a : f32,
    pub b : f32,
    pub c : f32,
    pub d : f32,
}

impl Function for Cubic{
    fn evaluate(&self, x : &f32) -> f32 {
        self.a + self.b * x + self.c * x * x + self.d * x * x * x
    }
}

*/