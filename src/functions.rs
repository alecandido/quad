pub trait Function {
    fn evaluate(&self, x : &f32) -> f32;
}
// PolynomialFunction : a + bx + cx^2 + dx^3 + ....
pub struct PolynomialFunction {
    pub parameters : Vec<f32>,
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
}




pub struct Sin{}

impl Function for Sin{
    fn evaluate(&self, x : &f32) -> f32 {
        x.sin()
    }
}





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