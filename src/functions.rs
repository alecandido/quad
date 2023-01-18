pub trait Function {
    fn evaluate(&self, x : &f32) -> f32;
}

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


pub struct Sin{}

impl Function for Sin{
    fn evaluate(&self, x : &f32) -> f32 {
        x.sin()
    }
}