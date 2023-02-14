use crate::*;
use std::sync::{Arc, Mutex};

pub struct FunctVector {
    pub components: Vec<Box<dyn Function + Send + Sync >>,
}
pub struct FunVector {
    pub components : Vec<Arc<Mutex<dyn Function + Send + Sync>>>,
}

pub struct FnVec {
    pub components : Vec<Box<dyn Fn(f64)->f64 + Send + Sync >>,
}

pub struct FnVecP {
    pub components : Vec<Arc<Mutex<dyn Fn(f64)->f64 + Send + Sync>>>,
}

pub struct IntVec{
    pub components : Vec<Arc<Mutex<dyn QuadIntegralMethod + Send + Sync>>>
}


pub fn random_polynomials3() -> FunctVector{
        FunctVector{
            components : vec![Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                              Box::new(PolynomialFunction::new(random_vector()))]
        }
    }

pub fn random_polynomials4() -> FunctVector{
    FunctVector{
        components : vec![Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector()))]
    }
}

pub fn random_polynomials5() -> FunctVector{
    FunctVector{
        components : vec![Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector()))]
    }
}

pub fn random_polynomials6() -> FunctVector{
    FunctVector{
        components : vec![Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector()))]
    }
}

pub fn random_polynomials7() -> FunctVector{
    FunctVector{
        components : vec![Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector()))]
    }
}

pub fn random_polynomials8() -> FunctVector{
    FunctVector{
        components : vec![Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector()))]
    }
}

pub fn random_polynomials9() -> FunctVector{
    FunctVector{
        components : vec![Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector()))]
    }
}

pub fn random_polynomials10() -> FunctVector{
    FunctVector{
        components : vec![Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector())),
                          Box::new(PolynomialFunction::new(random_vector())),Box::new(PolynomialFunction::new(random_vector()))]
    }
}

pub fn random3_polynomials() -> FunVector{
    FunVector{
        components : vec![Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector())))]
        }
}

pub fn random6_polynomials() -> FunVector{
    FunVector{
        components : vec![Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector())))]
    }
}

pub fn random9_polynomials() -> FunVector{
    FunVector{
        components : vec![Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector()))),
                          Arc::new(Mutex::new(PolynomialFunction::new(random_vector())))]
    }
}