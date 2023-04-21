use crate::*;
use std::sync::{Arc, Mutex};
use crate::quad_integral_method::*;


pub struct FunctVector {
    pub components: Vec<Box<dyn Function + Send + Sync >>,
}
#[derive(Clone)]
pub struct FunVector {
    pub components : Vec<Arc<Mutex<dyn Function + Send + Sync>>>,
}

pub struct FnVec {
    pub components : Vec<Box<dyn Fn(f64)->f64 + Send + Sync >>,
}

pub struct FnVec4 {
    pub components : [Box<dyn Fn(f64)->f64 + Send + Sync>;4]
}

pub struct FnVec8 {
    pub components : [Box<dyn Fn(f64)->f64 + Send + Sync>;8]
}

pub struct FnVec3 {
    pub components : [Box<dyn Fn(f64)->f64 + Send + Sync>;3]
}

pub struct FnVecRayon {
    pub components : Vec<Box<dyn Fn(f64)->f64 + Send + Sync >>,
}

#[derive(Clone)]
pub struct FnVecP {
    pub components : Vec<Arc<Mutex<dyn Fn(f64)->f64 + Send + Sync>>>,
}

#[derive(Clone)]
pub struct FnVecPa {
    pub components : Vec<Arc<dyn Fn(f64)->f64 + Send + Sync>>,
}

#[derive(Clone)]
pub struct FnPa {
    pub components : Arc<dyn Fn(f64)->f64 + Send + Sync>,
}

#[derive(Clone)]
pub struct ResVecPa {
    pub components : Vec<Arc<Mutex<qage_1dvec::Result>>>,
}
impl ResVecPa{
    pub fn push(&mut self,new_component: Arc<Mutex<qage_1dvec::Result>> ){
        self.components.push(new_component);
    }
    pub fn new_empty() -> Self{
        Self{components : vec![]}
    }
}
unsafe impl Send for ResVecPa{}
unsafe impl Sync for ResVecPa{}

#[derive(Clone)]
pub struct IntVec{
    pub components : Vec<Arc<Mutex<dyn QuadIntegralMethod + Send + Sync>>>
}



pub fn random_polynomials(number : i32) -> (FnVec,FnVecP){

        let mut comp : Vec<Box<dyn Fn(f64)->f64 + Send + Sync >> = vec![];
        let mut comp2 : Vec<Arc<Mutex<dyn Fn(f64)->f64 + Send + Sync>>> = vec![];
        for _i in 0..number{
            let f = |x:f64| PolynomialFunction::new(random_vector()).evaluate(&x);
            comp2.push(Arc::new(Mutex::new(f)));
            comp.push(Box::new(f));
        }
    let fun = FnVec{
        components : comp,
    };
    let fun2 = FnVecP{
        components : comp2,
    };
    (fun,fun2)
}

pub fn random_polynomials2(number : i32) -> (FnVec,FnVecP){

    let mut comp : Vec<Box<dyn Fn(f64)->f64 + Send + Sync >> = vec![];
    let mut comp2 : Vec<Arc<Mutex<dyn Fn(f64)->f64 + Send + Sync>>> = vec![];
    for _i in 0..number {
        let parameters = random_vector();
        let closure = move |x: f64| {
            let mut res = 0.0;
            for n in 0..parameters.len() {
                res += parameters[n] * x.powi(n as i32);
            }
            res
        };
        comp2.push(Arc::new(Mutex::new(closure.clone())));
        comp.push(Box::new(closure));
    }
        let fun = FnVec {
            components: comp,
        };
        let fun2 = FnVecP {
            components: comp2,
        };
    (fun,fun2)
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
