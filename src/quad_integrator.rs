use crate::*;
use std::thread;
use std::time::Instant;

pub trait QuadIntegrator{
    fn integrate(&self, integral_method : &impl QuadIntegralMethod, fun : & FnVec) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>) ;
    fn integrate_p(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>) ;
}

pub struct QngIntegrator {
    pub a : Vec<f64>,
    pub b : Vec<f64>,
    pub epsabs : Vec<f64>,
    pub epsrel : Vec<f64>,
}


impl QuadIntegrator for QngIntegrator {
    fn integrate(&self, integral_method : &impl QuadIntegralMethod, fun : & FnVec) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>){
        let mut res : Vec<f64> = vec![];
        let mut err : Vec<f64> = vec![];
        let mut number : Vec<i32> = vec![];
        let mut controll : Vec<i32> = vec![];
        let dimension = fun.components.len();
        let start = Instant::now();
        for k in 0..dimension{
            let time1 = Instant::now();
            let partial_result = integral_method.integrate(&fun.components[k],self.a[k],self.b[k],self.epsabs[k],self.epsrel[k]);
            res.push(partial_result.0);
            err.push(partial_result.1);
            number.push(partial_result.2);
            controll.push(partial_result.3);
            let time2 = Instant::now();
            println!("{:?}",time2-time1);
        }
        let end = Instant::now();
        println!("{:?}",end-start);
        (res,err,number,controll)

    }

    fn integrate_p(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>){
        let mut res : Vec<f64> = vec![];
        let mut err : Vec<f64> = vec![];
        let mut number : Vec<i32> = vec![];
        let mut controll : Vec<i32> = vec![];
        let mut handles = vec![];
        let dimension = fun.components.len();
        let start = Instant::now();
        for k in 0..dimension{
            let time1 = Instant::now();
            let pointer = fun.components[k].clone();
            let pointer2 = integral_method.components[k].clone();
            let a = self.a[k];
            let b = self.b[k];
            let epsabs = self.epsabs[k];
            let epsrel = self.epsrel[k];
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let integral = pointer2.lock().unwrap();
                integral.integrate(&*function,a,b,epsabs,epsrel)
            });
            let time2 = Instant::now();
            println!("{:?}",time2-time1);
            handles.push(handle);
        }

        for handle in handles {
            let time3 = Instant::now();
            let results = handle.join().unwrap();
            res.push(results.0);
            err.push(results.1);
            number.push(results.2);
            controll.push(results.3);
            let time4 = Instant::now();
            println!("{:?}",time4-time3);
        }

        let end = Instant::now();
        println!("{:?}",end-start);
        (res,err,number,controll)
    }
}
