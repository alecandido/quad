use crate::*;
use std::thread;
//  use std::time::Instant;
use threadpool::ThreadPool;
use crate::quad_integral_method::*;
use crate::quad_integrator_result::QuadIntegratorResult;

pub struct PIntegrator {
    pub a : f64,
    pub b : f64,
    pub epsabs : f64,
    pub epsrel : f64,
    pub integral_method: IntVec,
}

pub struct SIntegrator {
    pub a : f64,
    pub b : f64,
    pub epsabs : f64,
    pub epsrel : f64,
    pub integral_method: Box<dyn QuadIntegralMethod>,
}


impl SIntegrator {
    pub fn integrate(&self, fun: &FnVec) -> Vec<QuadIntegratorResult> {
        let mut res: Vec<QuadIntegratorResult> = vec![];
        let iter = fun.components.iter();
        for f in iter {
            let partial_result = self.integral_method.integrate(&f, self.a, self.b, self.epsabs, self.epsrel);
            res.push(partial_result);
        }
        res
    }
}

impl PIntegrator{

    pub fn integrate_rayon(&self, fun : & FnVecP) -> Vec<QuadIntegratorResult>{
        let dimension = fun.components.len();
        let res = Arc::new(Mutex::new(vec![QuadIntegratorResult::new_error(result_state::ResultState::Failure);dimension]));
        rayon::scope(|s|{
        for k in 0..dimension{
            let pointer = fun.components[k].clone();
            let pointer2 = self.integral_method.components[k].clone();
            let res = res.clone();
            let a = self.a;
            let b = self.b;
            let epsabs = self.epsabs;
            let epsrel = self.epsrel;
            s.spawn( move |_| {
                let function = pointer.lock().unwrap();
                let integral = pointer2.lock().unwrap();
                let results = integral.integrate(&*function,a,b,epsabs,epsrel);
                res.lock().unwrap()[k] = results;
            });
        }});
        let res = res.lock().unwrap().clone();
        res
    }

    pub fn integrate_handles(&self, fun : & FnVecP) -> Vec<QuadIntegratorResult>{
        let mut res : Vec<QuadIntegratorResult> = vec![];
        let mut handles = vec![];
        let dimension = fun.components.len();
        for k in 0..dimension{
            let pointer = fun.components[k].clone();
            let pointer2 = self.integral_method.components[k].clone();
            let a = self.a;
            let b = self.b;
            let epsabs = self.epsabs;
            let epsrel = self.epsrel;
            let handle = thread::spawn( move || {
                let function = pointer.lock().unwrap();
                let integral = pointer2.lock().unwrap();
                let res = integral.integrate(&*function,a,b,epsabs,epsrel);
                res

            });
            handles.push(handle);
        }

        for handle in handles {
            let results = handle.join().unwrap();
            res.push(results);
        }

        res
    }


    pub fn integrate_pool(&self, fun : & FnVecP) -> Vec<QuadIntegratorResult>{
        let dimension = fun.components.len();
        let res = Arc::new(Mutex::new(vec![QuadIntegratorResult::new_error(result_state::ResultState::Failure);dimension]));
        let pool = ThreadPool::new(dimension);
        for k in 0..dimension{
            let pointer = fun.components[k].clone();
            let pointer2 = self.integral_method.components[k].clone();
            let res = res.clone();
            let a = self.a;
            let b = self.b;
            let epsabs = self.epsabs;
            let epsrel = self.epsrel;
            pool.execute( move || {
                let function = pointer.lock().unwrap();
                let integral = pointer2.lock().unwrap();
                let results = integral.integrate(&*function,a,b,epsabs,epsrel);
                res.lock().unwrap()[k] = results;
            });
        }

        pool.join();
        let res = res.lock().unwrap().clone();
        res
    }


}

