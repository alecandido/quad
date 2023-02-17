use crate::*;
use std::thread;
//  use std::time::Instant;
use threadpool::ThreadPool;

pub trait QuadIntegrator{
    fn integrate(&self, integral_method : &impl QuadIntegralMethod, fun : & FnVec) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>) ;
    fn integrate_rayon(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>);
    fn integrate_p(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>) ;
    fn integrate_p2(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>) ;
    fn integrate_p3(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>) ;
}

pub struct QngIntegrator {
    pub a : f64,
    pub b : f64,
    pub epsabs : f64,
    pub epsrel : f64,
}


impl QuadIntegrator for QngIntegrator {
    fn integrate(&self, integral_method : &impl QuadIntegralMethod, fun : & FnVec) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>){
        let mut res : Vec<f64> = vec![];
        let mut err : Vec<f64> = vec![];
        let mut number : Vec<i32> = vec![];
        let mut controll : Vec<i32> = vec![];
        //  let start = Instant::now();
        let iter = fun.components.iter();
        //for k in 0..dimension{
        for f in iter{
            //  let time1 = Instant::now();
            let partial_result = integral_method.integrate(&f,self.a,self.b,self.epsabs,self.epsrel);
            res.push(partial_result.0);
            err.push(partial_result.1);
            number.push(partial_result.2);
            controll.push(partial_result.3);
            //  let time2 = Instant::now();
            //println!("{k} function : {:?}",time2-time1);
        }
        //  let end = Instant::now();
        //println!("Tot : {:?}",end-start);
        (res,err,number,controll)

    }

    fn integrate_rayon(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>){
        let dimension = fun.components.len();
        let res = Arc::new(Mutex::new(vec![0.0;dimension]));
        let err = Arc::new(Mutex::new(vec![0.0;dimension]));
        let number = Arc::new(Mutex::new(vec![0;dimension]));
        let controll = Arc::new(Mutex::new(vec![0;dimension]));
        //let start = Instant::now();
        rayon::scope(|s|{
        for k in 0..dimension{
            //  let time_forstart = Instant::now();
            //  let time_forstart2 = Instant::now();
            let pointer = fun.components[k].clone();
            let pointer2 = integral_method.components[k].clone();
            let res = res.clone();
            let err = err.clone();
            let number = number.clone();
            let controll = controll.clone();
            let a = self.a;
            let b = self.b;
            let epsabs = self.epsabs;
            let epsrel = self.epsrel;
            //  let time_inizvar = Instant::now();
            //  let j = k;
            s.spawn( move |_| {
                //  let time_poolstart = Instant::now();
                //  println!("{j}Thread opened in : {:?}", time_poolstart - time_forstart2);
                let function = pointer.lock().unwrap();
                let integral = pointer2.lock().unwrap();
                //  let time_poollock = Instant::now();
                //  println!("{j}Thread lock in : {:?}", time_poollock - time_forstart2);
                let results = integral.integrate(&*function,a,b,epsabs,epsrel);
                //  let time_poolint = Instant::now();
                //  println!("{j}Thread integrate in : {:?}", time_poolint - time_forstart2);
                res.lock().unwrap()[k] = results.0;
                err.lock().unwrap()[k] = results.1;
                number.lock().unwrap()[k] = results.2;
                controll.lock().unwrap()[k] = results.3;
                //  let time_poolresult = Instant::now();
                //  println!("{j} Thread end in : {:?}\n{:?},{:?}", time_poolresult - time_forstart2,rayon::current_num_threads(),rayon::current_thread_index());
            });
            //  let time_forend = Instant::now();
            //  println!("{k} For : Iniz : {:?}, Pool open : {:?}, Tot : {:?}",
            //          time_forstart - time_inizvar, time_forend - time_inizvar, time_forend - time_forstart);
        }});
        //  let time_forafter = Instant::now();
        //  let time_poolend = Instant::now();
        let res = res.lock().unwrap().clone();
        let err = err.lock().unwrap().clone();
        let number = number.lock().unwrap().clone();
        let controll = controll.lock().unwrap().clone();
        //  let time_end = Instant::now();
        //  println!("Join pool: {:?}, wrap result {:?}, Tot : {:?}", time_poolend - time_forafter,
        //          time_end - time_poolend, time_end - time_forafter);

        (res,err,number,controll)
    }

    fn integrate_p(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>){
        let mut res : Vec<f64> = vec![];
        let mut err : Vec<f64> = vec![];
        let mut number : Vec<i32> = vec![];
        let mut controll : Vec<i32> = vec![];
        let mut handles = vec![];
        let dimension = fun.components.len();
        //  let start = Instant::now();
        for k in 0..dimension{
            //  let time1 = Instant::now();
            let pointer = fun.components[k].clone();
            let pointer2 = integral_method.components[k].clone();
            //  let time11 = Instant::now();
            let a = self.a;
            let b = self.b;
            let epsabs = self.epsabs;
            let epsrel = self.epsrel;
            //  let time12 = Instant::now();
            let handle = thread::spawn( move || {
                //println!("thread creato");
                //  let time21 = Instant::now();
                let function = pointer.lock().unwrap();
                let integral = pointer2.lock().unwrap();
                //  let time22 = Instant::now();
                let res = integral.integrate(&*function,a,b,epsabs,epsrel);
                //  let time23 = Instant::now();
                //println!("in thread : {:?},{:?},{:?}",time22 - time21,time23-time22,time23-time21);
                res

            });
            //  let time13 = Instant::now();
            handles.push(handle);
            //  let time2 = Instant::now();
            //println!("{k} inizializz : {:?},{:?},{:?},{:?},{:?}",time11 - time1,time12-time11,time13-time12,time2-time13,time2-time1);
        }

        for handle in handles {
            //  let time3 = Instant::now();
            let results = handle.join().unwrap();
            //  let time4 = Instant::now();
            res.push(results.0);
            //  let time5 = Instant::now();
            //println!("unwrap : {:?} \n push : {:?}",time4 - time3,time5-time4);
            err.push(results.1);
            number.push(results.2);
            controll.push(results.3);

        }

        //  let end = Instant::now();
        //  println!("{:?}",end-start);
        (res,err,number,controll)
    }


    fn integrate_p2(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>){
        let dimension = fun.components.len();
        let res = Arc::new(Mutex::new(vec![0.0;dimension]));
        let err = Arc::new(Mutex::new(vec![0.0;dimension]));
        let number = Arc::new(Mutex::new(vec![0;dimension]));
        let controll = Arc::new(Mutex::new(vec![0;dimension]));
        let pool = ThreadPool::new(dimension);
        //  let start = Instant::now();
        for k in 0..dimension{
            //  let time_forstart = Instant::now();
            //  let time_forstart2 = Instant::now();
            let pointer = fun.components[k].clone();
            let pointer2 = integral_method.components[k].clone();
            let res = res.clone();
            let err = err.clone();
            let number = number.clone();
            let controll = controll.clone();
            let a = self.a;
            let b = self.b;
            let epsabs = self.epsabs;
            let epsrel = self.epsrel;
            //  let time_inizvar = Instant::now();
            //  let j = k;
            pool.execute( move || {
                //  let time_poolstart = Instant::now();
                //  println!("{j}Thread opened in : {:?}", time_poolstart - time_forstart2);
                let function = pointer.lock().unwrap();
                let integral = pointer2.lock().unwrap();
                //  let time_poollock = Instant::now();
                //  println!("{j}Thread lock in : {:?}", time_poollock - time_forstart2);
                let results = integral.integrate(&*function,a,b,epsabs,epsrel);
                //  let time_poolint = Instant::now();
                //  println!("{j}Thread integrate in : {:?}", time_poolint - time_forstart2);
                res.lock().unwrap()[k] = results.0;
                err.lock().unwrap()[k] = results.1;
                number.lock().unwrap()[k] = results.2;
                controll.lock().unwrap()[k] = results.3;
                //  let time_poolresult = Instant::now();
                //  println!("{j} Thread end in : {:?}", time_poolresult - time_forstart2);
            });
            //  let time_forend = Instant::now();
            //  println!("{k} For : Iniz : {:?}, Pool open : {:?}, Tot : {:?}",
            //          time_forstart - time_inizvar, time_forend - time_inizvar, time_forend - time_forstart);
        }

        //  let time_forafter = Instant::now();
        pool.join();
        //  let time_poolend = Instant::now();
        let res = res.lock().unwrap().clone();
        let err = err.lock().unwrap().clone();
        let number = number.lock().unwrap().clone();
        let controll = controll.lock().unwrap().clone();
        //  let time_end = Instant::now();
        //  println!("Join pool: {:?}, wrap result {:?}, Tot : {:?}", time_poolend - time_forafter,
        //          time_end - time_poolend, time_end - time_forafter);

        (res,err,number,controll)
    }

    fn integrate_p3(&self, integral_method : &IntVec, fun : & FnVecP) -> (Vec<f64>, Vec<f64>, Vec<i32>, Vec<i32>){
        let dimension = fun.components.len();
        let res = Arc::new(Mutex::new(vec![0.0;dimension]));
        let err = Arc::new(Mutex::new(vec![0.0;dimension]));
        let number = Arc::new(Mutex::new(vec![0;dimension]));
        let controll = Arc::new(Mutex::new(vec![0;dimension]));
        let dimension = fun.components.len();
        let pool = ThreadPool::new(dimension);
        //  let start = Instant::now();
        for k in 0..dimension{
            //  let time1 = Instant::now();
            let pointer = fun.components[k].clone();
            let pointer2 = integral_method.components[0].clone();
            let res = res.clone();
            let err = err.clone();
            let number = number.clone();
            let controll = controll.clone();
            //  let time11 = Instant::now();
            let a = self.a;
            let b = self.b;
            let epsabs = self.epsabs;
            let epsrel = self.epsrel;
            //  let time12 = Instant::now();
            pool.execute( move || {
                //println!("thread creato");
                //  let time21 = Instant::now();
                let function = pointer.lock().unwrap();
                let integral = pointer2.lock().unwrap();
                //  let time22 = Instant::now();
                let results = integral.integrate(&*function,a,b,epsabs,epsrel);
                //  let time23 = Instant::now();
                res.lock().unwrap()[k] = results.0;
                err.lock().unwrap()[k] = results.1;
                number.lock().unwrap()[k] = results.2;
                controll.lock().unwrap()[k] = results.3;
                //println!("in thread : {:?},{:?},{:?}",time22 - time21,time23-time22,time23-time21);
            });
        }

        pool.join();

        //  let end = Instant::now();
        //  println!("{:?}",end-start);
        let res = res.lock().unwrap().clone();
        let err = err.lock().unwrap().clone();
        let number = number.lock().unwrap().clone();
        let controll = controll.lock().unwrap().clone();
        (res,err,number,controll)
    }

}
