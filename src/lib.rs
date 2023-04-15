#![feature(portable_simd)]
extern crate nalgebra as na;
extern crate core;

pub mod functions;
pub mod integral_method;
pub mod integrator;
pub mod funct_vector;
pub mod vector_integral_method;
pub mod vector_intergrator;
pub mod gauss_legendre;
pub mod qng;
pub mod quad_integrator;
pub mod qk15;
pub mod qk21;
pub mod qk31;
pub mod qk41;
pub mod qk51;
pub mod qk61;
pub mod qk;
pub mod qage;
pub mod qpsrt;
pub mod qag_integration_result;
pub mod qag_integrator_result;
pub mod qng_integration_result;
pub mod result_state;
pub mod qng_integrator_result;
pub mod quad_integral_method;
pub mod quad_integration_result;
pub mod quad_integrator_result;
pub mod special_function;
pub mod qagse;
pub mod qelg;
pub mod qagpe;
pub mod qk15_1dvec;
pub mod qk61_1dvec;
pub mod qk61_1dvec2;
pub mod qk61_1dvec3;
pub mod qk61_1dvec4;
pub mod qsrt2;
pub mod qage2;
pub mod qage_1dvec;
pub mod qag_1dvec_integration_result;
pub mod qag_1dvec_integrator_result;
pub mod qage_1dvec2;
pub mod qage_1dvec_parall;
pub mod qag_1dvec_parall_integration_result;
pub mod qag_1dvec_parall_integrator_result;
pub mod qage_vec;
mod qag_vec_integration_result;
mod qag_vec_integrator_result;
mod qk61_vec;
mod qk61_simd;
mod qk61_vec_simd_faster;
mod qage_4vec_simd;
mod qk61_4vec_simd;
mod qag_vec4_integration_result;
mod qag_vec4_integrator_result;
mod qk21_simd;
mod qk21_3vec_simd;
mod qage_1dvec_parall2;
mod qage_1dvec_parall3;
mod qage_1dvec_parall_2thread;
mod qage_1dvec_parall_4thread;
mod qage_1dvec_parall_8thread;
mod qk51_simd;
mod qk41_simd;
mod qk15_simd;
mod qk31_simd;


use functions::*;
use integral_method::*;
use integrator::*;
use funct_vector::*;
use vector_integral_method::*;
use vector_intergrator::*;
use std::time::Instant;
use std::sync::{Arc, Mutex};
use qng::*;
use quad_integrator::*;
use qage::*;
use special_function::*;
use ndarray::*;

#[cfg(test)]
mod tests {
    use core::time;
    use std::thread;
    use std::time::Duration;
    //  use std::io::IntoInnerError;
    use crate::gauss_legendre::{Gauss,GaussLegendre};
    use crate::qage2::Qag2;
    use crate::qage_1dvec2::Qag_1dvec2;
    use crate::qage_1dvec::Qag_1dvec;
    use crate::qage_1dvec_parall2::Qag_1dvec_parall2;
    use crate::qage_1dvec_parall3::Qag_1dvec_parall3;
    use crate::qage_1dvec_parall::Qag_1dvec_parall;
    use crate::qage_1dvec_parall_2thread::Qag_1dvec_parall_2thread;
    use crate::qage_1dvec_parall_4thread::Qag_1dvec_parall_4thread;
    use crate::qage_1dvec_parall_8thread::Qag_1dvec_parall_8thread;
    use crate::qagpe::Qagp;
    use crate::qagse::Qags;
    use crate::qk15::Qk15;
    use crate::qk::{EPMACH, Qk};
    use crate::qk15_1dvec::Qk151DVec;
    use crate::qk61::Qk61;
    use crate::qk61_1dvec2::Qk611DVec2;
    use crate::qk61_1dvec3::Qk611DVec3;
    use crate::qk61_1dvec4::Qk611DVec4;
    use crate::qk61_1dvec::Qk611DVec;
    use crate::quad_integral_method::QuadIntegralMethod;
    use super::*;

    #[test]
    fn vectorization(){
        let f = |x:f64| x.cos();
        let (a,b) = (0.0,1.0);
        let max = 200;
        let qk61 = Qk61{};
        let qk61_vec = Qk611DVec {};
        let qk61_vec2 = Qk611DVec2 {};
        let qk61_vec3 = Qk611DVec3 {};
        let qk61_vec4 = Qk611DVec4 {};
        let mut results = vec![vec![],vec![],vec![],vec![],vec![]];
        let mut tot = 0.0;
        let (mut t1,mut t2, mut t3,mut t4, mut t5) = ( 0.0, 0.0, 0.0,0.0, 0.0);

        for i in 1..max+1 {
            let time1 = Instant::now();
            results[1].push(qk61.integrate(&f, a, b));
            let time2 = Instant::now();
            //println!("{i} : {:?}",res1);
            if i > 100 {
                t1 += (time2 - time1).as_secs_f64();
                tot += 1.0;
            }
        }

        for i in 1..max+1 {
            let time1 = Instant::now();
            results[0].push(qk61_vec3.integrate(&f, a, b));
            let time2 = Instant::now();
            //println!("{i} : {:?}",res1);
            if i > 100 {
                t4 += (time2 - time1).as_secs_f64();
            }
        }


        for i in 1..max+1 {
            let time1 = Instant::now();
            results[2].push(qk61_vec.integrate(&f, a, b));
            let time2 = Instant::now();
            //println!("{i} : {:?}",res1);
            if i > 100 {
                t2 += (time2 - time1).as_secs_f64();
            }
        }
        for i in 1..max+1 {
            let time1 = Instant::now();
            results[3].push(qk61_vec2.integrate(&f, a, b));
            let time2 = Instant::now();
            //println!("{i} : {:?}",res1);
            if i > 100 {
                t3 += (time2 - time1).as_secs_f64();
            }
        }
        for i in 1..max+1 {
            let time1 = Instant::now();
            results[4].push(qk61_vec4.integrate(&f, a, b));
            let time2 = Instant::now();
            //println!("{i} : {:?}",res1);
            if i > 100 {
                t5 += (time2 - time1).as_secs_f64();
            }
        }

        let mut n = 0;
        for i in 1..max{
            if results[0][i] == results[1][i] && results[1][i] == results[2][i] && results[2][i] == results[3][i]
                && results[3][i] == results[4][i] {
                n += 1;
            }
        }

        t1 = t1 / tot;
        t2 = t2 / tot;
        t3 = t3 / tot;
        t4 = t4 / tot;
        t5 = t5 / tot;
        let equalres = n as f64 / tot;
        println!("time v1: {:?}; time v2: {:?}; time v3: {:?}; time v4 : {:?}; time v5 : {:?}; equalres : {equalres}"
                 , t1, t2,t3,t4,t5);
    }


    #[test]
    fn qags_test(){
        let test : bool = true; // set to true if you want to have the results printed out
        let sigma : f64 = 2000.0;
        let pi : f64 = std::f64::consts::PI;
        let norm = (2.0 * pi).sqrt();
        let f = |x:f64|  1.0 / (sigma * norm ) * (-0.5 * ( x / sigma).powi(2) ).exp();
        let f1 = |x:f64| (-x.powi(2)).exp();
        let f2 = |x:f64| 1.0 / ( 1.0 - x ).sqrt();
        let a = 0.0;
        let b = 1.0;
        let epsabs = 1.0;
        let epsrel = 0.0;
        let limit = 100;
        let key = 2;
        let npts2 = 3;
        let qagp = Qagp{limit, npts2, points : vec![0.9] };
        let qags = Qags{limit};
        let qag = Qag {key,limit};
        let time1 = Instant::now();
        let res = qagp.qintegrate(&f2, a, b, epsabs, epsrel).unwrap();
        let time2 = Instant::now();
        if test {
            println!("{:?}",res);
            println!("time:{:?}",time2-time1);
        }


    }

    #[test]
    fn qag_test(){
        let test : bool = true; // set to true if you want to have the results printed out
        let sigma : f64 = 2000.0;
        let pi : f64 = std::f64::consts::PI;
        let norm = (2.0 * pi).sqrt();
        let f = |x:f64|  1.0 / (sigma * norm ) * (-0.5 * ( x / sigma).powi(2) ).exp();
        let f1 = |x:f64| (-x.powi(2)).exp();
        let f2 = |x:f64| 10.0 + x.powi(2) - 10.0 * ( 2.0 * pi * x ).cos();
        let f2 = |x:f64| x.sin();
        let a = 0.0;
        let b = 1000.0;
        let epsabs = 1.0e-8;
        let epsrel = 0.0;
        let key = 2;
        let limit = 10000000;
        let max = 100;
        let qag = Qag {key,limit};
        let qag2 = Qag2{key,limit};
        let qag_vec = Qag_1dvec {key,limit};
        let qag_vec2 = Qag_1dvec2 {key,limit};
        let qag_parall_2thread = Qag_1dvec_parall_2thread {key,limit};
        let qag_parall_4thread = Qag_1dvec_parall_4thread {key,limit};
        let qag_parall_8thread = Qag_1dvec_parall_8thread {key,limit};
        let qag_vec_par2 = Qag_1dvec_parall2{key,limit};
        let mut vec : Vec<Arc<dyn Fn(f64)->f64 + Send + Sync>> = vec![];
        let (mut t1,mut t2,mut t3) = (0.0,0.0,0.0);
        let (mut t4, mut t5) = (0.0,0.0);



        vec.push(Arc::new(f2.clone()));
        let fun = FnVecPa{ components : vec};
        let fun2 = FnPa{ components : Arc::new(f2.clone())};
              for i in 0..max{
                  let start = Instant::now();
                  let res = qag_parall_8thread.qintegrate(&fun2, a, b, epsabs, epsrel).unwrap();
                  if i > 10 { t5 += start.elapsed().as_secs_f64();}
                  let start = Instant::now();
                  let res2 = qag_vec_par2.qintegrate(&fun,a,b,epsabs,epsrel).unwrap();
                  if i > 10 { t1 += start.elapsed().as_secs_f64();}
                  let start = Instant::now();
                  let res3 = qag_parall_4thread.qintegrate(&fun2,a,b,epsabs,epsrel).unwrap();
                  if i > 10 { t4 += start.elapsed().as_secs_f64();}
                  let start = Instant::now();
                  let res4 = qag_vec2.qintegrate(&f2,a,b,epsabs,epsrel).unwrap();
                  if i > 10 { t2 += start.elapsed().as_secs_f64();}
                  let start = Instant::now();
                  let res5 = qag_parall_2thread.qintegrate(&fun2, a, b, epsabs, epsrel).unwrap();
                  if i > 10 { t3 += start.elapsed().as_secs_f64();}
                  if test && i == max-1 {
                      println!("{:?}",res4.neval);
                      println!("{:?}",res2.neval);
                      //println!("{:?}",res4.result);
                      //println!("{:?}",res4.abserr);
                      println!("{:?}",res5.neval);
                  }
              }
        t1 = t1 / ( max as f64 - 10.0);
        t2 = t2 / ( max as f64 - 10.0);
        t3 = t3 / ( max as f64 - 10.0);
        t4 = t4 / ( max as f64 - 10.0);
        t5 = t5 / ( max as f64 - 10.0);


        println!("best vec : {t2}; 2thread parall : {t3}; rayon parall : {t1}; 4thread parall : {t4};\
        8thread parall : {t5}");


    }


    #[test]
    fn qagg_test(){
        let test : bool = false; // set to true if you want to have the results printed out


        let (fun,fun2) = random_polynomials2(10);

        let qng = Qng{};

        let s_int = SIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : Box::new(qng.clone()),
        };


        let qng2 = IntVec{
            components : vec![Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone()))],
        };

        let qng3 = IntVec{
            components : vec![Arc::new(Mutex::new(qng.clone()))],
        };

        let p_int = PIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : qng2,
        };

        let p2_int = PIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : qng3,
        };
        let (mut p1,mut p2,mut p3) = (0,0,0);
        let (mut t1, mut t2, mut t3) = (0.0,0.0,0.0);
        let mut i = 0;
        let total = 1000;
        for _k in 0..total{
            let time1 = Instant::now();
            let res1 = s_int.integrate( &fun);
            let time2 = Instant::now();
            let res2 = p_int.integrate_rayon(&fun2);
            let time3 = Instant::now();
            let res3 = p_int.integrate_handles(&fun2);
            let time4 = Instant::now();
            let res4 = p_int.integrate_pool(&fun2);
            let time5 = Instant::now();
            if time3-time2 < time2-time1 {
                p1 += 1;
                t1 += (time3-time2).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time4-time3 < time2-time1 {
                p2 += 1;
                t2 += (time4-time3).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time5-time4 < time2-time1 {
                p3 += 1;
                t3 += (time5-time4).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time3-time2 < time2-time1 || time4-time3 < time2-time1 || time5-time4 < time2-time1 {
                i += 1;
            }
            if test{
                println!("{:?},{:?},{:?},{:?}",time2-time1,time3-time2,time4-time3,time5-time4 );
            }
        }
        println!("out of a total of {total} integration {i} time one of the parallel one was faster.\n \
        rayon : {p1}, on avg. {} of time\n handles: {p2}, on avg. {} of time \n \
        threadpool: {p3}, on avg. {} of time",t1/p1 as f64,t2/p2 as f64,t3/p3 as f64);

        let qag = Qag{key:6,limit:10};

        let s_int = SIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : Box::new(qag.clone()),
        };


        let qng2 = IntVec{
            components : vec![Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone()))],
        };

        let qng3 = IntVec{
            components : vec![Arc::new(Mutex::new(qag.clone()))],
        };

        let p_int = PIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : qng2,
        };

        let p2_int = PIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : qng3,
        };
        (p1,p2,p3) = (0,0,0);
        (t1,t2,t3) = (0.0,0.0,0.0);
        let mut i = 0;
        let total = 1000;
        for _k in 0..total{
            let time1 = Instant::now();
            let res1 = s_int.integrate( &fun);
            let time2 = Instant::now();
            let res2 = p_int.integrate_rayon(&fun2);
            let time3 = Instant::now();
            let res3 = p_int.integrate_handles(&fun2);
            let time4 = Instant::now();
            let res4 = p_int.integrate_pool(&fun2);
            let time5 = Instant::now();
            if time3-time2 < time2-time1 {
                p1 += 1;
                t1 += (time3-time2).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time4-time3 < time2-time1 {
                p2 += 1;
                t2 += (time4-time3).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time5-time4 < time2-time1 {
                p3 += 1;
                t3 += (time5-time4).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time3-time2 < time2-time1 || time4-time3 < time2-time1 || time5-time4 < time2-time1 {
                i += 1;
            }
        }
        println!("out of a total of {total} integration {i} time one of the parallel one was faster.\n \
        rayon : {p1}, on avg. {} of time\n handles: {p2}, on avg. {} of time \n \
        threadpool: {p3}, on avg. {} of time",t1/p1 as f64,t2/p2 as f64,t3/p3 as f64);





    }

    #[test]
    fn qng_test(){
        let test : bool = false; // set to true if you want to have the results printed out
        let f1 = |x:f64|  (x*x.sin().abs()).powf(1.5);
        let f2 = |x:f64|  (x.powi(3)*x.cos().abs()).powf(2.5);
        let f3 = |x:f64|  rastrigin(x);
        /*
        let f4 = |x:f64|  (x.powi(7)*x.sin()).powf(3.5);
        let f5 = |x:f64|  (x+5.0).ln() * x.sin().powf(3.3);
        let f6 = |x:f64|  x.exp();
        let f7 = |x:f64|  1.0 / (10.0 * (2.0 * std::f64::consts::PI).sqrt() ) * (-0.5 * ( x / 10.0).powi(2) ).exp();
        let f8 = |x:f64|  1.0 / (50.0 * (2.0 * std::f64::consts::PI).sqrt() ) * (-0.5 * ( x / 50.0).powi(2) ).exp();
        let f9 = |x:f64|  1.0 / (200.0 * (2.0 * std::f64::consts::PI).sqrt() ) * (-0.5 * ( x / 200.0).powi(2) ).exp();


         */
        let parameters4 = random_vector();
        let parameters4p = parameters4.clone();
        let parameters5 = random_vector();
        let parameters5p = parameters5.clone();
        let parameters6 = random_vector();
        let parameters6p = parameters6.clone();
        let parameters7 = random_vector();
        let parameters7p = parameters7.clone();
        let parameters8 = random_vector();
        let parameters8p = parameters8.clone();
        let parameters9 = random_vector();
        let parameters9p = parameters9.clone();


        let f4 = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters4.len(){
                tot += parameters4[i] * x.powi(i as i32);
            }
            tot
        };
        let f4p = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters4p.len(){
                tot += parameters4p[i] * x.powi(i as i32);
            }
            tot
        };
        let f5 = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters5.len(){
                tot += parameters5[i] * x.powi(i as i32);
            }
            tot
        };
        let f5p = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters5p.len(){
                tot += parameters5p[i] * x.powi(i as i32);
            }
            tot
        };
        let f6 = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters6.len(){
                tot += parameters6[i] * x.powi(i as i32);
            }
            tot
        };
        let f6p = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters6p.len(){
                tot += parameters6p[i] * x.powi(i as i32);
            }
            tot
        };
        let f7 = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters7.len(){
                tot += parameters7[i] * x.powi(i as i32);
            }
            tot
        };
        let f7p = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters7p.len(){
                tot += parameters7p[i] * x.powi(i as i32);
            }
            tot
        };
        let f8 = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters8.len(){
                tot += parameters8[i] * x.powi(i as i32);
            }
            tot
        };
        let f8p = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters8p.len(){
                tot += parameters8p[i] * x.powi(i as i32);
            }
            tot
        };
        let f9 = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters9.len(){
                tot += parameters9[i] * x.powi(i as i32);
            }
            tot
        };
        let f9p = move |x:f64| {
            let mut tot = 0.0;
            for i in 0..parameters9p.len(){
                tot += parameters9p[i] * x.powi(i as i32);
            }
            tot
        };




        let fun = FnVec  {
            components : vec![Box::new(f1),Box::new(f2),Box::new(f3),Box::new(f4),
                              Box::new(f5),Box::new(f6),Box::new(f7),Box::new(f8),
                              Box::new(f9)],
        };
        let fun2 = FnVecP{
            components : vec![Arc::new(Mutex::new(f1)),Arc::new(Mutex::new(f2)),
                              Arc::new(Mutex::new(f3)),Arc::new(Mutex::new(f4p)),
                              Arc::new(Mutex::new(f5p)),Arc::new(Mutex::new(f6p)),
                              Arc::new(Mutex::new(f7p)),Arc::new(Mutex::new(f8p)),
                              Arc::new(Mutex::new(f9p))],
        };


        let qng = Qng{};

        let s_int = SIntegrator {
            a : 0.0,
            b : 1000.0,
            epsabs : 0.0,
            epsrel : 2.0 * 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : Box::new(qng.clone()),
        };


        let qng2 = IntVec{
            components : vec![Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone())),Arc::new(Mutex::new(qng.clone())),
                              Arc::new(Mutex::new(qng.clone()))],
        };

        let qng3 = IntVec{
            components : vec![Arc::new(Mutex::new(qng.clone()))],
        };

        let p_int = PIntegrator {
            a : 0.0,
            b : 1000.0,
            epsabs : 0.0,
            epsrel : 2.0 * 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : qng2,
        };

        let p2_int = PIntegrator {
            a : 0.0,
            b : 1000.0,
            epsabs : 0.0,
            epsrel : 2.0 * 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : qng3,
        };
        let (mut p1,mut p2,mut p3) = (0,0,0);
        let (mut t1, mut t2, mut t3) = (0.0,0.0,0.0);
        let mut i = 0;
        let total = 1000;
        for _k in 0..total{
            let time1 = Instant::now();
            let res1 = s_int.integrate( &fun);
            let time2 = Instant::now();
            let res2 = p_int.integrate_rayon(&fun2);
            let time3 = Instant::now();
            let res3 = p_int.integrate_handles(&fun2);
            let time4 = Instant::now();
            let res4 = p_int.integrate_pool(&fun2);
            let time5 = Instant::now();
            if time3-time2 < time2-time1 {
                p1 += 1;
                t1 += (time3-time2).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time4-time3 < time2-time1 {
                p2 += 1;
                t2 += (time4-time3).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time5-time4 < time2-time1 {
                p3 += 1;
                t3 += (time5-time4).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time3-time2 < time2-time1 || time4-time3 < time2-time1 || time5-time4 < time2-time1 {
                i += 1;
            }
            if test{
                println!("{:?},{:?},{:?},{:?}",time2-time1,time3-time2,time4-time3,time5-time4 );
            }
        }
        println!("out of a total of {total} integration {i} time one of the parallel one was faster.\n \
        rayon : {p1}, on avg. {} of time\n handles: {p2}, on avg. {} of time \n \
        threadpool: {p3}, on avg. {} of time",t1/p1 as f64,t2/p2 as f64,t3/p3 as f64);



        let qag = Qag{key:6,limit:10};

        let s_int = SIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : Box::new(qag.clone()),
        };


        let qng2 = IntVec{
            components : vec![Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone())), Arc::new(Mutex::new(qag.clone())),
                              Arc::new(Mutex::new(qag.clone()))],
        };

        let qng3 = IntVec{
            components : vec![Arc::new(Mutex::new(qag.clone()))],
        };

        let p_int = PIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : qng2,
        };

        let p2_int = PIntegrator {
            a : 0.0,
            b : 1.0,
            epsabs : 0.0,
            epsrel : 0.5e-28_f64.max(50.0 * EPMACH),
            integral_method : qng3,
        };
        (p1,p2,p3) = (0,0,0);
        (t1,t2,t3) = (0.0,0.0,0.0);
        let mut i = 0;
        let total = 1000;
        for _k in 0..total{
            let time1 = Instant::now();
            let res1 = s_int.integrate( &fun);
            let time2 = Instant::now();
            let res2 = p_int.integrate_rayon(&fun2);
            let time3 = Instant::now();
            let res3 = p_int.integrate_handles(&fun2);
            let time4 = Instant::now();
            let res4 = p_int.integrate_pool(&fun2);
            let time5 = Instant::now();
            if time3-time2 < time2-time1 {
                p1 += 1;
                t1 += (time3-time2).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time4-time3 < time2-time1 {
                p2 += 1;
                t2 += (time4-time3).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time5-time4 < time2-time1 {
                p3 += 1;
                t3 += (time5-time4).as_secs_f64()/(time2-time1).as_secs_f64();
            }
            if time3-time2 < time2-time1 || time4-time3 < time2-time1 || time5-time4 < time2-time1 {
                i += 1;
            }
        }
        println!("out of a total of {total} integration {i} time one of the parallel one was faster.\n \
        rayon : {p1}, on avg. {} of time\n handles: {p2}, on avg. {} of time \n \
        threadpool: {p3}, on avg. {} of time",t1/p1 as f64,t2/p2 as f64,t3/p3 as f64);



    }

    #[test]
    fn gauss_legendre(){
        let test : bool = false; // set to true if you want to have the results printed out
        let precision : f64 = 0.001 ; // set the relative errors required
        let parameters = random_vector();
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        //let polynomial = Sin{i : 0} ;
        let simpson2 = Simpson2{};
        let integral_precision = FixedPrecision{precision};
        let exact_integral = exact_polynom_integral(&parameters);
        //let exact_integral = - 1.0_f64.cos() + 1.0 ;

        let gauss = GaussLegendre::new( 8);
        if test{
            println!("{}\n{:?}\n{:?}\n{:?}", parameters.len(), integral_precision.integrate(&simpson2,&polynomial),
                     exact_integral, gauss.integrate(&polynomial));
        }


    }



    //-------------------   RELATIVE ERRORS   -----------------

    #[test]
    fn parallel() {
        let test : bool = false; // set to true if you want to have the results printed out
        let integration_method = VecSimpson2{};
        let number_of_points = vec![10000;9];
        let number_of_points2 = vec![10000;9];

        let parameters1 = [random_vector(),random_vector(),random_vector()];
        let parameters2 = [random_vector(),random_vector(),random_vector()];
        let parameters3 = [random_vector(),random_vector(),random_vector()];
        let (polynomial1, polynomial2,polynomial3) =
            (PolynomialFunction::new(parameters1[0].to_vec()),
             PolynomialFunction::new(parameters1[1].to_vec()),
             PolynomialFunction::new(parameters1[2].to_vec()) );
        let (polynomial4, polynomial5,polynomial6) =
            (PolynomialFunction::new(parameters2[0].to_vec()),
             PolynomialFunction::new(parameters2[1].to_vec()),
             PolynomialFunction::new(parameters2[2].to_vec()) );
        let (polynomial7, polynomial8,polynomial9) =
            (PolynomialFunction::new(parameters3[0].to_vec()),
             PolynomialFunction::new(parameters3[1].to_vec()),
             PolynomialFunction::new(parameters3[2].to_vec()) );

        let collection : FunVector = FunVector{ components : vec![Arc::new(Mutex::new(polynomial1.clone())),
                                                                  Arc::new(Mutex::new(polynomial2.clone())),
                                                                  Arc::new(Mutex::new(polynomial3.clone())),
                                                                  Arc::new(Mutex::new(polynomial4.clone())),
                                                                  Arc::new(Mutex::new(polynomial5.clone())),
                                                                  Arc::new(Mutex::new(polynomial6.clone())),
                                                                  Arc::new(Mutex::new(polynomial7.clone())),
                                                                  Arc::new(Mutex::new(polynomial8.clone())),
                                                                  Arc::new(Mutex::new(polynomial9.clone()))]};

        let collection2 : FunctVector = FunctVector{ components : vec![Box::new(polynomial1),Box::new(polynomial2),Box::new(polynomial3),
                                                                       Box::new(polynomial4),Box::new(polynomial5),Box::new(polynomial6),
                                                                       Box::new(polynomial7),Box::new(polynomial8),Box::new(polynomial9)]};

        let time3 = Instant::now();
        let res_parall = integration_method.integrate_non_uniform_parallel(&collection, number_of_points);
        let time4 = Instant::now();
        let res = integration_method.integrate_non_uniform(&collection2, number_of_points2);
        let time5 = Instant::now();
        let parallel_time = time4 - time3;
        let non_parallel_time = time5 - time4;
        let proportion = parallel_time.as_secs_f64()/non_parallel_time.as_secs_f64();
        //assert!(res == res_parall && parallel_time < non_parallel_time );

        if test {
            println!("{:?},{:?}", res, res_parall);
            println!("Parallel function duration :{:?},\nNon parallel function duration :{:?}\n{:?}", parallel_time, non_parallel_time,proportion);
        }

    }


    #[test]
    fn error_comparison() {
        let test : bool = false; // set to true if you want to have the results printed out
        let precision : f64 = 0.5 ; // set the relative errors required

        let parameters = random_vector();
        let polynomial = PolynomialFunction::new(parameters.to_vec());
        //let polynomial = Sin{i : 0} ;
        let integral_precision = FixedPrecision{precision};
        let rectangular = Rectangular{};
        let trapezoidal = Trapezoiadal{};
        let simpson1 = Simpson1{};
        let simpson2 = Simpson2{};
        let exact_integral = exact_polynom_integral(&parameters) ;
        //let exact_integral = - 1.0_f64.cos() + 1.0 ;


        let rel_error_real_rectangular = ((integral_precision.integrate(&rectangular,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_trapezoidal = ((integral_precision.integrate(&trapezoidal,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_simpson1 = ((integral_precision.integrate(&simpson1,&polynomial) - exact_integral).abs() / exact_integral).abs() ;
        let rel_error_real_simpson2 = ((integral_precision.integrate(&simpson2,&polynomial) - exact_integral).abs() / exact_integral).abs() ;

        if test {
            println!("Exact integral : {}", exact_integral);
            println!("The random polynomial function has order : {}", parameters.len());
            println!("Rectangular: {} , with {} points ( relative errors bound : {}; real errors : {})",
                     integral_precision.integrate(&rectangular, &polynomial),
                     integral_precision.number_of_points(&rectangular, &polynomial),
                     integral_precision.relative_error(&rectangular, &polynomial),
                     rel_error_real_rectangular);
            println!("Trapezoidal: {}, with {} points ( relative errors bound : {}; real errors : {})",
                     integral_precision.integrate(&trapezoidal, &polynomial),
                     integral_precision.number_of_points(&trapezoidal, &polynomial),
                     integral_precision.relative_error(&trapezoidal, &polynomial),
                     rel_error_real_trapezoidal);
            println!("Simpson1: {}, with {} points ( relative errors bound : {}; real errors : {})",
                     integral_precision.integrate(&simpson1, &polynomial),
                     integral_precision.number_of_points(&simpson1, &polynomial),
                     integral_precision.relative_error(&simpson1, &polynomial),
                     rel_error_real_simpson1);
            println!("Simpson2: {}, with {} points ( relative errors bound : {}; real errors : {})",
                     integral_precision.integrate(&simpson2, &polynomial),
                     integral_precision.number_of_points(&simpson2, &polynomial),
                     integral_precision.relative_error(&simpson2, &polynomial),
                     rel_error_real_simpson2);
        }

        assert!(rel_error_real_rectangular < precision && rel_error_real_trapezoidal < precision
                && rel_error_real_simpson1 < precision && rel_error_real_simpson2 < precision)

    }

//----------------------- TEST OF FUNCTVECTORS


    #[test]
    fn funct_vector() {
        let test : bool = false; // set to true if you want to have the results printed out
        let precision :f64 = 0.5; // set the relative errors required

        let parameters = [random_vector(),random_vector(),random_vector()];
        let functions_order = [ parameters[0].len(),parameters[1].len(),parameters[2].len()];
        let (polynomial1, polynomial2,polynomial3) =
            (PolynomialFunction::new(parameters[0].to_vec()),
            PolynomialFunction::new(parameters[1].to_vec()),
            PolynomialFunction::new(parameters[2].to_vec()) );
        let collection : FunctVector = FunctVector{ components : vec![Box::new(polynomial1),Box::new(polynomial2),Box::new(polynomial3)]};
        let rectangular = VecRectangular{};
        let trapezoidal = VecTrapezoiadal{};
        let simpson1 = VecSimpson1{};
        let simpson2 = VecSimpson2{};
        let vec_integral_method = VecFixedPrecision{precision};
        let exact_integral = [exact_polynom_integral(&parameters[0]),
            exact_polynom_integral(&parameters[1]),exact_polynom_integral(&parameters[2])];

        // UNIFORM
        let result = [vec_integral_method.integrate_uniform(&rectangular,&collection),
            vec_integral_method.integrate_uniform(&trapezoidal,&collection),
            vec_integral_method.integrate_uniform(&simpson1,&collection),
            vec_integral_method.integrate_uniform(&simpson2,&collection)];
        let number_of_points = [vec_integral_method.number_of_points_uniform(&rectangular,&collection),
            vec_integral_method.number_of_points_uniform(&trapezoidal,&collection),
            vec_integral_method.number_of_points_uniform(&simpson1,&collection),
            vec_integral_method.number_of_points_uniform(&simpson2,&collection)];
        let relative_errors = [vec_integral_method.relative_errors_uniform(&rectangular, &collection),
            vec_integral_method.relative_errors_uniform(&trapezoidal, &collection),
            vec_integral_method.relative_errors_uniform(&simpson1, &collection),
            vec_integral_method.relative_errors_uniform(&simpson2, &collection)];
        let rel_error_real = [[((result[0][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
            ((result[0][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
            ((result[0][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result[1][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result[1][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result[1][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result[2][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result[2][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result[2][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result[3][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result[3][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result[3][2] - exact_integral[2]).abs() / exact_integral[2]).abs()]];


        // NON UNIFORM (nu)
        let result_nu = [vec_integral_method.integrate_non_uniform(&rectangular,&collection),
            vec_integral_method.integrate_non_uniform(&trapezoidal,&collection),
            vec_integral_method.integrate_non_uniform(&simpson1,&collection),
            vec_integral_method.integrate_non_uniform(&simpson2,&collection)];
        let number_of_points_nu = [vec_integral_method.number_of_points_non_uniform(&rectangular,&collection),
            vec_integral_method.number_of_points_non_uniform(&trapezoidal,&collection),
            vec_integral_method.number_of_points_non_uniform(&simpson1,&collection),
            vec_integral_method.number_of_points_non_uniform(&simpson2,&collection)];
        let relative_errors_nu = [vec_integral_method.relative_errors_non_uniform(&rectangular, &collection),
            vec_integral_method.relative_errors_non_uniform(&trapezoidal, &collection),
            vec_integral_method.relative_errors_non_uniform(&simpson1, &collection),
            vec_integral_method.relative_errors_non_uniform(&simpson2, &collection)];
        let rel_error_real_nu = [[((result_nu[0][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
            ((result_nu[0][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
            ((result_nu[0][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result_nu[1][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result_nu[1][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result_nu[1][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result_nu[2][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result_nu[2][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result_nu[2][2] - exact_integral[2]).abs() / exact_integral[2]).abs()],
            [((result_nu[3][0] - exact_integral[0]).abs() / exact_integral[0]).abs(),
                ((result_nu[3][1] - exact_integral[1]).abs() / exact_integral[1]).abs(),
                ((result_nu[3][2] - exact_integral[2]).abs() / exact_integral[2]).abs()]];

        if test {
            println!("UNIFORM:\nThe random polynomial functions have order : {:?}", functions_order);
            println!("Notation : [Rectangular, Trapezoidal, Simpson1, Simpson2]\nNumber of point used : {:?} ", number_of_points);
            println!("Integral result: [first function, second function, third function]\n{:?}", result);
            //println!("{:?}", rectangular.integrate_uniform_parallel(&collection,number_of_points[0]));
            println!("Relative errors bound:\n{:?}", relative_errors);
            println!("NON-UNIFORM:\nThe random polynomial functions have order : {:?}", functions_order);
            println!("Notation : [Rectangular, Trapezoidal, Simpson1, Simpson2]\nNumber of point used : {:?} ", number_of_points_nu);
            println!("Integral result: [first function, second function, third function]\n{:?}", result_nu);
            println!("Relative errors bound:\n{:?}", relative_errors_nu);
        }

        assert!(rel_error_real[0][0] < precision && rel_error_real[0][1] < precision &&
            rel_error_real[0][2] < precision &&
            rel_error_real[1][0] < precision && rel_error_real[1][1] < precision &&
            rel_error_real[1][2] < precision &&
            rel_error_real[2][0] < precision && rel_error_real[2][1] < precision &&
            rel_error_real[2][2] < precision &&
            rel_error_real[3][0] < precision && rel_error_real[3][1] < precision &&
            rel_error_real[3][2] < precision &&
            rel_error_real_nu[0][0] < precision && rel_error_real_nu[0][1] < precision &&
            rel_error_real_nu[0][2] < precision &&
            rel_error_real_nu[1][0] < precision && rel_error_real_nu[1][1] < precision &&
            rel_error_real_nu[1][2] < precision &&
            rel_error_real_nu[2][0] < precision && rel_error_real_nu[2][1] < precision &&
            rel_error_real_nu[2][2] < precision &&
            rel_error_real_nu[3][0] < precision && rel_error_real_nu[3][1] < precision &&
            rel_error_real_nu[3][2] < precision)
        }





}




