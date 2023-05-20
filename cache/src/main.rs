#![feature(portable_simd)]
extern crate nalgebra as na;

use std::simd::Simd;
use std::time::Instant;
use crate::na::SimdComplexField;
use simba::simd::f64x4;
use na::{matrix, Matrix3x1, SMatrix, SVector, vector};

fn main(){
    let m1: Matrix3x1<f64x4> = na::matrix![
                                        f64x4::new(1.0, 2.0, 3.0, 4.0);
                                        f64x4::new(5.0, 6.0, 7.0, 8.0);
                                        f64x4::new(9.0, 10.0, 11.0, 12.0);
                                   ];
    let mut m2: Matrix3x1<f64x4> = na::matrix![
                                        f64x4::new(1.0, 2.0, 3.0, 4.0);
                                        f64x4::new(5.0, 6.0, 7.0, 8.0);
                                        f64x4::new(9.0, 10.0, 11.0, 12.0);
                                    ];
    let m3 = m1+m2;
    println!("m3 is: {:?}", m3);

    let m4 = m2.transpose();
    println!("m2 shape is: {:?}", m2.shape());
    println!("m4 shape is: {:?}", m4.shape());
    println!("m1 shape is: {:?}", m1.shape());
    let m5 = m4*m1;
    println!("m5 is: {:?}", m5);


    let x = packed_simd_2::f64x4::new(1.0,1.0,1.0,1.0);
    let x = x * 2.0;
    println!("x: {:?}",x);


    let v1: SVector<f64x4,2> = na::vector![f64x4::new(1.0, 2.0, 3.0, 4.0),f64x4::new(1.0, 0.0, 3.0, 4.0)];
    let v3 = na::vector![1.0,2.0,3.0,4.0];
    let v4 = na::vector![1.0,2.0,3.0,4.0];
    let a = v1.x.simd_horizontal_sum();
    let t = f64x4::new(1.0, 2.0, 3.0, 4.0);
    let t = t.0.cos();
    println!{"{:?}",t};


    println!("a : {:?}",a);
    let v2: SVector<f64x4,2> = na::vector![f64x4::new(1.0, 3.0, 3.0, 4.0),f64x4::new(1.0, 2.0, 3.0, 4.0)];
    let res = (v1.dot(&v2));
    println!("{:?}",res.simd_horizontal_sum());




    let f1 = |x:f64| x.cos();
    let f2 = |x:f64| x.sin();
    let f3 = |x:f64| x.sin() + x.cos();
    let f4 = |x:f64| x.sin() - x.cos();

    let f = |x:f64| Simd::from_array([x.cos(),x.sin(),x.sin() + x.cos(),x.sin() - x.cos()]);

    let mut res ;
    let mut res2;

    for _k in 0..100{
        let start = Instant::now();
        for x in 0..100{
            res = Simd::from_array([f1(x as f64),f2(x as f64),f3(x as f64),f4(x as f64)]);
        }
        println!("normal : {:?}", start.elapsed());
        let start = Instant::now();
        for x in 0..100{
            res2 = f(x as f64);
        }
        println!("array : {:?}", start.elapsed());
    }



    let x = [1.0,2.0,3.0,4.0];
    let start = Instant::now();
    let res = [f1(x[0]),f1(x[1]),f1(x[2]),f1(x[3])];
    println!("norm : {:?}",start.elapsed());
    println!("{:?}",res);
    let start = Instant::now();
    let res = x.map(f1);
    println!("map : {:?}",start.elapsed());
    println!("{:?}",res);






}
