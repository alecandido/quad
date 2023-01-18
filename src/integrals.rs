use crate::functions::*;

pub fn rectangular_integral(funct : &impl Function, number_of_point : i32 ) -> f32 {
    let mut tot : f32 = funct.evaluate(&0.0) ;
    let m = number_of_point as f32 ;
    for i in 1..number_of_point {
        let j =  i as f32 ;
        tot += funct.evaluate( &(j/m) ) ;
    }
    tot / m
}

pub fn trapezoidal_integral(funct : &impl Function, number_of_point : i32 ) -> f32 {
    let mut tot : f32 = 0.5 * (funct.evaluate(&0.0) + funct.evaluate(&1.0)) ;
    let m = number_of_point as f32 ;
    for i in 1..number_of_point {
        let j =  i as f32 ;
        tot += funct.evaluate( &(j/m) ) ;
    }
    tot / m
}

pub fn simpson1(funct : &impl Function, number_of_point : i32 ) -> f32 { // number_of_point has to be even !!!
    let mut tot : f32 =  funct.evaluate(&0.0) - funct.evaluate(&1.0) ;
    let m = number_of_point as f32 ;
    for i in 1..(number_of_point/2 + 1 ) {
        let j =  i as f32 ;
        tot += 4.0 * funct.evaluate( &((2.0 * j - 1.0 )/m) ) + 2.0 * funct.evaluate(&((2.0 * j)/m) ) ;
    }
    tot / ( 3.0 * m )
}


pub fn simpson2(funct : &impl Function, number_of_point : i32 ) -> f32 { // number_of_point has to be even !!!
    let mut tot : f32 = funct.evaluate(&0.0) + funct.evaluate(&1.0) ;
    let m = number_of_point as f32 ;
    for i in 1..(number_of_point ) {
        let j =  i as f32 ;
        if i % 3 == 0 {
            tot += 2.0 * funct.evaluate(&(j/m)) ;
        }
        else {
            tot += 3.0 * funct.evaluate(&(j/m)) ;
        }
    }
    3.0 * tot / ( 8.0 * m )
}