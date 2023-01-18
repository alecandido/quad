pub mod functions;
pub mod integrals;

use functions::*;
use integrals::*;


fn main() {
    let line = Line{ a : 1.1 , b : 2.2 } ;
    let parabola = Parabola{ a : 3.4 , b : 6.5 , c : 7.8 } ;
    let cubic = Cubic{ a : 7.6 , b : 1.9 , c : 9.1 , d : 11.2 } ;
    let sin = Sin{} ;
    // println!("{}", line1.evaluate(&x)) ;
    for i in 1..10 {
        println!("{i} : Line  R =  {}, T = {} , S1 = {} , S2 = {}", rectangular_integral(&line, i * 10 ) ,
                 trapezoidal_integral(&line, i * 10 ), simpson1(&line, i * 10),
                simpson2(&line, i * 10));
        println!("Parabola  R =  {}, T = {} , S1 = {} , S2 = {}", rectangular_integral(&parabola, i * 10 ) ,
                 trapezoidal_integral(&parabola, i * 10 ), simpson1(&parabola, i * 10 ),
                simpson2(&parabola, i*10 ));
        println!("Cubic  R =  {}, T = {} , S1 = {} , S2 = {}", rectangular_integral(&cubic, i * 10 ),
                 trapezoidal_integral(&cubic, i * 10 ), simpson1(&cubic, i * 10 ),
                simpson2(&cubic, i*10 ) );
        println!("Sin  R =  {}, T = {} , S1 = {} , S2 = {}", rectangular_integral(&sin, i * 10 ),
                 trapezoidal_integral(&sin, i * 10 ), simpson1(&sin, i * 10 ),
                 simpson2(&sin, i*10 ) );
    }
}

