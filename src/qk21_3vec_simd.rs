use std::iter::zip;
use std::simd::{Simd, SimdFloat};
use std::time::Instant;
use crate::funct_vector::FnVec3;
use crate::qk::*;


pub struct Qk21Vec3_Simd {}
///     Parameters:
///
///     On entry:
///         f   :   f64
///                 function
///
///         a   :   f64
///                 lower limit of integration
///
///         b   :   f64
///                 upper limit of integration
///
///     On return:
///         result  :   f64
///                     approximation to the integral i result is computed by applying
///                     the 61-point kronrod rule (resk) obtained by optimal addition
///                     of abscissae to the 30-point gauss rule(resg).
///
///         abserr  :   f64
///                     estimate of the modulus of the absolute error, which should not
///                     exceed abs(i-result)
///
///         resabs  :   f64
///                     approximation to the integral j
///
///         resasc  :   f64
///                     approximation to the integral of abs(f-i/(b-a)) over (a,b)
///
///     The abscissae and weights are given for the interval (-1,1).
///     Because of symmetry only the positive abscissae and their
///     corresponding weights are given.
///
///         xgk     :   abscissae of the 61-point kronrod rule
///                     xgk(2), xgk(4), ...  abscissae of the 30-point
///                     gauss rule
///                     xgk(1), xgk(3), ...  abscissae which are optimally
///                     added to the 30-point gauss rule
///
///         wgk     :   weights of the 61-point kronrod rule
///
///         wg      :   weights of the 30-point gauss rule
///
///
///     Gauss quadrature weights and kronrod quadrature abscissae and weights
///     as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
///     bell labs, nov. 1981.
///
///
///


const XGK: Simd<f64,64> = Simd::from_array([-0.995657163025808080735527280689003, -0.973906528517171720077964012084452,
    -0.930157491355708226001207180059508, -0.865063366688984510732096688423493,
    -0.780817726586416897063717578345042, -0.679409568299024406234327365114874,
    -0.562757134668604683339000099272694, -0.433395394129247190799265943165784,
    -0.294392862701460198131126603103866, -0.148874338981631210884826001129720,
    0.000000000000000000000000000000000, 0.148874338981631210884826001129720,
    0.294392862701460198131126603103866, 0.433395394129247190799265943165784,
    0.562757134668604683339000099272694, 0.679409568299024406234327365114874,
    0.780817726586416897063717578345042, 0.865063366688984510732096688423493,
    0.930157491355708226001207180059508, 0.973906528517171720077964012084452,
    0.995657163025808080735527280689003, -0.995657163025808080735527280689003,
    -0.973906528517171720077964012084452,
    -0.930157491355708226001207180059508, -0.865063366688984510732096688423493,
    -0.780817726586416897063717578345042, -0.679409568299024406234327365114874,
    -0.562757134668604683339000099272694, -0.433395394129247190799265943165784,
    -0.294392862701460198131126603103866, -0.148874338981631210884826001129720,
    0.000000000000000000000000000000000, 0.148874338981631210884826001129720,
    0.294392862701460198131126603103866, 0.433395394129247190799265943165784,
    0.562757134668604683339000099272694, 0.679409568299024406234327365114874,
    0.780817726586416897063717578345042, 0.865063366688984510732096688423493,
    0.930157491355708226001207180059508, 0.973906528517171720077964012084452,
    0.995657163025808080735527280689003, -0.995657163025808080735527280689003,
    -0.973906528517171720077964012084452,
    -0.930157491355708226001207180059508, -0.865063366688984510732096688423493,
    -0.780817726586416897063717578345042, -0.679409568299024406234327365114874,
    -0.562757134668604683339000099272694, -0.433395394129247190799265943165784,
    -0.294392862701460198131126603103866, -0.148874338981631210884826001129720,
    0.000000000000000000000000000000000, 0.148874338981631210884826001129720,
    0.294392862701460198131126603103866, 0.433395394129247190799265943165784,
    0.562757134668604683339000099272694, 0.679409568299024406234327365114874,
    0.780817726586416897063717578345042, 0.865063366688984510732096688423493,
    0.930157491355708226001207180059508, 0.973906528517171720077964012084452,
    0.995657163025808080735527280689003, 0.0]);


const WGK: Simd<f64,64> = Simd::from_array([0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
    0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
    0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
    0.123491976262065851077958109831074, 0.134709217311473325928054001771707,
    0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
    0.149445554002916905664936468389821, 0.147739104901338491374841515972068,
    0.142775938577060080797094273138717, 0.134709217311473325928054001771707,
    0.123491976262065851077958109831074, 0.109387158802297641899210590325805,
    0.093125454583697605535065465083366, 0.075039674810919952767043140916190,
    0.054755896574351996031381300244580, 0.032558162307964727478818972459390,
    0.011694638867371874278064396062192, 0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
    0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
    0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
    0.123491976262065851077958109831074, 0.134709217311473325928054001771707,
    0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
    0.149445554002916905664936468389821, 0.147739104901338491374841515972068,
    0.142775938577060080797094273138717, 0.134709217311473325928054001771707,
    0.123491976262065851077958109831074, 0.109387158802297641899210590325805,
    0.093125454583697605535065465083366, 0.075039674810919952767043140916190,
    0.054755896574351996031381300244580, 0.032558162307964727478818972459390,
    0.011694638867371874278064396062192, 0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
    0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
    0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
    0.123491976262065851077958109831074, 0.134709217311473325928054001771707,
    0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
    0.149445554002916905664936468389821, 0.147739104901338491374841515972068,
    0.142775938577060080797094273138717, 0.134709217311473325928054001771707,
    0.123491976262065851077958109831074, 0.109387158802297641899210590325805,
    0.093125454583697605535065465083366, 0.075039674810919952767043140916190,
    0.054755896574351996031381300244580, 0.032558162307964727478818972459390,
    0.011694638867371874278064396062192,0.0]);

const WG: Simd<f64,64> = Simd::from_array([0.0, 0.066671344308688137593568809893332, 0.0, 0.149451349150580593145776339657697,
    0.0, 0.219086362515982043995534934228163, 0.0, 0.269266719309996355091226921569469,
    0.0, 0.295524224714752870173892994651338, 0.0, 0.295524224714752870173892994651338,
    0.0, 0.269266719309996355091226921569469, 0.0, 0.219086362515982043995534934228163,
    0.0, 0.149451349150580593145776339657697, 0.0, 0.066671344308688137593568809893332,
    0.0, 0.0, 0.066671344308688137593568809893332, 0.0, 0.149451349150580593145776339657697,
    0.0, 0.219086362515982043995534934228163, 0.0, 0.269266719309996355091226921569469,
    0.0, 0.295524224714752870173892994651338, 0.0, 0.295524224714752870173892994651338,
    0.0, 0.269266719309996355091226921569469, 0.0, 0.219086362515982043995534934228163,
    0.0, 0.149451349150580593145776339657697, 0.0, 0.066671344308688137593568809893332,
    0.0, 0.0, 0.066671344308688137593568809893332, 0.0, 0.149451349150580593145776339657697,
    0.0, 0.219086362515982043995534934228163, 0.0, 0.269266719309996355091226921569469,
    0.0, 0.295524224714752870173892994651338, 0.0, 0.295524224714752870173892994651338,
    0.0, 0.269266719309996355091226921569469, 0.0, 0.219086362515982043995534934228163,
    0.0, 0.149451349150580593145776339657697, 0.0, 0.066671344308688137593568809893332,
    0.0, 0.0]);



impl Qk21Vec3_Simd {
    fn integrate(&self, f: &FnVec3, a: f64, b: f64, ) -> ([f64;3], [f64;3], [f64;3], [f64;3]) {
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let fv = fvec_simd(f, centr, hlgth);
        let resk = extract_scalar_product(&(fv * WGK).to_array());
        let reskh = [resk[0] * 0.5,resk[1] * 0.5,resk[2] * 0.5];
        let reskhs = Simd::from_array([reskh[0],reskh[0],reskh[0],reskh[0],
            reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],
            reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],reskh[0],
            reskh[1],reskh[1],reskh[1],reskh[1], reskh[1],reskh[1],reskh[1],reskh[1],reskh[1],
            reskh[1],reskh[1],reskh[1], reskh[1],reskh[1],reskh[1],reskh[1],reskh[1],reskh[1],
            reskh[1],reskh[1],reskh[1], reskh[2],reskh[2],reskh[2],reskh[2],
            reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],
            reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],reskh[2],0.0]);
        let mut resabs = extract_scalar_product(&(fv.abs() * WGK).to_array());
        let resg = extract_scalar_product(&(fv * WG).to_array());
        let mut resasc = extract_scalar_product(&(WGK * ( fv - reskhs).abs()).to_array());

        let result = [resk[0] * hlgth,resk[1] * hlgth,resk[2] * hlgth];
        resabs = [resabs[0] * dhlgth,resabs[1] * dhlgth,resabs[2] * dhlgth];
        resasc = [resasc[0] * dhlgth,resasc[1] * dhlgth,resasc[2] * dhlgth];
        let mut abserr = [((resk[0] - resg[0]) * hlgth).abs(),
            ((resk[1] - resg[1]) * hlgth).abs(),((resk[2] - resg[2]) * hlgth).abs()];

        for (i,j) in zip(&resasc,&mut abserr){
            if *i != 0.0 && *j != 0.0 {
                *j = *i * 1.0_f64.min((200.0 * *j / *i).powf(1.5));
            }
        }
        for (i,j) in zip(&resabs,&mut abserr){
            if *i > UFLOW / (50.0 * EPMACH) {
                *j = (*j).max((EPMACH * 50.0) * *i);
            }
        }




        //  if resasc != 0.0 && abserr != 0.0 {
        //      abserr = resasc * 1.0_f64.min((200.0 * abserr / resasc).powf(1.5));
        //  }
        //  if resabs > UFLOW / (50.0 * EPMACH) {
        //      abserr = abserr.max((EPMACH * 50.0) * resabs);
        //  }

        (result, abserr, resabs, resasc)
    }
}



pub fn extract_scalar_product(a : &[f64;64]) -> [f64;3]{
    let (mut sum1,mut sum2,mut sum3) = (0.0,0.0,0.0);
    for k in 0..21{
        sum1 += a[k];
        sum2 += a[k+21];
        sum3 += a[k+42];
    }
    [sum1,sum2,sum3]
}


pub fn fvec_simd(fun : &FnVec3, centr : f64, hlgth : f64) -> Simd<f64,64>{
    let f = &fun.components;
    Simd::from_array([f[0](centr + hlgth * XGK[0]), f[0](centr + hlgth * XGK[1]), f[0](centr + hlgth * XGK[2]),
        f[0](centr + hlgth * XGK[3]), f[0](centr + hlgth * XGK[4]), f[0](centr + hlgth * XGK[5]),
        f[0](centr + hlgth * XGK[6]), f[0](centr + hlgth * XGK[7]), f[0](centr + hlgth * XGK[8]),
        f[0](centr + hlgth * XGK[9]), f[0](centr + hlgth * XGK[10]),f[0](centr + hlgth * XGK[11]),
        f[0](centr + hlgth * XGK[12]), f[0](centr + hlgth * XGK[13]), f[0](centr + hlgth * XGK[14]),
        f[0](centr + hlgth * XGK[15]), f[0](centr + hlgth * XGK[16]), f[0](centr + hlgth * XGK[17]),
        f[0](centr + hlgth * XGK[18]), f[0](centr + hlgth * XGK[19]), f[0](centr + hlgth * XGK[20]),
        f[1](centr + hlgth * XGK[0]), f[1](centr + hlgth * XGK[1]), f[1](centr + hlgth * XGK[2]),
        f[1](centr + hlgth * XGK[3]), f[1](centr + hlgth * XGK[4]), f[1](centr + hlgth * XGK[5]),
        f[1](centr + hlgth * XGK[6]), f[1](centr + hlgth * XGK[7]), f[1](centr + hlgth * XGK[8]),
        f[1](centr + hlgth * XGK[9]), f[1](centr + hlgth * XGK[10]),f[1](centr + hlgth * XGK[11]),
        f[1](centr + hlgth * XGK[12]), f[1](centr + hlgth * XGK[13]), f[1](centr + hlgth * XGK[14]),
        f[1](centr + hlgth * XGK[15]), f[1](centr + hlgth * XGK[16]), f[1](centr + hlgth * XGK[17]),
        f[1](centr + hlgth * XGK[18]), f[1](centr + hlgth * XGK[19]), f[1](centr + hlgth * XGK[20]),
        f[2](centr + hlgth * XGK[0]), f[2](centr + hlgth * XGK[1]), f[2](centr + hlgth * XGK[2]),
        f[2](centr + hlgth * XGK[3]), f[2](centr + hlgth * XGK[4]), f[2](centr + hlgth * XGK[5]),
        f[2](centr + hlgth * XGK[6]), f[2](centr + hlgth * XGK[7]), f[2](centr + hlgth * XGK[8]),
        f[2](centr + hlgth * XGK[9]), f[2](centr + hlgth * XGK[10]),f[2](centr + hlgth * XGK[11]),
        f[2](centr + hlgth * XGK[12]), f[2](centr + hlgth * XGK[13]), f[2](centr + hlgth * XGK[14]),
        f[2](centr + hlgth * XGK[15]), f[2](centr + hlgth * XGK[16]), f[2](centr + hlgth * XGK[17]),
        f[2](centr + hlgth * XGK[18]), f[2](centr + hlgth * XGK[19]), f[2](centr + hlgth * XGK[20]), 0.0])
}


#[cfg(test)]
mod tests {
    use std::simd::Simd;
    use std::time::Instant;
    use crate::funct_vector::FnVec3;
    use crate::qk21::Qk21;
    use crate::qk21_1dvec_simd::Qk211DVec_Simd;
    use crate::qk21_3vec_simd::Qk21Vec3_Simd;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.sin();
        let a = 0.0;
        let b = 1.0;
        let qks = Qk21Vec3_Simd{};
        let qk = Qk211DVec_Simd{};
        let fun = FnVec3{ components : [Box::new(f),Box::new(f),Box::new(f)]};

        for k in 0..10000 {
            let start = Instant::now();
            let res2 = qk.integrate(&f, a, b);
            let res2 = qk.integrate(&f, a, b);
            let res2 = qk.integrate(&f, a, b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            let res1 = qks.integrate(&fun, a, b);
            println!("simd {:?}", start.elapsed());
        }
        //println!("{:?}",res1);
        //println!("{:?}",res1);


        let sim = Simd::from_array([1.0,2.0,3.0,4.0]);
        let arr = sim.to_array();
        let arr2 = &arr[1..3];
        let sum : f64 = arr[1..3].iter().sum();
        println!("{:?}",sim);
        println!("{:?}",arr2);
        println!("{:?}",sum);
    }
}


