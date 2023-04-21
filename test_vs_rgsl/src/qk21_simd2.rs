use std::iter::zip;
use std::simd::{Simd, SimdFloat};
use std::time::Instant;
use crate::qk::*;


pub struct Qk21Simd2 {}
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


pub const XGK: [f64;21] = [-0.995657163025808080735527280689003, -0.973906528517171720077964012084452,
    -0.930157491355708226001207180059508, -0.865063366688984510732096688423493,
    -0.780817726586416897063717578345042, -0.679409568299024406234327365114874,
    -0.562757134668604683339000099272694, -0.433395394129247190799265943165784,
    -0.294392862701460198131126603103866, -0.148874338981631210884826001129720,
    0.000000000000000000000000000000000, 0.148874338981631210884826001129720,
    0.294392862701460198131126603103866, 0.433395394129247190799265943165784,
    0.562757134668604683339000099272694, 0.679409568299024406234327365114874,
    0.780817726586416897063717578345042, 0.865063366688984510732096688423493,
    0.930157491355708226001207180059508, 0.973906528517171720077964012084452,
    0.995657163025808080735527280689003];


pub const WGK1: Simd<f64,16> = Simd::from_array([0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
    0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
    0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
    0.123491976262065851077958109831074, 0.134709217311473325928054001771707,
    0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
    0.149445554002916905664936468389821, 0.147739104901338491374841515972068,
    0.142775938577060080797094273138717, 0.134709217311473325928054001771707,
    0.123491976262065851077958109831074, 0.109387158802297641899210590325805]);

pub const WGK2: Simd<f64,4> = Simd::from_array([0.093125454583697605535065465083366, 0.075039674810919952767043140916190,
    0.054755896574351996031381300244580, 0.032558162307964727478818972459390]);

pub const WGK3 : f64 = 0.011694638867371874278064396062192;


pub const WG1: Simd<f64,16> = Simd::from_array([0.0, 0.066671344308688137593568809893332, 0.0, 0.149451349150580593145776339657697,
    0.0, 0.219086362515982043995534934228163, 0.0, 0.269266719309996355091226921569469,
    0.0, 0.295524224714752870173892994651338, 0.0, 0.295524224714752870173892994651338,
    0.0, 0.269266719309996355091226921569469, 0.0, 0.219086362515982043995534934228163]);

pub const WG2: Simd<f64,4> = Simd::from_array(
    [0.0, 0.149451349150580593145776339657697, 0.0, 0.066671344308688137593568809893332]);



impl Qk for Qk21Simd2 {
    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, ) -> (f64, f64, f64, f64) {
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let (fv1,fv2,fv3) = fvec_simd(f, centr, hlgth);
        let resk = (fv1 * WGK1).reduce_sum() + ( fv2 * WGK2).reduce_sum() + fv3 * WGK3;
        let reskh = resk * 0.5;
        let (reskhs1,reskhs2) = (Simd::from_array([reskh;16]),Simd::from_array([reskh;4]));
        let mut resabs = (fv1.abs() * WGK1).reduce_sum() + ( fv2.abs() * WGK2).reduce_sum() + fv3.abs() * WGK3;
        let resg = (fv1 * WG1).reduce_sum() + ( fv2 * WG2).reduce_sum();
        let mut resasc = ( WGK1 * ( fv1 - reskhs1).abs()).reduce_sum() +
            ( WGK2 * ( fv2 - reskhs2).abs()).reduce_sum() + ( WGK3 * ( fv3 - reskh).abs());

        let result = resk * hlgth;
        resabs = resabs * dhlgth;
        resasc = resasc * dhlgth;
        let mut abserr = ((resk - resg) * hlgth).abs();
        if resasc != 0.0 && abserr != 0.0 {
            abserr = resasc * 1.0_f64.min((200.0 * abserr / resasc).powf(1.5));
        }
        if resabs > UFLOW / (50.0 * EPMACH) {
            abserr = abserr.max((EPMACH * 50.0) * resabs);
        }

        (result, abserr, resabs, resasc)
    }
}


pub fn fvec_simd(f : &dyn Fn(f64)->f64, centr : f64, hlgth : f64) -> (Simd<f64,16>,Simd<f64,4>,f64){
    (Simd::from_array([f(centr + hlgth * XGK[0]), f(centr + hlgth * XGK[1]), f(centr + hlgth * XGK[2]),
        f(centr + hlgth * XGK[3]), f(centr + hlgth * XGK[4]), f(centr + hlgth * XGK[5]),
        f(centr + hlgth * XGK[6]), f(centr + hlgth * XGK[7]), f(centr + hlgth * XGK[8]),
        f(centr + hlgth * XGK[9]), f(centr + hlgth * XGK[10]),f(centr + hlgth * XGK[11]),
        f(centr + hlgth * XGK[12]), f(centr + hlgth * XGK[13]), f(centr + hlgth * XGK[14]),
        f(centr + hlgth * XGK[15])]),
     Simd::from_array([f(centr + hlgth * XGK[16]), f(centr + hlgth * XGK[17]),
         f(centr + hlgth * XGK[18]), f(centr + hlgth * XGK[19])]),
     f(centr + hlgth * XGK[20]))
}





