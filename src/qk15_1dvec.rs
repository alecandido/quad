extern crate nalgebra as na;

use std::time::Instant;
use na::{SVector, vector};
use crate::qk::*;


pub struct Qk151DVec {}
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
///                     the 15-point kronrod rule (resk) obtained by optimal addition
///                     of abscissae to the7-point gauss rule(resg).
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
///         xgk     :   abscissae of the 15-point kronrod rule
///                     xgk(2), xgk(4), ...  abscissae of the 7-point
///                     gauss rule
///                     xgk(1), xgk(3), ...  abscissae which are optimally
///                     added to the 7-point gauss rule
///
///         wgk     :   weights of the 15-point kronrod rule
///
///         wg      :   weights of the 7-point gauss rule
///
///
///     Gauss quadrature weights and kronrod quadrature abscissae and weights
///     as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
///     bell labs, nov. 1981.
///
///
///


const XGK: SVector<f64,15> = vector![- 0.991455371120812639206854697526329,
            -0.949107912342758524526189684047851, -0.864864423359769072789712788640926,
            -0.741531185599394439863864773280788, -0.586087235467691130294144838258730,
            -0.405845151377397166906606412076961, -0.207784955007898467600689403773245,
            0.000000000000000000000000000000000, 0.207784955007898467600689403773245,
            0.405845151377397166906606412076961, 0.586087235467691130294144838258730,
            0.741531185599394439863864773280788, 0.864864423359769072789712788640926,
            0.949107912342758524526189684047851, 0.991455371120812639206854697526329];

const WGK: SVector<f64,15> = vector![0.022935322010529224963732008058970,
            0.063092092629978553290700663189204, 0.104790010322250183839876322541518,
            0.140653259715525918745189590510238, 0.169004726639267902826583426598550,
            0.190350578064785409913256402421014, 0.204432940075298892414161999234649,
            0.209482141084727828012999174891714, 0.204432940075298892414161999234649,
            0.190350578064785409913256402421014, 0.169004726639267902826583426598550,
            0.140653259715525918745189590510238, 0.104790010322250183839876322541518,
            0.063092092629978553290700663189204, 0.022935322010529224963732008058970];


const WG: SVector<f64,15> = vector![0.0, 0.129484966168869693270611432679082, 0.0,
            0.279705391489276667901467771423780, 0.0, 0.381830050505118944950369775488975, 0.0,
            0.417959183673469387755102040816327, 0.0, 0.381830050505118944950369775488975, 0.0,
            0.279705391489276667901467771423780, 0.0, 0.129484966168869693270611432679082, 0.0];


impl Qk for Qk151DVec {
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64) -> (f64, f64, f64, f64){
        //let start = Instant::now();
        let hlgth : f64 = 0.5*(b-a);
        let dhlgth : f64 = hlgth.abs();
        let centr : f64 = 0.5 * (b+a);

        //println!("first initialization :{:?}", start.elapsed());

        let fv = fvec(f,centr,hlgth);
        let resk = WGK.dot(&fv);
        let mut resabs = WGK.dot(&fv.abs());
        let resg = WG.dot(&fv);



        //println!("cdot :{:?}", start.elapsed());

        let reskh = resk * 0.5;
        let reskh_vec = vector![reskh,reskh,reskh,reskh,reskh,reskh,reskh,reskh,
        reskh,reskh,reskh,reskh,reskh,reskh,reskh];
        let lv = (fv - reskh_vec).abs();

        //println!("other vector initialization :{:?}", start.elapsed());

        let mut resasc = WGK.dot(&lv);


        //println!("second cdot :{:?}", start.elapsed());

        let result = resk * hlgth;
        resabs = resabs * dhlgth;
        resasc = resasc * dhlgth;
        let mut abserr = ((resk-resg) * hlgth).abs();
        if resasc != 0.0 && abserr != 0.0 {
            abserr = resasc * 1.0_f64.min( (200.0 * abserr/resasc).powf(1.5));
        }
        if resabs > UFLOW /(50.0 * EPMACH) {
            abserr = abserr.max((EPMACH * 50.0) * resabs);
        }

        //println!("last computation :{:?}", start.elapsed());

        (result,abserr,resabs,resasc)
    }
}



pub fn fvec(f : &dyn Fn(f64)->f64, centr : f64, hlgth : f64) -> SVector<f64,15>{
    vector![f(centr + hlgth * XGK[0]), f(centr + hlgth * XGK[1]), f(centr + hlgth * XGK[2]),
        f(centr + hlgth * XGK[3]), f(centr + hlgth * XGK[4]), f(centr + hlgth * XGK[5]),
        f(centr + hlgth * XGK[6]), f(centr + hlgth * XGK[7]), f(centr + hlgth * XGK[8]),
        f(centr + hlgth * XGK[9]), f(centr + hlgth * XGK[10]),f(centr + hlgth * XGK[11]),
        f(centr + hlgth * XGK[12]), f(centr + hlgth * XGK[13]), f(centr + hlgth * XGK[14])]
}








