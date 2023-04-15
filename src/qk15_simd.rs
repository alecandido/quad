use std::simd::{Simd, SimdFloat};
use crate::qk::*;


pub struct Qk15Simd {}
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


const XGK: Simd<f64,16> = Simd::from_array([- 0.991455371120812639206854697526329,
            -0.949107912342758524526189684047851, -0.864864423359769072789712788640926,
            -0.741531185599394439863864773280788, -0.586087235467691130294144838258730,
            -0.405845151377397166906606412076961, -0.207784955007898467600689403773245,
            0.000000000000000000000000000000000, 0.207784955007898467600689403773245,
            0.405845151377397166906606412076961, 0.586087235467691130294144838258730,
            0.741531185599394439863864773280788, 0.864864423359769072789712788640926,
            0.949107912342758524526189684047851, 0.991455371120812639206854697526329,
            0.000000000000000000000000000000000]);

const WGK: Simd<f64,16> = Simd::from_array([0.022935322010529224963732008058970,
            0.063092092629978553290700663189204, 0.104790010322250183839876322541518,
            0.140653259715525918745189590510238, 0.169004726639267902826583426598550,
            0.190350578064785409913256402421014, 0.204432940075298892414161999234649,
            0.209482141084727828012999174891714, 0.204432940075298892414161999234649,
            0.190350578064785409913256402421014, 0.169004726639267902826583426598550,
            0.140653259715525918745189590510238, 0.104790010322250183839876322541518,
            0.063092092629978553290700663189204, 0.022935322010529224963732008058970,
            0.000000000000000000000000000000000]);


const WG: Simd<f64,16> = Simd::from_array([0.0, 0.129484966168869693270611432679082, 0.0,
            0.279705391489276667901467771423780, 0.0, 0.381830050505118944950369775488975, 0.0,
            0.417959183673469387755102040816327, 0.0, 0.381830050505118944950369775488975, 0.0,
            0.279705391489276667901467771423780, 0.0, 0.129484966168869693270611432679082, 0.0,
            0.000000000000000000000000000000000]);


impl Qk for Qk15Simd {
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64) -> (f64, f64, f64, f64){
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let fv = fvec_simd(f, centr, hlgth);
        let resk = (fv * WGK).reduce_sum();
        let reskh = resk * 0.5;
        let reskhs = Simd::from_array([reskh;16]);
        let mut resabs = (fv.abs() * WGK).reduce_sum();
        let resg = ( fv * WG).reduce_sum();
        let mut resasc = ( WGK * ( fv - reskhs).abs()).reduce_sum();

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



pub fn fvec_simd(f : &dyn Fn(f64)->f64, centr : f64, hlgth : f64) -> Simd<f64,16>{
    Simd::from_array([f(centr + hlgth * XGK[0]), f(centr + hlgth * XGK[1]), f(centr + hlgth * XGK[2]),
        f(centr + hlgth * XGK[3]), f(centr + hlgth * XGK[4]), f(centr + hlgth * XGK[5]),
        f(centr + hlgth * XGK[6]), f(centr + hlgth * XGK[7]), f(centr + hlgth * XGK[8]),
        f(centr + hlgth * XGK[9]), f(centr + hlgth * XGK[10]),f(centr + hlgth * XGK[11]),
        f(centr + hlgth * XGK[12]), f(centr + hlgth * XGK[13]), f(centr + hlgth * XGK[14]),
        0.000000000000000000000000000000000])
}

#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qk15::Qk15;
    use crate::qk15_simd::Qk15Simd;
    use crate::qk21::Qk21;
    use crate::qk21_simd::Qk21Simd;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.sin();
        let a = 0.0;
        let b = 1.0;
        let qks = Qk15Simd {};
        let qk = Qk15{};

        let mut res2 = (0.0,0.0,0.0,0.0);
        let mut res1 = res2.clone();

        for k in 0..100 {
            let start = Instant::now();
            res2 = qk.integrate(&f, a, b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            res1 = qks.integrate(&f, a, b);
            println!("simd {:?}", start.elapsed());
        }
        println!("{:?}",res1);
        println!("{:?}",res2);
    }
}






