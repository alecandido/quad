use crate::qk::*;

pub struct Qk21{}

//           f      : f64
//                     function
//
//           a      : f64
//                    lower limit of integration
//
//           b      : f64
//                    upper limit of integration
//
//         on return
//              result : f64
//                       approximation to the integral i result is computed by applying
//                       the 21-point kronrod rule (resk) obtained by optimal addition
//                       of abscissae to the 10-point gauss rule(resg).
//
//              abserr : f64
//                       estimate of the modulus of the absolute error, which should not
//                       exceed abs(i-result)
//
//              resabs : f64
//                       approximation to the integral j
//
//              resasc : f64
//                       approximation to the integral of abs(f-i/(b-a)) over (a,b)
//
//           The abscissae and weights are given for the interval (-1,1).
//           Because of symmetry only the positive abscissae and their
//           corresponding weights are given.
//
//           xgk    : abscissae of the 21-point kronrod rule
//                    xgk(2), xgk(4), ...  abscissae of the 10-point
//                    gauss rule
//                    xgk(1), xgk(3), ...  abscissae which are optimally
//                    added to the 10-point gauss rule
//
//           wgk    : weights of the 21-point kronrod rule
//
//           wg     : weights of the 10-point gauss rule
//
//
//           Gauss quadrature weights and kronrod quadrature abscissae and weights
//           as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
//           bell labs, nov. 1981.
//
//
//

const XGK : [f64;11] = [0.995657163025808080735527280689003, 0.973906528517171720077964012084452,
                        0.930157491355708226001207180059508, 0.865063366688984510732096688423493,
                        0.780817726586416897063717578345042, 0.679409568299024406234327365114874,
                        0.562757134668604683339000099272694, 0.433395394129247190799265943165784,
                        0.294392862701460198131126603103866, 0.148874338981631210884826001129720,
                        0.000000000000000000000000000000000];

const WGK : [f64;11] = [0.011694638867371874278064396062192, 0.032558162307964727478818972459390,
                       0.054755896574351996031381300244580, 0.075039674810919952767043140916190,
                       0.093125454583697605535065465083366, 0.109387158802297641899210590325805,
                       0.123491976262065851077958109831074, 0.134709217311473325928054001771707,
                       0.142775938577060080797094273138717, 0.147739104901338491374841515972068,
                       0.149445554002916905664936468389821];


const WG : [f64;5] = [0.066671344308688137593568809893332, 0.149451349150580593145776339657697,
                      0.219086362515982043995534934228163, 0.269266719309996355091226921569469,
                      0.295524224714752870173892994651338];


impl Qk for Qk21 {
    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, ) -> (f64, f64, f64, f64) {
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let mut fv1: Vec<f64> = vec![0.0; 10];
        let mut fv2: Vec<f64> = vec![0.0; 10];

        //compute the 21-point kronrod approximation to
        //the integral, and estimate the absolute error.

        let mut resg = 0.0;
        let fc : f64 = f(centr);
        let mut resk = WGK[10] * fc;
        let mut resabs = resk.abs();

        for j in 1..6 {
            let jtw = 2 * j;
            let absc = hlgth * XGK[jtw - 1];
            let fval1 : f64 = f(centr - absc);
            let fval2 : f64 = f(centr + absc);
            fv1[jtw - 1] = fval1;
            fv2[jtw - 1] = fval2;
            let fsum = fval1 + fval2;
            resg += WG[j - 1] * fsum;
            resk += WGK[jtw - 1] * fsum;
            resabs += WGK[jtw - 1] * (fval1.abs() + fval2.abs());
        }

        for j in 1..6 {
            let jtwm1 = 2 * j - 1;
            let absc = hlgth * XGK[jtwm1 - 1];
            let fval1 : f64 = f(centr - absc);
            let fval2 : f64 = f(centr + absc);
            fv1[jtwm1 - 1] = fval1;
            fv2[jtwm1 - 1] = fval2;
            let fsum = fval1 + fval2;
            resk += WGK[jtwm1 - 1] * fsum;
            resabs += WGK[jtwm1 - 1] * (fval1.abs() + fval2.abs());
        }

        let reskh = resk * 0.5;
        let mut resasc = WGK[10] * (fc - reskh).abs();

        for j in 1..11 {
            resasc += WGK[j - 1] * ((fv1[j - 1] - reskh).abs() + (fv2[j - 1] - reskh).abs());
        }

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
