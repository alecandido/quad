use crate::qk::*;

pub struct Qk15{}

///           f      : f64
///                     function
///
///           a      : f64
///                    lower limit of integration
///
///           b      : f64
///                    upper limit of integration
///
///         on return
///              result : f64
///                       approximation to the integral i result is computed by applying
///                       the 15-point kronrod rule (resk) obtained by optimal addition
///                       of abscissae to the7-point gauss rule(resg).
///
///              abserr : f64
///                       estimate of the modulus of the absolute error, which should not
///                       exceed abs(i-result)
///
///              resabs : f64
///                       approximation to the integral j
///
///              resasc : f64
///                       approximation to the integral of abs(f-i/(b-a)) over (a,b)
///
///           The abscissae and weights are given for the interval (-1,1).
///           Because of symmetry only the positive abscissae and their
///           corresponding weights are given.
///
///           xgk    : abscissae of the 15-point kronrod rule
///                    xgk(2), xgk(4), ...  abscissae of the 7-point
///                    gauss rule
///                    xgk(1), xgk(3), ...  abscissae which are optimally
///                    added to the 7-point gauss rule
///
///           wgk    : weights of the 15-point kronrod rule
///
///           wg     : weights of the 7-point gauss rule
///
///
///           Gauss quadrature weights and kronrod quadrature abscissae and weights
///           as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
///           bell labs, nov. 1981.
///
///
///

const XGK : [f64;8] = [0.991455371120812639206854697526329, 0.949107912342758524526189684047851,
                       0.864864423359769072789712788640926, 0.741531185599394439863864773280788,
                       0.586087235467691130294144838258730, 0.405845151377397166906606412076961,
                       0.207784955007898467600689403773245, 0.000000000000000000000000000000000];

const WGK : [f64;8] = [0.022935322010529224963732008058970, 0.063092092629978553290700663189204,
                       0.104790010322250183839876322541518, 0.140653259715525918745189590510238,
                       0.169004726639267902826583426598550, 0.190350578064785409913256402421014,
                       0.204432940075298892414161999234649, 0.209482141084727828012999174891714];

const WG : [f64;4] = [0.129484966168869693270611432679082, 0.279705391489276667901467771423780,
                      0.381830050505118944950369775488975, 0.417959183673469387755102040816327];

impl Qk for Qk15{
    fn integrate(&self,f : &dyn Fn(f64)->f64, a : f64, b : f64,) -> (f64, f64, f64, f64){
        let hlgth : f64 = 0.5*(b-a);
        let dhlgth : f64 = hlgth.abs();
        let centr : f64 = 0.5 * (b+a);

        let mut fv1 : Vec<f64> = vec![0.0;8];
        let mut fv2 : Vec<f64> = vec![0.0;8];

        //      compute the 15-point kronrod approximation to
        //      the integral, and estimate the absolute error.

        let fc : f64 = f(centr);
        let mut resg = fc * WG[3];
        let mut resk = fc * WGK[7];
        let mut resabs = resk.abs();

        for j in 1..4{
            let jtw = j * 2;
            let absc = hlgth * XGK[jtw-1];
            let fval1 : f64 = f(centr - absc);
            let fval2 : f64 = f(centr + absc);
            fv1[jtw-1] = fval1;
            fv2[jtw-1] = fval2;
            let fsum = fval1+fval2;
            resg += WG[j-1] * fsum;
            resk += WGK[jtw-1] * fsum;
            resabs += WGK[jtw-1] * (fval1.abs() + fval2.abs());
        }

        for j in 1..5{
            let jtwm1 = j*2-1;
            let absc = hlgth * XGK[jtwm1-1];
            let fval1 : f64 = f(centr - absc);
            let fval2 : f64 = f(centr + absc);
            fv1[jtwm1-1] = fval1;
            fv2[jtwm1-1] = fval2;
            let fsum = fval1 + fval2;
            resk += WGK[jtwm1-1] * fsum;
            resabs += WGK[jtwm1-1] * (fval1.abs() + fval2.abs());
        }
        let reskh = resk * 0.5;
        let mut resasc = WGK[7] * (fc - reskh).abs();

        for j in 1..8{
            resasc += WGK[j-1] * ((fv1[j-1] - reskh).abs() + (fv2[j-1] - reskh).abs());
        }

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

        (result,abserr,resabs,resasc)
    }
    }


