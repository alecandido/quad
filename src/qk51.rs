use crate::qk::*;

pub struct Qk51{}

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
//                       the 51-point kronrod rule (resk) obtained by optimal addition
//                       of abscissae to the 25-point gauss rule(resg).
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
//           xgk    : abscissae of the 51-point kronrod rule
//                    xgk(2), xgk(4), ...  abscissae of the 25-point
//                    gauss rule
//                    xgk(1), xgk(3), ...  abscissae which are optimally
//                    added to the 25-point gauss rule
//
//           wgk    : weights of the 51-point kronrod rule
//
//           wg     : weights of the 25-point gauss rule
//
//
//           Gauss quadrature weights and kronrod quadrature abscissae and weights
//           as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
//           bell labs, nov. 1981.
//
//
//

const XGK : [f64;26] = [0.999262104992609834193457486540341, 0.995556969790498097908784946893902,
                        0.988035794534077247637331014577406, 0.976663921459517511498315386479594,
                        0.961614986425842512418130033660167, 0.942974571228974339414011169658471,
                        0.920747115281701561746346084546331, 0.894991997878275368851042006782805,
                        0.865847065293275595448996969588340, 0.833442628760834001421021108693570,
                        0.797873797998500059410410904994307, 0.759259263037357630577282865204361,
                        0.717766406813084388186654079773298, 0.673566368473468364485120633247622,
                        0.626810099010317412788122681624518, 0.577662930241222967723689841612654,
                        0.526325284334719182599623778158010, 0.473002731445714960522182115009192,
                        0.417885382193037748851814394594572, 0.361172305809387837735821730127641,
                        0.303089538931107830167478909980339, 0.243866883720988432045190362797452,
                        0.183718939421048892015969888759528, 0.122864692610710396387359818808037,
                        0.061544483005685078886546392366797, 0.000000000000000000000000000000000];

const WGK : [f64;26] = [0.001987383892330315926507851882843, 0.005561932135356713758040236901066,
                        0.009473973386174151607207710523655, 0.013236229195571674813656405846976,
                        0.016847817709128298231516667536336, 0.020435371145882835456568292235939,
                        0.024009945606953216220092489164881, 0.027475317587851737802948455517811,
                        0.030792300167387488891109020215229, 0.034002130274329337836748795229551,
                        0.037116271483415543560330625367620, 0.040083825504032382074839284467076,
                        0.042872845020170049476895792439495, 0.045502913049921788909870584752660,
                        0.047982537138836713906392255756915, 0.050277679080715671963325259433440,
                        0.052362885806407475864366712137873, 0.054251129888545490144543370459876,
                        0.055950811220412317308240686382747, 0.057437116361567832853582693939506,
                        0.058689680022394207961974175856788, 0.059720340324174059979099291932562,
                        0.060539455376045862945360267517565, 0.061128509717053048305859030416293,
                        0.061471189871425316661544131965264, 0.061580818067832935078759824240066];

const WG : [f64;13] = [0.011393798501026287947902964113235, 0.026354986615032137261901815295299,
                       0.040939156701306312655623487711646, 0.054904695975835191925936891540473,
                       0.068038333812356917207187185656708, 0.080140700335001018013234959669111,
                       0.091028261982963649811497220702892, 0.100535949067050644202206890392686,
                       0.108519624474263653116093957050117, 0.114858259145711648339325545869556,
                       0.119455763535784772228178126512901, 0.122242442990310041688959518945852,
                       0.123176053726715451203902873079050];


impl Qk for Qk51 {
    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, ) -> (f64, f64, f64, f64) {
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let mut fv1: Vec<f64> = vec![0.0; 25];
        let mut fv2: Vec<f64> = vec![0.0; 25];

        //compute the 51-point kronrod approximation to
        //the integral, and estimate the absolute error.

        let fc : f64 = f(centr);
        let mut resg = WG[12] * fc;
        let mut resk = WGK[25] * fc;
        let mut resabs = resk.abs();


        for j in 1..13 {
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


        for j in 1..14 {
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
        let mut resasc = WGK[25] * (fc - reskh).abs();

        for j in 1..26 {
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
