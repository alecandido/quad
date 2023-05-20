use crate::qk::*;
use na::{DVector, dvector, SVector, vector};


pub struct Qk61VecNa {}
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

const XGK : [f64;31] = [0.999484410050490637571325895705811, 0.996893484074649540271630050918695,
    0.991630996870404594858628366109486, 0.983668123279747209970032581605663,
    0.973116322501126268374693868423707, 0.960021864968307512216871025581798,
    0.944374444748559979415831324037439, 0.926200047429274325879324277080474,
    0.905573307699907798546522558925958, 0.882560535792052681543116462530226,
    0.857205233546061098958658510658944, 0.829565762382768397442898119732502,
    0.799727835821839083013668942322683, 0.767777432104826194917977340974503,
    0.733790062453226804726171131369528, 0.697850494793315796932292388026640,
    0.660061064126626961370053668149271, 0.620526182989242861140477556431189,
    0.579345235826361691756024932172540, 0.536624148142019899264169793311073,
    0.492480467861778574993693061207709, 0.447033769538089176780609900322854,
    0.400401254830394392535476211542661, 0.352704725530878113471037207089374,
    0.304073202273625077372677107199257, 0.254636926167889846439805129817805,
    0.204525116682309891438957671002025, 0.153869913608583546963794672743256,
    0.102806937966737030147096751318001, 0.051471842555317695833025213166723,
    0.000000000000000000000000000000000];

const WGK : [f64;31] =[0.001389013698677007624551591226760, 0.003890461127099884051267201844516,
    0.006630703915931292173319826369750, 0.009273279659517763428441146892024,
    0.011823015253496341742232898853251, 0.014369729507045804812451432443580,
    0.016920889189053272627572289420322, 0.019414141193942381173408951050128,
    0.021828035821609192297167485738339, 0.024191162078080601365686370725232,
    0.026509954882333101610601709335075, 0.028754048765041292843978785354334,
    0.030907257562387762472884252943092, 0.032981447057483726031814191016854,
    0.034979338028060024137499670731468, 0.036882364651821229223911065617136,
    0.038678945624727592950348651532281, 0.040374538951535959111995279752468,
    0.041969810215164246147147541285970, 0.043452539701356069316831728117073,
    0.044814800133162663192355551616723, 0.046059238271006988116271735559374,
    0.047185546569299153945261478181099, 0.048185861757087129140779492298305,
    0.049055434555029778887528165367238, 0.049795683427074206357811569379942,
    0.050405921402782346840893085653585, 0.050881795898749606492297473049805,
    0.051221547849258772170656282604944, 0.051426128537459025933862879215781,
    0.051494729429451567558340433647099,];

const WG : [f64;15] = [0.007968192496166605615465883474674, 0.018466468311090959142302131912047,
    0.028784707883323369349719179611292, 0.038799192569627049596801936446348,
    0.048402672830594052902938140422808, 0.057493156217619066481721689402056,
    0.065974229882180495128128515115962, 0.073755974737705206268243850022191,
    0.080755895229420215354694938460530, 0.086899787201082979802387530715126,
    0.092122522237786128717632707087619, 0.096368737174644259639468626351810,
    0.099593420586795267062780282103569, 0.101762389748405504596428952168554,
    0.102852652893558840341285636705415];


impl Qk61VecNa {
    fn integrate<F>(&self, f: F, a: f64, b: f64) -> () //(DVector<f64>, DVector<f64>, DVector<f64>, DVector<f64>)
    where F : Fn(f64) -> DVector<f64>{
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let a =  dvector![f(XGK[0]),f(XGK[1])];
        let b = dvector![2.0,1.0];
        let c = a*b;
        println!("mmm {:?}",c);

        /*

        let mut fv1 = vec![];
        let mut fv2 = vec![];

        //compute the 61-point kronrod approximation to
        //the integral, and estimate the absolute error.

        let mut resg = dvector![];
        let fc   = f(centr);
        let mut resk = WGK[30] * fc;
        let mut resabs = resk.abs();

        for j in 1..16 {
            let jtw = 2 * j;
            let absc = hlgth * XGK[jtw - 1];
            let fval1 : DVector<f64> = f(centr - absc);
            let fval2 : DVector<f64> = f(centr + absc);
            fv1.push(fval1.clone());
            fv2.push(fval2.clone());
            let fsum  : DVector<f64> = fval1 + fval2;
            resg += WG[j - 1] * fsum;
            resk += WGK[jtw - 1] * fsum;
            resabs += WGK[jtw - 1] * (fval1.abs() + fval2.abs());


        }

        for j in 1..16 {
            let jtwm1 = 2 * j - 1;
            let absc = hlgth * XGK[jtwm1 - 1];
            let fval1 = f(centr - absc);
            let fval2 = f(centr + absc);
            fv1.push(fval1);
            fv2.push(fval2);
            let fsum = fval1 + fval2;
            resk += WGK[jtwm1 - 1] * fsum;
            resabs += WGK[jtwm1 - 1] * (fval1.abs() + fval2.abs());
        }

        let reskh = 0.5 * resk;
        let mut resasc = WGK[30] * (fc - reskh).abs();

        for j in 1..16 {

            resasc += WGK[2*j - 1] * ((fv1[j - 1] - reskh).abs() + (fv2[j - 1] -
                reskh).abs()) + WGK[2*(j - 1)] * ((fv1[j + 14] - reskh).abs() + (fv2[j + 14] -
                reskh).abs())  ;


        }

        let result = resk * hlgth;
        resabs = resabs * dhlgth;
        resasc = resasc * dhlgth;

        let mut abserr = ((resk - resg) * hlgth).abs();
        for k in 0..4{
            if resasc[k] != 0.0 && abserr[k] != 0.0 {
                abserr[k] = resasc[k] * 1.0_f64.min((200.0 * abserr[k] / resasc[k]).powf(1.5));
            }
            if resabs[k] > UFLOW / (50.0 * EPMACH) {
                abserr[k] = abserr[k].max((EPMACH * 50.0) * resabs[k]);
            }
        }


        (result, abserr, resabs, resasc)

         */
    }
}


#[cfg(test)]
mod tests {
    use std::time::Instant;
    use na::{dvector, vector};
    use crate::qk61::Qk61;
    use crate::qk61_4vec::Qk614Vec;
    use crate::qk61_4vec_na::Qk614VecNa;
    use crate::qk61_vec_scipy::Qk61VecNa;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f1 = |x:f64| x.cos();
        let f2 = |x:f64| x.sin();
        let f3 = |x:f64| f1(x) + f2(x);
        let f4 = |x:f64| - 2.0 * f3(x);
        let f = |x:f64| [f1(x),f2(x),f3(x),f4(x)];
        let ff = |x:f64| vector![f1(x),f2(x),f3(x),f4(x)];
        let ff2 = |x:f64| dvector![f1(x),f2(x),f3(x),f4(x)];

        let a = 0.0;
        let b = 1.0;
        let qk = Qk61 {};
        let qk_vec = Qk614Vec {};
        let qk_vec_na = Qk614VecNa {};
        let qk_vec_scipy = Qk61VecNa{};



        let mut res1 = (0.0,0.0,0.0,0.0);
        let mut res2 = res1.clone();
        let mut res3 = res1.clone();
        let mut res4 = res1.clone();
        let mut res_vec = ([0.0;4],[0.0;4],[0.0;4],[0.0;4]);
        let mut res_vec_na = (vector![0.0,0.0,0.0,0.0],
                              vector![0.0,0.0,0.0,0.0],vector![0.0,0.0,0.0,0.0],
                              vector![0.0,0.0,0.0,0.0]);
        let mut res_vec_na2 ;

        for k in 0..100 {
            let start = Instant::now();
            res1 = qk.integrate(&f1, a, b);
            res2 = qk.integrate(&f2,a,b);
            res3 = qk.integrate(&f3,a,b);
            res4 = qk.integrate(&f4,a,b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            res_vec = qk_vec.integrate(&f, a, b);
            println!("vec {:?}", start.elapsed());
            let start = Instant::now();
            //res_vec_na = qk_vec_na.integrate(&ff, a, b);
            println!("vec na {:?}", start.elapsed());
            let start = Instant::now();
            res_vec_na2 = qk_vec_scipy.integrate(&ff2, a, b);
            println!("vec na2 {:?}", start.elapsed());

        }
        println!("normal {:?},{:?},{:?},{:?}",res1,res2,res3,res4);
        println!("vec {:?}",res_vec);
        println!("vec na {:?}",res_vec_na);
    }
}

