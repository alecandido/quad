use std::iter::zip;
use std::simd::{Simd, SimdFloat};
use std::time::Instant;
use crate::qk::*;

pub struct Qk611DVec_Simd {}
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


const XGK: Simd<f64,64> = Simd::from_array([-0.999484410050490637571325895705811, -0.996893484074649540271630050918695,
    -0.991630996870404594858628366109486, -0.983668123279747209970032581605663,
    -0.973116322501126268374693868423707, -0.960021864968307512216871025581798,
    -0.944374444748559979415831324037439, -0.926200047429274325879324277080474,
    -0.905573307699907798546522558925958, -0.882560535792052681543116462530226,
    -0.857205233546061098958658510658944, -0.829565762382768397442898119732502,
    -0.799727835821839083013668942322683, -0.767777432104826194917977340974503,
    -0.733790062453226804726171131369528, -0.697850494793315796932292388026640,
    -0.660061064126626961370053668149271, -0.620526182989242861140477556431189,
    -0.579345235826361691756024932172540, -0.536624148142019899264169793311073,
    -0.492480467861778574993693061207709, -0.447033769538089176780609900322854,
    -0.400401254830394392535476211542661, -0.352704725530878113471037207089374,
    -0.304073202273625077372677107199257, -0.254636926167889846439805129817805,
    -0.204525116682309891438957671002025, -0.153869913608583546963794672743256,
    -0.102806937966737030147096751318001, -0.051471842555317695833025213166723,
    0.000000000000000000000000000000000, 0.051471842555317695833025213166723,
    0.102806937966737030147096751318001, 0.153869913608583546963794672743256,
    0.204525116682309891438957671002025, 0.254636926167889846439805129817805,
    0.304073202273625077372677107199257, 0.352704725530878113471037207089374,
    0.400401254830394392535476211542661, 0.447033769538089176780609900322854,
    0.492480467861778574993693061207709, 0.536624148142019899264169793311073,
    0.579345235826361691756024932172540, 0.620526182989242861140477556431189,
    0.660061064126626961370053668149271, 0.697850494793315796932292388026640,
    0.733790062453226804726171131369528, 0.767777432104826194917977340974503,
    0.799727835821839083013668942322683, 0.829565762382768397442898119732502,
    0.857205233546061098958658510658944, 0.882560535792052681543116462530226,
    0.905573307699907798546522558925958, 0.926200047429274325879324277080474,
    0.944374444748559979415831324037439, 0.960021864968307512216871025581798,
    0.973116322501126268374693868423707, 0.983668123279747209970032581605663,
    0.991630996870404594858628366109486, 0.996893484074649540271630050918695,
    0.999484410050490637571325895705811, 0.0, 0.0, 0.0]);


const WGK: Simd<f64,64> = Simd::from_array([0.001389013698677007624551591226760, 0.003890461127099884051267201844516,
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
    0.051494729429451567558340433647099, 0.051426128537459025933862879215781,
    0.051221547849258772170656282604944, 0.050881795898749606492297473049805,
    0.050405921402782346840893085653585, 0.049795683427074206357811569379942,
    0.049055434555029778887528165367238, 0.048185861757087129140779492298305,
    0.047185546569299153945261478181099, 0.046059238271006988116271735559374,
    0.044814800133162663192355551616723, 0.043452539701356069316831728117073,
    0.041969810215164246147147541285970, 0.040374538951535959111995279752468,
    0.038678945624727592950348651532281, 0.036882364651821229223911065617136,
    0.034979338028060024137499670731468, 0.032981447057483726031814191016854,
    0.030907257562387762472884252943092, 0.028754048765041292843978785354334,
    0.026509954882333101610601709335075, 0.024191162078080601365686370725232,
    0.021828035821609192297167485738339, 0.019414141193942381173408951050128,
    0.016920889189053272627572289420322, 0.014369729507045804812451432443580,
    0.011823015253496341742232898853251, 0.009273279659517763428441146892024,
    0.006630703915931292173319826369750, 0.003890461127099884051267201844516,
    0.001389013698677007624551591226760, 0.0, 0.0, 0.0]);

const WG: Simd<f64,64> = Simd::from_array([0.0, 0.007968192496166605615465883474674, 0.0, 0.018466468311090959142302131912047,
    0.0, 0.028784707883323369349719179611292, 0.0, 0.038799192569627049596801936446348,
    0.0, 0.048402672830594052902938140422808, 0.0, 0.057493156217619066481721689402056,
    0.0, 0.065974229882180495128128515115962, 0.0, 0.073755974737705206268243850022191,
    0.0, 0.080755895229420215354694938460530, 0.0, 0.086899787201082979802387530715126,
    0.0, 0.092122522237786128717632707087619, 0.0, 0.096368737174644259639468626351810,
    0.0, 0.099593420586795267062780282103569, 0.0, 0.101762389748405504596428952168554,
    0.0, 0.102852652893558840341285636705415, 0.0,
    0.102852652893558840341285636705415, 0.0, 0.101762389748405504596428952168554,
    0.0, 0.099593420586795267062780282103569, 0.0, 0.096368737174644259639468626351810,
    0.0, 0.092122522237786128717632707087619, 0.0, 0.086899787201082979802387530715126,
    0.0, 0.080755895229420215354694938460530, 0.0, 0.073755974737705206268243850022191,
    0.0, 0.065974229882180495128128515115962, 0.0, 0.057493156217619066481721689402056,
    0.0, 0.048402672830594052902938140422808, 0.0, 0.038799192569627049596801936446348,
    0.0, 0.028784707883323369349719179611292, 0.0, 0.018466468311090959142302131912047,
    0.0, 0.007968192496166605615465883474674, 0.0, 0.0, 0.0, 0.0]);






impl Qk611DVec_Simd {
    pub fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, ) -> (f64, f64, f64, f64) {
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let fv = fvec_simd(f, centr, hlgth);
        let resk = (fv * WGK).reduce_sum();
        let reskh = resk * 0.5;
        let reskhs = Simd::from_array([reskh;64]);
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


pub fn fvec_simd(f : &dyn Fn(f64)->f64, centr : f64, hlgth : f64) -> Simd<f64,64>{
    Simd::from_array([f(centr + hlgth * XGK[0]), f(centr + hlgth * XGK[1]), f(centr + hlgth * XGK[2]),
        f(centr + hlgth * XGK[3]), f(centr + hlgth * XGK[4]), f(centr + hlgth * XGK[5]),
        f(centr + hlgth * XGK[6]), f(centr + hlgth * XGK[7]), f(centr + hlgth * XGK[8]),
        f(centr + hlgth * XGK[9]), f(centr + hlgth * XGK[10]),f(centr + hlgth * XGK[11]),
        f(centr + hlgth * XGK[12]), f(centr + hlgth * XGK[13]), f(centr + hlgth * XGK[14]),
        f(centr + hlgth * XGK[15]), f(centr + hlgth * XGK[16]), f(centr + hlgth * XGK[17]),
        f(centr + hlgth * XGK[18]), f(centr + hlgth * XGK[19]), f(centr + hlgth * XGK[20]),
        f(centr + hlgth * XGK[21]), f(centr + hlgth * XGK[22]), f(centr + hlgth * XGK[23]),
        f(centr + hlgth * XGK[24]), f(centr + hlgth * XGK[25]),f(centr + hlgth * XGK[26]),
        f(centr + hlgth * XGK[27]), f(centr + hlgth * XGK[28]), f(centr + hlgth * XGK[29]),
        f(centr + hlgth * XGK[30]), f(centr + hlgth * XGK[31]), f(centr + hlgth * XGK[32]),
        f(centr + hlgth * XGK[33]), f(centr + hlgth * XGK[34]), f(centr + hlgth * XGK[35]),
        f(centr + hlgth * XGK[36]), f(centr + hlgth * XGK[37]), f(centr + hlgth * XGK[38]),
        f(centr + hlgth * XGK[39]), f(centr + hlgth * XGK[40]),f(centr + hlgth * XGK[41]),
        f(centr + hlgth * XGK[42]), f(centr + hlgth * XGK[43]), f(centr + hlgth * XGK[44]),
        f(centr + hlgth * XGK[45]), f(centr + hlgth * XGK[46]), f(centr + hlgth * XGK[47]),
        f(centr + hlgth * XGK[48]), f(centr + hlgth * XGK[49]), f(centr + hlgth * XGK[50]),
        f(centr + hlgth * XGK[51]), f(centr + hlgth * XGK[52]), f(centr + hlgth * XGK[53]),
        f(centr + hlgth * XGK[54]), f(centr + hlgth * XGK[55]),f(centr + hlgth * XGK[56]),
        f(centr + hlgth * XGK[57]), f(centr + hlgth * XGK[58]), f(centr + hlgth * XGK[59]),
        f(centr + hlgth * XGK[60]), 0.0, 0.0, 0.0])
}


#[cfg(test)]
mod tests {
    use std::simd::Simd;
    use std::time::Instant;
    use crate::qk61_simd::Qk611DVec_Simd;
    use crate::qk::Qk;
    use rgsl::integration::qk61;
    use crate::qk61::Qk61;
    use crate::qk61_blas::Qk61Blas;
    use cblas::ddot;
    use crate::qk61_simd2::Qk61Simd2;
    use crate::qk61_simd3::Qk61Simd3;

    #[test]
    fn test(){
        unsafe{
            let f = |x:f64| x.cos();
            let a = 0.0;
            let b = 1000.0;
            let qks = Qk611DVec_Simd{};
            let qk = Qk61{};
            let qk_blas = Qk61Blas{};
            let qks2 = Qk61Simd2{};
            let qks3 = Qk61Simd3{};

            let mut rgsl_res = (0.0,0.0,0.0,0.0);
            let mut my_res = rgsl_res.clone();
            let mut res = rgsl_res.clone();
            let mut res_blas = rgsl_res.clone();
            let mut my_res2 = rgsl_res.clone();
            let mut my_res3 = rgsl_res.clone();
            let mut res_new = rgsl_res.clone();

            for k in 0..100 {
                let start = Instant::now();
                rgsl_res = qk61(f,a,b);
                println!("rgsl {:?}", start.elapsed());
                let start = Instant::now();
                my_res = qks.integrate(&f, a, b);
                println!("simd {:?}", start.elapsed());

                let start = Instant::now();
                res_new = qk.integrate(&f,a,b);
                println!("normal {:?}", start.elapsed());

                let start = Instant::now();
                let (mut result,mut abserr,mut resabs, mut resasc) = (0.0,0.0,0.0,0.0);
                qk.integrate2(f,a,b,&mut result, &mut abserr, &mut resabs, &mut resasc);
                println!("normal new {:?}", start.elapsed());

                let start = Instant::now();
                res_blas = qk_blas.integrate(&f,a,b);
                println!("blas {:?}", start.elapsed());
//
                let start = Instant::now();
                my_res2 = qks2.integrate(&f,a,b);
                println!("simd2 {:?}", start.elapsed());
//
                let start = Instant::now();
                my_res3 = qks3.integrate(&f,a,b);
                println!("simd3 {:?}", start.elapsed());

            }
            println!("simd : {:?}", my_res);
            println!("rgsl : {:?}",rgsl_res);
            println!("normal new : {:?}", res_new);
            println!("blas : {:?}", res_blas);
            println!("simd2 : {:?}", my_res2);
            println!("simd3 : {:?}", my_res3);


            /*
            let mut x = [1.0,2.0,3.0];
            let mut y = [2.0,3.0,4.0];
            cblas::daxpy(3,-1.0, & x, 1, &mut y, 1);
            println!("{:?}", y);

             */
        }

    }


}