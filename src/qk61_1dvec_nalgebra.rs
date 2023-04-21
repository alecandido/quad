extern crate nalgebra as na;
/*
use std::iter::zip;
use std::simd::{Simd, SimdFloat};
use na::SVector;
use crate::qk::*;
use simba::simd::f64x4;
use crate::na::SimdComplexField;

pub struct Qk61Nalgebra {}
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


pub const XGK: [f64;61] = [-0.999484410050490637571325895705811, -0.996893484074649540271630050918695,
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
    0.999484410050490637571325895705811];



pub const WGK1: [[f64;4];8] = [[0.001389013698677007624551591226760,
    0.006630703915931292173319826369750, 0.011823015253496341742232898853251,
    0.016920889189053272627572289420322],
    [0.021828035821609192297167485738339,0.026509954882333101610601709335075,
        0.030907257562387762472884252943092,0.034979338028060024137499670731468],
    [0.038678945624727592950348651532281,
    0.041969810215164246147147541285970, 0.044814800133162663192355551616723,
    0.047185546569299153945261478181099],
    [0.049055434555029778887528165367238,
    0.050405921402782346840893085653585, 0.051221547849258772170656282604944,
    0.051494729429451567558340433647099],
    [0.051221547849258772170656282604944,
    0.050405921402782346840893085653585, 0.049055434555029778887528165367238,
    0.047185546569299153945261478181099],
    [0.044814800133162663192355551616723,
    0.041969810215164246147147541285970, 0.038678945624727592950348651532281,
    0.034979338028060024137499670731468],
    [0.030907257562387762472884252943092,
    0.026509954882333101610601709335075, 0.021828035821609192297167485738339,
    0.016920889189053272627572289420322],
    [0.011823015253496341742232898853251,
    0.006630703915931292173319826369750, 0.001389013698677007624551591226760, 0.0]];



pub const WGK2: [[f64;4];8] = [[0.003890461127099884051267201844516,
    0.009273279659517763428441146892024, 0.014369729507045804812451432443580,
    0.019414141193942381173408951050128],
    [0.024191162078080601365686370725232,
    0.028754048765041292843978785354334, 0.032981447057483726031814191016854,
    0.036882364651821229223911065617136],
    [0.040374538951535959111995279752468,
    0.043452539701356069316831728117073, 0.046059238271006988116271735559374,
    0.048185861757087129140779492298305],
    [0.049795683427074206357811569379942,
    0.050881795898749606492297473049805, 0.051426128537459025933862879215781,
    0.051426128537459025933862879215781],
    [0.050881795898749606492297473049805,
    0.049795683427074206357811569379942, 0.048185861757087129140779492298305,
    0.046059238271006988116271735559374],
    [0.043452539701356069316831728117073,
    0.040374538951535959111995279752468, 0.036882364651821229223911065617136,
    0.032981447057483726031814191016854],
    [0.028754048765041292843978785354334,
    0.024191162078080601365686370725232, 0.019414141193942381173408951050128,
    0.014369729507045804812451432443580],
    [0.009273279659517763428441146892024,
    0.003890461127099884051267201844516, 0.000000000000000000000000000000000,
    0.000000000000000000000000000000000]];

pub const WG:  [[f64;4];8] = [[0.007968192496166605615465883474674,0.018466468311090959142302131912047,
    0.028784707883323369349719179611292, 0.038799192569627049596801936446348],
    [0.048402672830594052902938140422808, 0.057493156217619066481721689402056,
    0.065974229882180495128128515115962, 0.073755974737705206268243850022191],
    [0.080755895229420215354694938460530, 0.086899787201082979802387530715126,
    0.092122522237786128717632707087619, 0.096368737174644259639468626351810],
    [0.099593420586795267062780282103569, 0.101762389748405504596428952168554,
    0.102852652893558840341285636705415, 0.102852652893558840341285636705415],
    [0.101762389748405504596428952168554, 0.099593420586795267062780282103569,
    0.096368737174644259639468626351810, 0.092122522237786128717632707087619],
    [0.086899787201082979802387530715126, 0.080755895229420215354694938460530,
    0.073755974737705206268243850022191, 0.065974229882180495128128515115962],
    [0.057493156217619066481721689402056, 0.048402672830594052902938140422808,
    0.038799192569627049596801936446348, 0.028784707883323369349719179611292],
    [0.018466468311090959142302131912047, 0.007968192496166605615465883474674, 0.0, 0.0]];






impl Qk for Qk61Nalgebra {
    fn integrate(&self, f: &dyn Fn(f64) -> f64, a: f64, b: f64, ) -> (f64, f64, f64, f64) {
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let wgk1 = na::vector![f64x4::from_slice_unaligned(&WGK1[0]),f64x4::from_slice_unaligned(&WGK1[1]),
        f64x4::from_slice_unaligned(&WGK1[2]),f64x4::from_slice_unaligned(&WGK1[3]),
            f64x4::from_slice_unaligned(&WGK1[4]),f64x4::from_slice_unaligned(&WGK1[5]),
            f64x4::from_slice_unaligned(&WGK1[6]),f64x4::from_slice_unaligned(&WGK1[7])];

        let wgk2 = na::vector![f64x4::from_slice_unaligned(&WGK2[0]),f64x4::from_slice_unaligned(&WGK2[1]),
        f64x4::from_slice_unaligned(&WGK2[2]),f64x4::from_slice_unaligned(&WGK2[3]),
            f64x4::from_slice_unaligned(&WGK2[4]),f64x4::from_slice_unaligned(&WGK2[5]),
            f64x4::from_slice_unaligned(&WGK2[6]),f64x4::from_slice_unaligned(&WGK2[7])];

        let wg = na::vector![f64x4::from_slice_unaligned(&WG[0]),f64x4::from_slice_unaligned(&WG[1]),
        f64x4::from_slice_unaligned(&WG[2]),f64x4::from_slice_unaligned(&WG[3]),
            f64x4::from_slice_unaligned(&WG[4]),f64x4::from_slice_unaligned(&WG[5]),
            f64x4::from_slice_unaligned(&WG[6]),f64x4::from_slice_unaligned(&WG[7])];

        let (fv1,fv2) = fvec_simd(f, centr, hlgth);
        let resk = fv1.dot(&wgk1).simd_horizontal_sum() + fv2.dot(&wgk2).simd_horizontal_sum();
        let reskh = resk * 0.5;
        let reskhs1 = [reskh;4];
        let reskhs = na::vector![f64x4:: from_slice_unaligned(&reskhs1),f64x4:: from_slice_unaligned(&reskhs1),
        f64x4:: from_slice_unaligned(&reskhs1),f64x4:: from_slice_unaligned(&reskhs1),
        f64x4:: from_slice_unaligned(&reskhs1),f64x4:: from_slice_unaligned(&reskhs1),
        f64x4:: from_slice_unaligned(&reskhs1),f64x4:: from_slice_unaligned(&reskhs1)];
        let mut resabs = wgk1.dot(&fv1.abs()).simd_horizontal_sum() + wgk2.dot(&fv2.abs()).simd_horizontal_sum();
        let resg = fv2.dot(&wg).simd_horizontal_sum();
        let mut resasc = wgk1.dot( &( fv1 - reskhs).abs()).simd_horizontal_sum() +
            wgk2.dot(&( fv2 - reskhs).abs()).simd_horizontal_sum();

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


pub fn fvec_simd(f : &dyn Fn(f64)->f64, centr : f64, hlgth : f64) -> (SVector<f64x4,8>,SVector<f64x4,8>){
    (na::vector![f64x4::new(f(centr + hlgth * XGK[0]), f(centr + hlgth * XGK[2]),
        f(centr + hlgth * XGK[4]), f(centr + hlgth * XGK[6])),
        f64x4::new(f(centr + hlgth * XGK[8]),
        f(centr + hlgth * XGK[10]), f(centr + hlgth * XGK[12]), f(centr + hlgth * XGK[14])),
        f64x4::new(f(centr + hlgth * XGK[16]), f(centr + hlgth * XGK[18]), f(centr + hlgth * XGK[20]),
        f(centr + hlgth * XGK[22])),
        f64x4::new(f(centr + hlgth * XGK[24]),f(centr + hlgth * XGK[26]),
        f(centr + hlgth * XGK[28]), f(centr + hlgth * XGK[30])),
        f64x4::new(f(centr + hlgth * XGK[32]),
        f(centr + hlgth * XGK[34]), f(centr + hlgth * XGK[36]), f(centr + hlgth * XGK[38])),
        f64x4::new(f(centr + hlgth * XGK[40]), f(centr + hlgth * XGK[42]), f(centr + hlgth * XGK[44]),
        f(centr + hlgth * XGK[46])),
        f64x4::new(f(centr + hlgth * XGK[48]), f(centr + hlgth * XGK[50]),
        f(centr + hlgth * XGK[52]), f(centr + hlgth * XGK[54])),
        f64x4::new(f(centr + hlgth * XGK[56]),
        f(centr + hlgth * XGK[58]), f(centr + hlgth * XGK[60]), 0.0)],
     na::vector![f64x4::new(f(centr + hlgth * XGK[1]), f(centr + hlgth * XGK[3]),
         f(centr + hlgth * XGK[5]), f(centr + hlgth * XGK[7])),
         f64x4::new(f(centr + hlgth * XGK[9]),
         f(centr + hlgth * XGK[11]), f(centr + hlgth * XGK[13]), f(centr + hlgth * XGK[15])),
         f64x4::new(f(centr + hlgth * XGK[17]), f(centr + hlgth * XGK[19]), f(centr + hlgth * XGK[21]),
         f(centr + hlgth * XGK[23])),
                 f64x4::new(f(centr + hlgth * XGK[25]), f(centr + hlgth * XGK[27]),
         f(centr + hlgth * XGK[29]), f(centr + hlgth * XGK[31])),
                     f64x4::new(f(centr + hlgth * XGK[33]),
         f(centr + hlgth * XGK[35]), f(centr + hlgth * XGK[37]), f(centr + hlgth * XGK[39])),
         f64x4::new(f(centr + hlgth * XGK[41]), f(centr + hlgth * XGK[43]), f(centr + hlgth * XGK[45]),
         f(centr + hlgth * XGK[47])),
                             f64x4::new(f(centr + hlgth * XGK[49]), f(centr + hlgth * XGK[51]),
         f(centr + hlgth * XGK[53]), f(centr + hlgth * XGK[55])),
                                 f64x4::new(f(centr + hlgth * XGK[57]),
         f(centr + hlgth * XGK[59]), 0.0, 0.0)])

}


#[cfg(test)]
mod tests {
    use std::simd::Simd;
    use std::time::Instant;
    use crate::qk61::Qk61;
    use crate::qk61_1dvec3::Qk611DVec3;
    use crate::qk61_simd2::Qk61Simd2;
    use crate::qk61_simd::Qk61Simd;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.sin();
        let a = 0.0;
        let b = 1.0;
        let qks = Qk61Simd {};
        let qk = Qk61{};
        let qks2 = Qk61Simd2{};

        let mut res2 = (0.0,0.0,0.0,0.0);
        let mut res3 = (0.0,0.0,0.0,0.0);
        let mut res1 = res2.clone();

        for k in 0..100 {
            let start = Instant::now();
            res2 = qk.integrate(&f, a, b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            res1 = qks.integrate(&f, a, b);
            println!("simd {:?}", start.elapsed());
            let start = Instant::now();
            res3 = qks2.integrate(&f, a, b);
            println!("simd2 {:?}", start.elapsed());
        }
        println!("{:?}",res1);
        println!("{:?}",res2);
        println!("{:?}",res3);


    }
}

 */