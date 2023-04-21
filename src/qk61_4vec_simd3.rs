use std::iter::zip;
use packed_simd_2::f64x4;
use crate::funct_vector::FnVec4;
use crate::qk::*;

pub struct Qk61Vec4Simd3 {}
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


pub const WGK1: [f64;31] = [0.001389013698677007624551591226760,
    0.006630703915931292173319826369750, 0.011823015253496341742232898853251,
    0.016920889189053272627572289420322, 0.021828035821609192297167485738339,
    0.026509954882333101610601709335075, 0.030907257562387762472884252943092,
    0.034979338028060024137499670731468, 0.038678945624727592950348651532281,
    0.041969810215164246147147541285970, 0.044814800133162663192355551616723,
    0.047185546569299153945261478181099, 0.049055434555029778887528165367238,
    0.050405921402782346840893085653585, 0.051221547849258772170656282604944,
    0.051494729429451567558340433647099, 0.051221547849258772170656282604944,
    0.050405921402782346840893085653585, 0.049055434555029778887528165367238,
    0.047185546569299153945261478181099, 0.044814800133162663192355551616723,
    0.041969810215164246147147541285970, 0.038678945624727592950348651532281,
    0.034979338028060024137499670731468, 0.030907257562387762472884252943092,
    0.026509954882333101610601709335075, 0.021828035821609192297167485738339,
    0.016920889189053272627572289420322, 0.011823015253496341742232898853251,
    0.006630703915931292173319826369750, 0.001389013698677007624551591226760];

pub const WGK2: [f64;31] = [0.003890461127099884051267201844516,
    0.009273279659517763428441146892024, 0.014369729507045804812451432443580,
    0.019414141193942381173408951050128, 0.024191162078080601365686370725232,
    0.028754048765041292843978785354334, 0.032981447057483726031814191016854,
    0.036882364651821229223911065617136, 0.040374538951535959111995279752468,
    0.043452539701356069316831728117073, 0.046059238271006988116271735559374,
    0.048185861757087129140779492298305, 0.049795683427074206357811569379942,
    0.050881795898749606492297473049805, 0.051426128537459025933862879215781,
    0.051426128537459025933862879215781, 0.050881795898749606492297473049805,
    0.049795683427074206357811569379942, 0.048185861757087129140779492298305,
    0.046059238271006988116271735559374, 0.043452539701356069316831728117073,
    0.040374538951535959111995279752468, 0.036882364651821229223911065617136,
    0.032981447057483726031814191016854, 0.028754048765041292843978785354334,
    0.024191162078080601365686370725232, 0.019414141193942381173408951050128,
    0.014369729507045804812451432443580, 0.009273279659517763428441146892024,
    0.003890461127099884051267201844516, 0.000000000000000000000000000000000];

pub const WG: [f64;31] = [0.007968192496166605615465883474674,0.018466468311090959142302131912047,
    0.028784707883323369349719179611292, 0.038799192569627049596801936446348,
    0.048402672830594052902938140422808, 0.057493156217619066481721689402056,
    0.065974229882180495128128515115962, 0.073755974737705206268243850022191,
    0.080755895229420215354694938460530, 0.086899787201082979802387530715126,
    0.092122522237786128717632707087619, 0.096368737174644259639468626351810,
    0.099593420586795267062780282103569, 0.101762389748405504596428952168554,
    0.102852652893558840341285636705415, 0.102852652893558840341285636705415,
    0.101762389748405504596428952168554, 0.099593420586795267062780282103569,
    0.096368737174644259639468626351810, 0.092122522237786128717632707087619,
    0.086899787201082979802387530715126, 0.080755895229420215354694938460530,
    0.073755974737705206268243850022191, 0.065974229882180495128128515115962,
    0.057493156217619066481721689402056, 0.048402672830594052902938140422808,
    0.038799192569627049596801936446348, 0.028784707883323369349719179611292,
    0.018466468311090959142302131912047, 0.007968192496166605615465883474674, 0.0];


impl Qk61Vec4Simd3 {
    pub(crate) fn integrate(&self, fun: &FnVec4, a: f64, b: f64, ) -> (f64x4,f64x4,f64x4,f64x4 ) {
        let hlgth = f64x4::splat(0.5 * (b - a));
        let dhlgth = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let mut resk = f64x4::splat(0.0);
        let mut resabs = f64x4::splat(0.0);
        let mut resg =  f64x4::splat(0.0);
        let mut resasc =  f64x4::splat(0.0);


        let f = fvec_simd2(fun,centr,hlgth.extract(0));



        for k in 0..31{
            resk += f.0[k] * WGK1[k] + f.1[k] * WGK2[k];
            resabs += f.0[k].abs() * WGK1[k] + f.1[k].abs() * WGK2[k];

        }


/*
        for k in 0..4 {
            let (fv1,fv2) = fvec_simd(&fun.components[k], centr, hlgth[0]);
            resk[k] = (fv1 * WGK1).reduce_sum() + (fv2 * WGK2).reduce_sum();
            let reskh = resk[k] * 0.5;
            let reskhs = Simd::from_array([reskh;32]);
            resabs[k] = (fv1.abs() * WGK1).reduce_sum() + (fv2.abs() * WGK2).reduce_sum();
            resg[k] = ( fv2 * WG).reduce_sum();
            resasc[k] = ( WGK1 * ( fv1 - reskhs).abs()).reduce_sum() +
                ( WGK2 * ( fv2 - reskhs).abs()).reduce_sum();
        }


 */
        let result = resk * hlgth;
        resabs = resabs * dhlgth;
        resasc = resasc * dhlgth;
        let mut abserr = ((resk - resg) * hlgth).abs();
        for k in 0..4 {
            if (resasc.extract(k), abserr.extract(k)) != (0.0, 0.0) {
                abserr.replace(k,resasc.extract(k) * 1.0_f64.min((200.0 * abserr.extract(k) / resasc.extract(k)).powf(1.5)));
            }

            if resabs.extract(k) > UFLOW / (50.0 * EPMACH) {
                abserr.replace(k, abserr.extract(k).max((EPMACH * 50.0) * resabs.extract(k)));
            }
        }

        (result, abserr, resabs, resasc)
    }
}
pub fn fvec_simd2(fun: &FnVec4, centr : f64, hlgth : f64) -> ([f64x4;31],[f64x4;31]){
    ([f64x4::new(fun.components[0](centr + hlgth * XGK[0]),
                            fun.components[1](centr + hlgth * XGK[0]), fun.components[2](centr + hlgth * XGK[0]),
                            fun.components[3](centr + hlgth * XGK[0])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[2]),
               fun.components[1](centr + hlgth * XGK[2]), fun.components[2](centr + hlgth * XGK[2]),
               fun.components[3](centr + hlgth * XGK[2])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[4]),
               fun.components[1](centr + hlgth * XGK[4]), fun.components[2](centr + hlgth * XGK[4]),
               fun.components[3](centr + hlgth * XGK[4])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[6]),
               fun.components[1](centr + hlgth * XGK[6]), fun.components[2](centr + hlgth * XGK[6]),
               fun.components[3](centr + hlgth * XGK[6])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[8]),
               fun.components[1](centr + hlgth * XGK[8]), fun.components[2](centr + hlgth * XGK[8]),
               fun.components[3](centr + hlgth * XGK[8])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[10]),
               fun.components[1](centr + hlgth * XGK[10]), fun.components[2](centr + hlgth * XGK[10]),
               fun.components[3](centr + hlgth * XGK[10])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[12]),
               fun.components[1](centr + hlgth * XGK[12]), fun.components[2](centr + hlgth * XGK[12]),
               fun.components[3](centr + hlgth * XGK[12])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[14]),
               fun.components[1](centr + hlgth * XGK[14]), fun.components[2](centr + hlgth * XGK[14]),
               fun.components[3](centr + hlgth * XGK[14])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[16]),
               fun.components[1](centr + hlgth * XGK[16]), fun.components[2](centr + hlgth * XGK[16]),
               fun.components[3](centr + hlgth * XGK[16])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[18]),
               fun.components[1](centr + hlgth * XGK[18]), fun.components[2](centr + hlgth * XGK[18]),
               fun.components[3](centr + hlgth * XGK[18])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[20]),
               fun.components[1](centr + hlgth * XGK[20]), fun.components[2](centr + hlgth * XGK[20]),
               fun.components[3](centr + hlgth * XGK[20])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[22]),
               fun.components[1](centr + hlgth * XGK[22]), fun.components[2](centr + hlgth * XGK[22]),
               fun.components[3](centr + hlgth * XGK[22])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[24]),
               fun.components[1](centr + hlgth * XGK[24]), fun.components[2](centr + hlgth * XGK[24]),
               fun.components[3](centr + hlgth * XGK[24])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[26]),
               fun.components[1](centr + hlgth * XGK[26]), fun.components[2](centr + hlgth * XGK[26]),
               fun.components[3](centr + hlgth * XGK[26])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[28]),
               fun.components[1](centr + hlgth * XGK[28]), fun.components[2](centr + hlgth * XGK[28]),
               fun.components[3](centr + hlgth * XGK[28])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[30]),
               fun.components[1](centr + hlgth * XGK[30]), fun.components[2](centr + hlgth * XGK[30]),
               fun.components[3](centr + hlgth * XGK[30])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[32]),
               fun.components[1](centr + hlgth * XGK[32]), fun.components[2](centr + hlgth * XGK[32]),
               fun.components[3](centr + hlgth * XGK[32])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[34]),
               fun.components[1](centr + hlgth * XGK[34]), fun.components[2](centr + hlgth * XGK[34]),
               fun.components[3](centr + hlgth * XGK[34])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[36]),
               fun.components[1](centr + hlgth * XGK[36]), fun.components[2](centr + hlgth * XGK[36]),
               fun.components[3](centr + hlgth * XGK[36])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[38]),
               fun.components[1](centr + hlgth * XGK[38]), fun.components[2](centr + hlgth * XGK[38]),
               fun.components[3](centr + hlgth * XGK[38])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[40]),
               fun.components[1](centr + hlgth * XGK[40]), fun.components[2](centr + hlgth * XGK[40]),
               fun.components[3](centr + hlgth * XGK[40])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[42]),
               fun.components[1](centr + hlgth * XGK[42]), fun.components[2](centr + hlgth * XGK[42]),
               fun.components[3](centr + hlgth * XGK[42])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[44]),
               fun.components[1](centr + hlgth * XGK[44]), fun.components[2](centr + hlgth * XGK[44]),
               fun.components[3](centr + hlgth * XGK[44])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[46]),
               fun.components[1](centr + hlgth * XGK[46]), fun.components[2](centr + hlgth * XGK[46]),
               fun.components[3](centr + hlgth * XGK[46])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[48]),
               fun.components[1](centr + hlgth * XGK[48]), fun.components[2](centr + hlgth * XGK[48]),
               fun.components[3](centr + hlgth * XGK[48])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[50]),
               fun.components[1](centr + hlgth * XGK[50]), fun.components[2](centr + hlgth * XGK[50]),
               fun.components[3](centr + hlgth * XGK[50])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[52]),
               fun.components[1](centr + hlgth * XGK[52]), fun.components[2](centr + hlgth * XGK[52]),
               fun.components[3](centr + hlgth * XGK[52])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[54]),
               fun.components[1](centr + hlgth * XGK[54]), fun.components[2](centr + hlgth * XGK[54]),
               fun.components[3](centr + hlgth * XGK[54])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[56]),
               fun.components[1](centr + hlgth * XGK[56]), fun.components[2](centr + hlgth * XGK[56]),
               fun.components[3](centr + hlgth * XGK[56])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[58]),
               fun.components[1](centr + hlgth * XGK[58]), fun.components[2](centr + hlgth * XGK[58]),
               fun.components[3](centr + hlgth * XGK[58])),
    f64x4::new(fun.components[0](centr + hlgth * XGK[60]),
               fun.components[1](centr + hlgth * XGK[60]), fun.components[2](centr + hlgth * XGK[60]),
               fun.components[3](centr + hlgth * XGK[60]))],




    [f64x4::new(fun.components[0](centr + hlgth * XGK[1]),
                fun.components[1](centr + hlgth * XGK[1]), fun.components[2](centr + hlgth * XGK[1]),
                fun.components[3](centr + hlgth * XGK[1])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[3]),
                   fun.components[1](centr + hlgth * XGK[3]), fun.components[2](centr + hlgth * XGK[3]),
                   fun.components[3](centr + hlgth * XGK[3])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[5]),
                   fun.components[1](centr + hlgth * XGK[5]), fun.components[2](centr + hlgth * XGK[5]),
                   fun.components[3](centr + hlgth * XGK[5])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[7]),
                   fun.components[1](centr + hlgth * XGK[7]), fun.components[2](centr + hlgth * XGK[7]),
                   fun.components[3](centr + hlgth * XGK[7])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[9]),
                   fun.components[1](centr + hlgth * XGK[9]), fun.components[2](centr + hlgth * XGK[9]),
                   fun.components[3](centr + hlgth * XGK[9])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[11]),
                   fun.components[1](centr + hlgth * XGK[11]), fun.components[2](centr + hlgth * XGK[11]),
                   fun.components[3](centr + hlgth * XGK[11])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[13]),
                   fun.components[1](centr + hlgth * XGK[13]), fun.components[2](centr + hlgth * XGK[13]),
                   fun.components[3](centr + hlgth * XGK[13])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[15]),
                   fun.components[1](centr + hlgth * XGK[15]), fun.components[2](centr + hlgth * XGK[15]),
                   fun.components[3](centr + hlgth * XGK[15])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[17]),
                   fun.components[1](centr + hlgth * XGK[17]), fun.components[2](centr + hlgth * XGK[17]),
                   fun.components[3](centr + hlgth * XGK[17])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[19]),
                   fun.components[1](centr + hlgth * XGK[19]), fun.components[2](centr + hlgth * XGK[19]),
                   fun.components[3](centr + hlgth * XGK[19])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[21]),
                   fun.components[1](centr + hlgth * XGK[21]), fun.components[2](centr + hlgth * XGK[21]),
                   fun.components[3](centr + hlgth * XGK[21])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[23]),
                   fun.components[1](centr + hlgth * XGK[23]), fun.components[2](centr + hlgth * XGK[23]),
                   fun.components[3](centr + hlgth * XGK[23])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[25]),
                   fun.components[1](centr + hlgth * XGK[25]), fun.components[2](centr + hlgth * XGK[25]),
                   fun.components[3](centr + hlgth * XGK[25])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[27]),
                   fun.components[1](centr + hlgth * XGK[27]), fun.components[2](centr + hlgth * XGK[27]),
                   fun.components[3](centr + hlgth * XGK[27])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[29]),
                   fun.components[1](centr + hlgth * XGK[29]), fun.components[2](centr + hlgth * XGK[29]),
                   fun.components[3](centr + hlgth * XGK[29])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[31]),
                   fun.components[1](centr + hlgth * XGK[31]), fun.components[2](centr + hlgth * XGK[31]),
                   fun.components[3](centr + hlgth * XGK[31])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[33]),
                   fun.components[1](centr + hlgth * XGK[33]), fun.components[2](centr + hlgth * XGK[33]),
                   fun.components[3](centr + hlgth * XGK[33])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[35]),
                   fun.components[1](centr + hlgth * XGK[35]), fun.components[2](centr + hlgth * XGK[35]),
                   fun.components[3](centr + hlgth * XGK[35])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[37]),
                   fun.components[1](centr + hlgth * XGK[37]), fun.components[2](centr + hlgth * XGK[37]),
                   fun.components[3](centr + hlgth * XGK[37])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[39]),
                   fun.components[1](centr + hlgth * XGK[39]), fun.components[2](centr + hlgth * XGK[39]),
                   fun.components[3](centr + hlgth * XGK[39])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[41]),
                   fun.components[1](centr + hlgth * XGK[41]), fun.components[2](centr + hlgth * XGK[41]),
                   fun.components[3](centr + hlgth * XGK[41])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[43]),
                   fun.components[1](centr + hlgth * XGK[43]), fun.components[2](centr + hlgth * XGK[43]),
                   fun.components[3](centr + hlgth * XGK[43])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[45]),
                   fun.components[1](centr + hlgth * XGK[45]), fun.components[2](centr + hlgth * XGK[45]),
                   fun.components[3](centr + hlgth * XGK[45])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[47]),
                   fun.components[1](centr + hlgth * XGK[47]), fun.components[2](centr + hlgth * XGK[47]),
                   fun.components[3](centr + hlgth * XGK[47])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[49]),
                   fun.components[1](centr + hlgth * XGK[49]), fun.components[2](centr + hlgth * XGK[49]),
                   fun.components[3](centr + hlgth * XGK[49])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[51]),
                   fun.components[1](centr + hlgth * XGK[51]), fun.components[2](centr + hlgth * XGK[51]),
                   fun.components[3](centr + hlgth * XGK[51])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[53]),
                   fun.components[1](centr + hlgth * XGK[53]), fun.components[2](centr + hlgth * XGK[53]),
                   fun.components[3](centr + hlgth * XGK[53])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[55]),
                   fun.components[1](centr + hlgth * XGK[55]), fun.components[2](centr + hlgth * XGK[55]),
                   fun.components[3](centr + hlgth * XGK[55])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[57]),
                   fun.components[1](centr + hlgth * XGK[57]), fun.components[2](centr + hlgth * XGK[57]),
                   fun.components[3](centr + hlgth * XGK[57])),
        f64x4::new(fun.components[0](centr + hlgth * XGK[59]),
                   fun.components[1](centr + hlgth * XGK[59]), fun.components[2](centr + hlgth * XGK[59]),
                   fun.components[3](centr + hlgth * XGK[59])),
        f64x4::splat(0.0)])


}


/*
pub fn fvec_simd2(fun: &FnVec4, centr : f64, hlgth : f64) -> Vec<f64x4>{
    let mut vec = vec![];
    for x in XGK {
        vec.push(f64x4::new(fun.components[0](centr + hlgth * x),
            fun.components[1](centr + hlgth * x), fun.components[2](centr + hlgth * x),
            fun.components[3](centr + hlgth * x)));
    }
    vec.push(f64x4::splat(0.0));
    vec

}


 */
#[cfg(test)]
mod tests {
    use std::simd::f64x4;
    use std::time::Instant;
    use crate::funct_vector::FnVec4;
    use crate::qk61::Qk61;
    use crate::qk61_1dvec3::Qk611DVec3;
    use crate::qk61_4vec_simd2::Qk61Vec4Simd2;
    use crate::qk61_4vec_simd3::Qk61Vec4Simd3;
    use crate::qk61_simd::Qk61Simd;
    use crate::qk61_4vec_simd::Qk61Vec4Simd;
    use crate::qk61_simd2::Qk61Simd2;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.cos();

        let a = 0.0;
        let b = 1.0;
        let qks = Qk61Vec4Simd{};
        let qks2 = Qk61Vec4Simd3{};
        let qk = Qk61Simd {};
        let qk2 = Qk61Simd2 {};
        let fun = FnVec4{ components : [Box::new(f),Box::new(f),Box::new(f),Box::new(f)]};

        let null = f64x4::splat(0.0);
        let null2 = packed_simd_2::f64x4::splat(0.0);
        let mut res_simd = (null,null,null,null);
        let mut res_simd2 = (null2,null2,null2,null2);

        for k in 0..100 {
            let start = Instant::now();
            let res1 = qk.integrate(&f, a, b);
            let res2 = qk.integrate(&f,a,b);
            let res3 = qk.integrate(&f,a,b);
            let res4 = qk.integrate(&f,a,b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            let res1 = qk2.integrate(&f, a, b);
            let res2 = qk2.integrate(&f,a,b);
            let res3 = qk2.integrate(&f,a,b);
            let res4 = qk2.integrate(&f,a,b);
            println!("normal2 {:?}", start.elapsed());
            let start = Instant::now();
            res_simd = qks.integrate(&fun, a, b);
            println!("simd {:?}", start.elapsed());
            let start = Instant::now();
            res_simd2 = qks2.integrate(&fun, a, b);
            println!("simd2 {:?}", start.elapsed());
        }
        println!("simd {:?}",res_simd);
        println!("simd2 {:?}",res_simd2);
    }
}


