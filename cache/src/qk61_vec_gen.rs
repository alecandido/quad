use std::iter::zip;
use std::time::Instant;
use crate::qk::*;

pub struct Qk61VecGen {}
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

const XGK : [f64;61] = [-0.999484410050490637571325895705811, -0.996893484074649540271630050918695,
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

const WGK : [f64;61] = [0.001389013698677007624551591226760, 0.003890461127099884051267201844516,
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
    0.001389013698677007624551591226760];

const WG :  [f64;61] = [0.0, 0.007968192496166605615465883474674, 0.0, 0.018466468311090959142302131912047,
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
    0.0, 0.007968192496166605615465883474674, 0.0];

const WG2 :  [f64;31] = [0.007968192496166605615465883474674, 0.018466468311090959142302131912047,
    0.028784707883323369349719179611292,  0.038799192569627049596801936446348,
    0.048402672830594052902938140422808,  0.057493156217619066481721689402056,
    0.065974229882180495128128515115962,  0.073755974737705206268243850022191,
    0.080755895229420215354694938460530,  0.086899787201082979802387530715126,
    0.092122522237786128717632707087619,  0.096368737174644259639468626351810,
    0.099593420586795267062780282103569,  0.101762389748405504596428952168554,
    0.102852652893558840341285636705415,  0.0,
    0.102852652893558840341285636705415,  0.101762389748405504596428952168554,
    0.099593420586795267062780282103569,  0.096368737174644259639468626351810,
    0.092122522237786128717632707087619,  0.086899787201082979802387530715126,
    0.080755895229420215354694938460530,  0.073755974737705206268243850022191,
    0.065974229882180495128128515115962,  0.057493156217619066481721689402056,
    0.048402672830594052902938140422808,  0.038799192569627049596801936446348,
    0.028784707883323369349719179611292,  0.018466468311090959142302131912047,
    0.007968192496166605615465883474674];



impl Qk61VecGen {
    pub(crate) fn integrate(&self, f: &dyn Fn(f64) -> Vec<f64>, a: f64, b: f64, ) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {

        //  let start = Instant::now();
        let hlgth: f64 = 0.5 * (b - a);
        let dhlgth: f64 = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        //  println!("init var : {:?}",start.elapsed());

        //let fv = XGK.map( |x| f(centr + hlgth * x) );
        let fv = XGK.map( |x| f(centr + hlgth * x) );

        //  println!("init f : {:?}",start.elapsed());
        //let fv = fvec(f,centr,hlgth);

        let (mut resk, mut resabs,resg) = scalar_products3(&fv);

        //  println!("scalar prod 1/3 : {:?}",start.elapsed());

        //let reskh = resk * 0.5;

        let mut resasc = scalar_product_diff_abs(&fv,&resk,&WGK);

        //  println!("scalar prod 4 : {:?}",start.elapsed());



        //  println!("rescale and norm : {:?}",start.elapsed());

        let mut abserr = err(&resk,&resg,dhlgth);
        rescale_vec(&mut resk,hlgth);


        for k in 0..abserr.len(){
            if resasc[k] != 0.0 && abserr[k] != 0.0 {
                abserr[k] = resasc[k] * 1.0_f64.min((200.0 * abserr[k] / resasc[k]).powf(1.5));
            }
            if resabs[k] > UFLOW / (50.0 * EPMACH) {
                abserr[k] = abserr[k].max((EPMACH * 50.0) * resabs[k]);
            }
        }

        //  println!("error : {:?}",start.elapsed());

        (resk, abserr, resabs, resasc)
    }
}


pub fn fvec(f : &dyn Fn(f64)->f64, centr : f64, hlgth : f64) -> [f64;61]{
    [f(centr + hlgth * XGK[0]), f(centr + hlgth * XGK[1]), f(centr + hlgth * XGK[2]),
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
        f(centr + hlgth * XGK[60])]
}




pub fn scalar_products3(w: &[Vec<f64>]) -> (Vec<f64>,Vec<f64>,Vec<f64>){
    let mut res1 = vec![];
    let mut res2 = vec![];
    let mut res3 = vec![];
    for k in 0..w[0].len(){
        let (mut sum, mut sumneg,mut sum2) = (0.0,0.0,0.0);
        for (i,j) in zip(zip(&WGK, w),&WG){
            sum += i.0 * i.1[k];
            if i.0 * i.1[k] < 0.0 {
                sumneg += i.0 * i.1[k];
            }
            if *j != 0.0 {
                sum2 += i.1[k] * j;
            }
        }
        res1.push(sum);
        res2.push(sum - 2.0 * sumneg);
        res3.push(sum2);
    }
    (res1,res2,res3)
}





pub fn scalar_product_diff_abs(v: &[Vec<f64>], w: &[f64], z: &[f64]) -> Vec<f64>{
    let mut res = vec![];
    for i in 0..v[0].len(){
        let mut sum = 0.0;
        for k in 0..v.len(){
            sum += z[k] * ( v[k][i] - 0.5 * w[i]).abs();
        }
        res.push(sum)
    }
    res

}

pub fn rescale_vec(v : &mut [f64], s : f64){
    for comp in v{
        *comp = *comp * s;
    }
}

pub fn err(v : &[f64], w : &[f64] , s : f64) -> Vec<f64>{
    let mut err = vec![];
    for (i,j) in zip(v,w){
        err.push((i - j).abs() * s) ;
    }
    err
}

pub fn norm_s(v : &[f64], s : f64) -> f64{
    let mut norm = 0.0;
    for comp in v{
        norm += (comp * s).powi(2);
    }
    norm = norm.sqrt();
    norm
}

#[cfg(test)]
mod tests {
    use std::time::Instant;
    use crate::qk61::Qk61;
    use crate::qk61_4vec::Qk614Vec;
    use crate::qk61_vec_gen::Qk61VecGen;
    use crate::qk61_vec_norm::Qk61VecNorm;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f1 = |x:f64| x.cos();
        let f2 = |x:f64| x.sin();
        let f3 = |x:f64| f1(x) + f2(x);
        let f4 = |x:f64| - 2.0 * f3(x);

        let a = 0.0;
        let b = 1000.0;
        let max = 2;
        let qk = Qk61 {};
        let qk_vec = Qk614Vec {};
        let qk_vec_gen_norm = Qk61VecNorm {};
        let qk_vec_gen = Qk61VecGen{};
        let f = |x:f64| [f1(x),f2(x),f3(x),f4(x)];
        let fvec  = |x:f64| vec![f1(x),f2(x),f3(x),f4(x)];


        let mut res1 = (0.0,0.0,0.0,0.0);
        let mut res2 = res1.clone();
        let mut res3 = res1.clone();
        let mut res4 = res1.clone();
        let mut res_vec = ([0.0;4],[0.0;4],[0.0;4],[0.0;4]);
        let mut res_vec_gen_norm = (vec![], 0.0, 0.0);
        let mut res_vec_gen = (vec![],vec![],vec![],vec![]);

        for k in 0..max {
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
            res_vec_gen_norm = qk_vec_gen_norm.integrate(&fvec, a, b);
            println!("vec gen norm {:?}", start.elapsed());
            let start = Instant::now();
            res_vec_gen = qk_vec_gen.integrate(&fvec, a, b);
            println!("vec gen {:?}", start.elapsed());

            if k == max-1 {
                println!("normal {:?},{:?},{:?},{:?}",res1,res2,res3,res4);
                println!("vec {:?}",res_vec);
                println!("vec gen norm {:?}", res_vec_gen_norm);
                println!("vec gen {:?}", res_vec_gen);
            }
        }

    }
}


