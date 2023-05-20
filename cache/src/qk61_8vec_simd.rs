use std::simd::{f64x8, Simd, SimdFloat};
use crate::funct_vector::FnVec8;
use crate::qk::*;
use crate::qk61_simd2::*;

pub struct Qk61Vec8Simd {}

impl Qk61Vec8Simd {
    pub(crate) fn integrate(&self, fun: &FnVec8, a: f64, b: f64, ) -> (Simd<f64, 8>, Simd<f64, 8>, Simd<f64, 8>, Simd<f64, 8>, ) {
        let hlgth = f64x8::splat(0.5 * (b - a));
        let dhlgth = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let mut resk = f64x8::splat(0.0);
        let mut resabs = f64x8::splat(0.0);
        let mut resg =  f64x8::splat(0.0);
        let mut resasc =  f64x8::splat(0.0);


        for k in 0..8 {
            let (fv1,fv2) = fvec_simd(&fun.components[k], centr, hlgth[0]);
            resk[k] = (fv1 * WGK1).reduce_sum() + (fv2 * WGK2).reduce_sum();
            let reskh = resk[k] * 0.5;
            let reskhs = Simd::from_array([reskh;32]);
            resabs[k] = (fv1.abs() * WGK1).reduce_sum() + (fv2.abs() * WGK2).reduce_sum();
            resg[k] = ( fv2 * WG).reduce_sum();
            resasc[k] = ( WGK1 * ( fv1 - reskhs).abs()).reduce_sum() +
                ( WGK2 * ( fv2 - reskhs).abs()).reduce_sum();
        }

        let result = resk * hlgth;
        resabs = resabs * dhlgth;
        resasc = resasc * dhlgth;
        let mut abserr = ((resk - resg) * hlgth).abs();
        for k in 0..4 {
            if (resasc[k], abserr[k]) != (0.0, 0.0) {
                abserr[k] = resasc[k] * 1.0_f64.min((200.0 * abserr[k] / resasc[k]).powf(1.5));
            }

            if resabs[k] > UFLOW / (50.0 * EPMACH) {
                abserr[k] = abserr[k].max((EPMACH * 50.0) * resabs[k]);
            }
        }

        (result, abserr, resabs, resasc)
    }
}

#[cfg(test)]
mod tests {
    use std::simd::{f64x4, f64x8};
    use std::time::Instant;
    use crate::funct_vector::{FnVec4, FnVec8};
    use crate::qk61::Qk61;
    use crate::qk61_1dvec3::Qk611DVec3;
    use crate::qk61_4vec_simd2::Qk61Vec4Simd2;
    use crate::qk61_simd::Qk61Simd;
    use crate::qk61_4vec_simd::Qk61Vec4Simd;
    use crate::qk61_8vec_simd::Qk61Vec8Simd;
    use crate::qk61_simd2::Qk61Simd2;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.cos();

        let a = 0.0;
        let b = 1.0;
        let qks2 = Qk61Vec8Simd{};
        let qk2 = Qk61Simd2 {};
        let fun = FnVec8{ components : [Box::new(f),Box::new(f),Box::new(f),Box::new(f),
            Box::new(f),Box::new(f),Box::new(f),Box::new(f)]};

        let null = f64x8::splat(0.0);
        let mut res1 = (0.0,0.0,0.0,0.0);
        let mut res_simd2 = (null,null,null,null);

        for k in 0..100 {
            let start = Instant::now();
            res1 = qk2.integrate(&f, a, b);
            let res2 = qk2.integrate(&f,a,b);
            let res3 = qk2.integrate(&f,a,b);
            let res4 = qk2.integrate(&f,a,b);
            let res1 = qk2.integrate(&f, a, b);
            let res2 = qk2.integrate(&f,a,b);
            let res3 = qk2.integrate(&f,a,b);
            let res4 = qk2.integrate(&f,a,b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            res_simd2 = qks2.integrate(&fun, a, b);
            println!("simd {:?}", start.elapsed());
        }
        println!("normal {:?}",res1);
        println!("simd2 {:?}",res_simd2);
    }
}

/*

pub(crate) fn integrate(&self, fun: &FnVec4, a: f64, b: f64, ) -> (Simd<f64, 4>, Simd<f64, 4>, Simd<f64, 4>, Simd<f64, 4>, ) {
    let qk61_1d = Qk61Simd2 {};
    /*
    let (mut result_vec,mut abserr_vec,mut resabs_vec,mut resasc_vec) =
        ([0.0;4],[0.0;4],[0.0;4],[0.0;4]);
    for k in 0..4 {
        (result_vec[k], abserr_vec[k], resabs_vec[k], resasc_vec[k]) = qk61_1d.integrate(&fun.components[k],a,b);
    }
    (Simd::from_array(result_vec),Simd::from_array(abserr_vec),
     Simd::from_array(resabs_vec),Simd::from_array(resasc_vec))
      */

    let (mut result_vec, mut abserr_vec, mut resabs_vec, mut resasc_vec) =
        (f64x4::splat(0.0), f64x4::splat(0.0), f64x4::splat(0.0), f64x4::splat(0.0));
    (result_vec[0], abserr_vec[0], resabs_vec[0], resasc_vec[0]) =
        qk61_1d.integrate(&fun.components[0], a, b);
    (result_vec[1], abserr_vec[1], resabs_vec[1], resasc_vec[1]) =
        qk61_1d.integrate(&fun.components[1], a, b);
    (result_vec[2], abserr_vec[2], resabs_vec[2], resasc_vec[2]) =
        qk61_1d.integrate(&fun.components[2], a, b);
    (result_vec[3], abserr_vec[3], resabs_vec[3], resasc_vec[3]) =
        qk61_1d.integrate(&fun.components[3], a, b);
    (result_vec, abserr_vec, resabs_vec, resasc_vec)
}

 */
