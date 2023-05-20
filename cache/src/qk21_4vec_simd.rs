use std::simd::{f64x4, Simd, SimdFloat};
use crate::funct_vector::FnVec4;
use crate::qk::*;
use crate::qk21_simd::*;

pub struct Qk21Vec4Simd {}

impl Qk21Vec4Simd {
    pub(crate) fn integrate(&self, fun: &FnVec4, a: f64, b: f64, ) -> (Simd<f64, 4>, Simd<f64, 4>, Simd<f64, 4>, Simd<f64, 4>, ) {
        let qk21 = Qk21Simd{};
        let mut result = [0.0;4];
        let mut abserr = [0.0;4];
        let mut resabs = [0.0;4];
        let mut resasc = [0.0;4];

        for k in 0..4{
            (result[k],abserr[k],resabs[k],resasc[k]) = qk21.integrate(&fun.components[k],a,b);
        }

        (Simd::from_array(result), Simd::from_array(abserr),
         Simd::from_array(resabs), Simd::from_array(resasc))
    }
    pub(crate) fn integrate2(&self, fun: &FnVec4, a: f64, b: f64, ) -> (Simd<f64, 4>, Simd<f64, 4>, Simd<f64, 4>, Simd<f64, 4>, ) {
        let hlgth = f64x4::splat(0.5 * (b - a));
        let dhlgth = hlgth.abs();
        let centr: f64 = 0.5 * (b + a);

        let mut resk = f64x4::splat(0.0);
        let mut resabs = f64x4::splat(0.0);
        let mut resg =  f64x4::splat(0.0);
        let mut resasc =  f64x4::splat(0.0);


        for k in 0..4 {
            let fv = fvec_simd(&fun.components[k], centr, hlgth[k]);
            resk[k] = (fv * WGK).reduce_sum();
            let reskh = resk[k] * 0.5;
            let reskhs = Simd::from_array([reskh; 32]);
            resabs[k] = (fv.abs() * WGK).reduce_sum();
            resg[k] = (fv * WG).reduce_sum();
            resasc[k] = (WGK * (fv - reskhs).abs()).reduce_sum();
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
    use std::time::Instant;
    use crate::funct_vector::FnVec4;
    use crate::qk21_4vec_simd::Qk21Vec4Simd;
    use crate::qk21_simd::Qk21Simd;
    use crate::qk41_simd::Qk41Simd;
    use crate::qk61::Qk61;
    use crate::qk61_1dvec3::Qk611DVec3;
    use crate::qk61_simd::Qk61Simd;
    use crate::qk61_4vec_simd::Qk61Vec4Simd;
    use crate::qk::Qk;

    #[test]
    fn test(){
        let f = |x:f64| x.cos();

        let a = 0.0;
        let b = 1.0;
        let qks = Qk21Vec4Simd {};
        let qk = Qk21Simd {};
        let fun = FnVec4{ components : [Box::new(f),Box::new(f),Box::new(f),Box::new(f)]};

        for k in 0..100 {
            let start = Instant::now();
            let res1 = qk.integrate(&f, a, b);
            let res2 = qk.integrate(&f,a,b);
            let res3 = qk.integrate(&f,a,b);
            let res4 = qk.integrate(&f,a,b);
            println!("normal {:?}", start.elapsed());
            let start = Instant::now();
            let res_simd = qks.integrate(&fun, a, b);
            println!("simd {:?}", start.elapsed());
            let start = Instant::now();
            let res_simd2 = qks.integrate2(&fun, a, b);
            println!("simd2 {:?}", start.elapsed());
        }
        //println!("{:?}",res1);
        //println!("{:?}",res1);
    }
}
