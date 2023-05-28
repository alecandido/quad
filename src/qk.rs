use crate::constants::*;

pub fn qk_quadrature<const N:usize,const M : usize,F>(f: F, a: f64, b: f64, xgk : &[f64;M], wgk : &[f64], wg : &[f64])
                                                      -> ([f64; N], f64, f64)
where F : Fn(f64) -> [f64; N]{
    let hlgth: f64 = 0.5 * (b - a);
    let dhlgth: f64 = hlgth.abs();
    let centr: f64 = 0.5 * (b + a);

    let mut fv1 = [[0.0; N]; M];
    let mut fv2 = [[0.0; N]; M];


    let fc : [f64; N]  = f(centr);
    let mut resg = [0.0; N];
    let mut resk = fc.clone();
    let mut resabs = [0.0;N];
    if M % 2 == 1{
        resg = fc.clone();
    }
    for k in 0..N {
        resk[k] *= wgk[M];
        resabs[k] = resk[k].abs();
        if M % 2 == 1 {
            resg[k] *= wg[(M+1)/2-1];
        }
    }


    for j in 1..M/2 + 1 {
        let jtw = 2 * j;
        let absc = hlgth * xgk[jtw - 1];
        let fval1 : [f64; N] = f(centr - absc);
        let fval2 : [f64; N] = f(centr + absc);
        fv1[jtw - 1] = fval1;
        fv2[jtw - 1] = fval2;
        let mut fsum : [f64; N] = [0.0; N];
        for k in 0..N {
            fsum[k] = fval1[k] + fval2[k];
        }
        for k in 0..N {
            resg[k] += wg[j - 1] * fsum[k];
            resk[k] += wgk[jtw - 1] * fsum[k];
            resabs[k] += wgk[jtw - 1] * (fval1[k].abs() + fval2[k].abs());
        }

    }

    for j in 1..(M+1)/2 + 1 {
        let jtwm1 = 2 * j - 1;
        let absc = hlgth * xgk[jtwm1 - 1];
        let fval1 : [f64; N] = f(centr - absc);
        let fval2 : [f64; N] = f(centr + absc);
        fv1[jtwm1 - 1] = fval1;
        fv2[jtwm1 - 1] = fval2;
        let mut fsum : [f64; N] = [0.0; N];
        for k in 0..N {
            fsum[k] = fval1[k] + fval2[k];
        }
        for k in 0..N {
            resk[k] += wgk[jtwm1 - 1] * fsum[k];
            resabs[k] += wgk[jtwm1 - 1] * (fval1[k].abs() + fval2[k].abs());
        }
    }

    let mut reskh =resk.clone();
    for k  in 0..N {
        reskh[k] *= 0.5;
    }

    let mut resasc = [0.0; N];
    for k in 0..N {
        resasc[k] = wgk[M] * ( fc[k] - reskh[k] ).abs();
    }

    for j in 1..M + 1 {
        for k in 0..N {
            resasc[k] += wgk[j - 1] * ((fv1[j - 1][k] - reskh[k]).abs() + (fv2[j - 1][k] -
                reskh[k]).abs());
        }

    }

    let mut  result = resk.clone();

    for k in 0..N {
        result[k] *= hlgth;
        resabs[k] *= dhlgth;
        resasc[k] *= dhlgth;
    }

    let mut abserr = 0.0;
    let mut resabs_scalar = 0.0;
    let mut resasc_scalar = 0.0;

    for k in 0..N {
        abserr +=  (((resk[k] - resg[k]) * hlgth).abs()).powi(2);
        resabs_scalar += resabs[k].powi(2);
        resasc_scalar += resasc[k].powi(2);
    }

    abserr = abserr.sqrt();
    resabs_scalar = resabs_scalar.sqrt();
    resasc_scalar = resasc_scalar.sqrt();


    if resasc_scalar != 0.0 && abserr != 0.0 {
        abserr = resasc_scalar * 1.0_f64.min((200.0 * abserr / resasc_scalar).powf(1.5));
    }

    let round_error = 50.0 * EPMACH * resabs_scalar;

    if round_error > UFLOW {
        abserr = abserr.max(round_error);
    }


    (result, abserr, round_error)



}








#[cfg(test)]
mod tests {
    use crate::qk61_vec_norm2::*;
    use crate::qk::qk_quadrature;

    #[test]
    fn test(){
        let f = |x:f64| [x.cos(),x.sin()];
        let qk = Qk61VecNorm2{};

        let res = qk_quadrature(&f,0.0,1.0,&XGK,&WGK,&WG);
        let res2 = qk.integrate(&f,0.0,1.0);

        println!("{:?}",res);
        println!("{:?}",res2);

    }
}