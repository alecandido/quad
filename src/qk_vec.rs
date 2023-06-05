//use std::time::Instant;
use crate::constants::*;

pub fn qk_quadrature_vec<const M : usize,F>(f: F, a: f64, b: f64, xgk : &[f64;M], wgk : &[f64], wg : &[f64])
                                            -> (Vec<f64>, f64, f64)
    where F : Fn(f64) -> Vec<f64>{
    let hlgth: f64 = 0.5 * (b - a);
    let dhlgth: f64 = hlgth.abs();
    let centr: f64 = 0.5 * (b + a);

    //let start = Instant::now();

    let fc = f(centr);
    let dim = fc.len();
    let mut fv1 = Vec::with_capacity(dim*M);
    let mut fv2 = Vec::with_capacity(dim*M);
    let mut resg : Vec<f64> = Vec::with_capacity(dim);
    let mut resk : Vec<f64> = fc.clone();
    let mut resabs = Vec::with_capacity(dim);

    //println!("First initialization : {:?}",start.elapsed());

    if M % 2 == 1 {
        resg.extend( &fc);
    }
    else {
        resg.extend(& vec![0.0; dim]);
    }
    for k in 0..dim {
        //resk.push( fc[k] * wgk[M]);
        resk[k] *= wgk[M];
        resabs.push(resk[k].abs());
        if M % 2 == 1 {
            resg[k] *= wg[(M+1)/2-1];
        }
    }

    //println!("Second initialization : {:?}",start.elapsed());

    for j in 1..M/2 + 1 {
        let jtw1 = 2 * j - 1;
        let jtw2 = 2 * j;

        let absc1 = hlgth * xgk[jtw1 - 1];
        let absc2 = hlgth * xgk[jtw2 - 1];

        fv1.append(&mut f(centr - absc1));
        fv1.append(&mut f(centr - absc2));
        fv2.append(&mut f(centr + absc1));
        fv2.append(&mut f(centr + absc2));

        //  let mut fsum1 : Vec<f64> = Vec::with_capacity(dim);
        //  let mut fsum2 : Vec<f64> = Vec::with_capacity(dim);
        //  fsum1.extend(&fv1[(jtw1-1) * dim..jtw1 * dim]);
        //  fsum2.extend(&fv1[(jtw2-1) * dim..jtw2 * dim]);

        let mut fsum1 = fv1[(jtw1-1) * dim..jtw1 * dim].to_vec();
        let mut fsum2 = fv1[(jtw2-1) * dim..jtw2 * dim].to_vec();

        for k in 0..dim {
            //fsum1.push(fv1[(jtw1-1) * dim + k] + fv2[(jtw1-1) * dim + k]);
            //fsum2.push(fv1[(jtw2-1) * dim + k] + fv2[(jtw2-1) * dim + k]);
            fsum1[k] += fv2[(jtw1-1) * dim + k];
            fsum2[k] += fv2[(jtw2-1) * dim + k];
        }

        for k in 0..dim {
            resg[k] += wg[j - 1] * fsum2[k];
            resk[k] += wgk[jtw1 - 1] * fsum1[k] + wgk[jtw2 - 1] * fsum2[k] ;
            resabs[k] += wgk[jtw1 - 1] * (fv1[(jtw1-1) * dim + k].abs() + fv2[(jtw1-1) * dim + k].abs())
                + wgk[jtw2 - 1] * ( fv1[(jtw2-1) * dim + k].abs() + fv2[(jtw2-1) * dim + k].abs());
        }
    }

    if M/2 != (M+1)/2 {
        let jtw1 = M;
        let absc = hlgth * xgk[jtw1 - 1];
        fv1.append(&mut f(centr - absc));
        fv2.append(&mut f(centr + absc));

        let mut fsum = fv1[(jtw1-1) * dim..jtw1 * dim].to_vec();

        for k in 0..dim {
            fsum[k] +=  fv2[(jtw1-1) * dim + k];
        }
        for k in 0..dim {
            resk[k] += wgk[jtw1 - 1] * fsum[k] ;
            resabs[k] += wgk[jtw1 - 1] * (fv1[(jtw1-1) * dim + k].abs() + fv2[(jtw1-1) * dim + k].abs());
        }
    }

    //println!("Only for : {:?}",start.elapsed());

    let mut reskh =resk.clone();
    for k  in 0..dim {
        reskh[k] *= 0.5;
    }

    let mut resasc = vec![0.0; dim];
    for k in 0..dim {
        resasc[k] = wgk[M] * ( fc[k] - reskh[k] ).abs();
    }

    for j in 1..M + 1 {
        for k in 0..dim {
            resasc[k] += wgk[j - 1] * ((fv1[(j - 1)*dim + k] - reskh[k]).abs() + (fv2[(j - 1)*dim + k] -
                reskh[k]).abs());
        }
    }

    let mut  result = resk.clone();

    for k in 0..dim {
        result[k] *= hlgth;
        resabs[k] *= dhlgth;
        resasc[k] *= dhlgth;
    }

    let mut abserr = 0.0;
    let mut resabs_scalar = 0.0;
    let mut resasc_scalar = 0.0;

    for k in 0..dim {
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

    //println!("Return : {:?}",start.elapsed());

    (result, abserr, round_error)
}


#[cfg(test)]
mod tests {
    use std::time::Instant;

    #[test]
    fn arrayvsvec(){
        for i in 0..10 {
            let time = Instant::now();
            let mut a = [0.0; 1000];
            for k in 0..1000 {
                a[k] = k as f64;
            }
            println!("array : {:?}", time.elapsed());

            let time = Instant::now();
            let mut b = Vec::with_capacity(1000);
            for k in 0..1000 {
                b.push(k as f64);
            }
            println!("vec : {:?}", time.elapsed());

            let time = Instant::now();
            let mut c = [0.0; 1000].to_vec();
            for k in 0..1000 {
                c[k] = k as f64;
            }
            println!("vec 2 : {:?}", time.elapsed());

            print!("sum1 : {:?}", a.iter().sum::<f64>());
            print!("sum2: {:?}", b.iter().sum::<f64>());
            print!("sum3: {:?}", c.iter().sum::<f64>());
            println!("{i}");

        }

    }
}

