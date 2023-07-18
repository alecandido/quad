use crate::constants::*;

pub fn qk_quadrature_ndarray<const M: usize, F>(
    f: F,
    a: f64,
    b: f64,
    xgk: &[f64; M],
    wgk: &[f64],
    wg: &[f64],
) -> (ndarray::Array1::<f64>, f64, f64)
where
    F: Fn(f64) -> Vec<f64>,
{
    let hlgth: f64 = 0.5 * (b - a);
    let dhlgth: f64 = hlgth.abs();
    let centr: f64 = 0.5 * (b + a);
    let fc = f(centr);
    let dim = fc.len();
    let mut fv1 = Vec::with_capacity(dim * M);
    let mut fv2 = Vec::with_capacity(dim * M);
    let mut resg = {
        if M % 2 == 1 {
            ndarray::Array1::<f64>::from(fc.clone()) * wg[(M + 1) / 2 - 1]
        } else {
            ndarray::Array1::<f64>::zeros(dim)
        }
    };
    let mut resk = ndarray::Array1::<f64>::from(fc.clone()) * wgk[M];
    let mut resabs = Vec::with_capacity(dim);

    for k in 0..dim {
        resabs.push(resk[k].abs());
    }

    let mut resabs = ndarray::Array1::<f64>::from(resabs);

    for j in 1..M / 2 + 1 {
        let jtw1 = 2 * j - 1;
        let jtw2 = 2 * j;

        let absc1 = hlgth * xgk[jtw1 - 1];
        let absc2 = hlgth * xgk[jtw2 - 1];

        fv1.append(&mut f(centr - absc1));
        fv1.append(&mut f(centr - absc2));
        fv2.append(&mut f(centr + absc1));
        fv2.append(&mut f(centr + absc2));

        let mut fsum1 = fv1[(jtw1 - 1) * dim..jtw1 * dim].to_vec();
        let mut fsum2 = fv1[(jtw2 - 1) * dim..jtw2 * dim].to_vec();

        for k in 0..dim {
            fsum1[k] += fv2[(jtw1 - 1) * dim + k];
            fsum2[k] += fv2[(jtw2 - 1) * dim + k];
        }

        let fsum1 = ndarray::Array1::<f64>::from(fsum1);
        let fsum2 = ndarray::Array1::<f64>::from(fsum2);

        resg += &(fsum2.clone() * wg[j - 1]);
        resk += &(fsum1 * wgk[jtw1 - 1]);
        resk += &(fsum2 * wgk[jtw2 - 1]);

        for k in 0..dim {
            resabs[k] += wgk[jtw1 - 1]
                * (fv1[(jtw1 - 1) * dim + k].abs() + fv2[(jtw1 - 1) * dim + k].abs())
                + wgk[jtw2 - 1]
                    * (fv1[(jtw2 - 1) * dim + k].abs() + fv2[(jtw2 - 1) * dim + k].abs());
        }
    }

    if M / 2 != (M + 1) / 2 {
        let jtw1 = M;
        let absc = hlgth * xgk[jtw1 - 1];
        fv1.append(&mut f(centr - absc));
        fv2.append(&mut f(centr + absc));

        let mut fsum = fv1[(jtw1 - 1) * dim..jtw1 * dim].to_vec();

        for k in 0..dim {
            fsum[k] += fv2[(jtw1 - 1) * dim + k];
        }

        let mut fsum = ndarray::Array1::<f64>::from(fsum);
        resk += &(fsum * wgk[jtw1 - 1]);

        for k in 0..dim {
            resabs[k] +=
                wgk[jtw1 - 1] * (fv1[(jtw1 - 1) * dim + k].abs() + fv2[(jtw1 - 1) * dim + k].abs());
        }
    }

    let mut reskh = resk.clone();

    for k in 0..dim {
        reskh[k] *= 0.5;
    }

    let mut resasc = vec![0.0; dim];
    for k in 0..dim {
        resasc[k] = wgk[M] * (fc[k] - reskh[k]).abs();
    }

    for j in 1..M + 1 {
        for k in 0..dim {
            resasc[k] += wgk[j - 1]
                * ((fv1[(j - 1) * dim + k] - reskh[k]).abs()
                    + (fv2[(j - 1) * dim + k] - reskh[k]).abs());
        }
    }

    let mut result = resk.clone();
    result *= hlgth;

    for k in 0..dim {
        resabs[k] *= dhlgth;
        resasc[k] *= dhlgth;
    }

    let mut abserr = 0.0;
    let mut resabs_scalar = 0.0;
    let mut resasc_scalar = 0.0;

    for k in 0..dim {
        abserr += (((resk[k] - resg[k]) * hlgth).abs()).powi(2);
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
