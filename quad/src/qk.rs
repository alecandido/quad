use crate::constants::*;
use ndarray::Axis;

pub fn qk_quadrature<const M: usize, F>(
    f: F,
    a: f64,
    b: f64,
    xgk: &[f64; M],
    wgk: &[f64],
    wg: &[f64],
) -> (ndarray::Array1<f64>, f64, f64)
where
    F: Fn(f64) -> ndarray::Array1<f64>,
{
    let hlgth: f64 = 0.5 * (b - a);
    let dhlgth: f64 = hlgth.abs();
    let centr: f64 = 0.5 * (b + a);
    let fc = f(centr);
    let dim = fc.len();
    let mut fv1 = ndarray::Array1::<f64>::zeros(0);
    let mut fv2 = ndarray::Array1::<f64>::zeros(0);
    let mut resg = {
        if M % 2 == 1 {
            &fc * wg[(M + 1) / 2 - 1]
        } else {
            ndarray::Array1::<f64>::zeros(dim)
        }
    };
    let mut resk = &fc * wgk[M];
    let mut resabs = resk.map(|x| x.abs());

    for j in 1..M / 2 + 1 {
        let jtw1 = 2 * j - 1;
        let jtw2 = 2 * j;

        let absc1 = hlgth * xgk[jtw1 - 1];
        let absc2 = hlgth * xgk[jtw2 - 1];

        let f11 = f(centr - absc1);
        let f12 = f(centr - absc2);
        let f21 = f(centr + absc1);
        let f22 = f(centr + absc2);

        fv1.append(Axis(0), f11.view());
        fv1.append(Axis(0), f12.view());
        fv2.append(Axis(0), f21.view());
        fv2.append(Axis(0), f22.view());

        //resabs += &(&(f11.map(|x| x.abs()) + &(f21.map(|x| x.abs()) ) ) * wgk[jtw1 -1]);
        //resabs += &(&(f12.map(|x| x.abs()) + &(f22.map(|x| x.abs()) ) ) * wgk[jtw1 -1]);

        for k in 0..dim {
            resabs[k] += wgk[jtw1 - 1] * (f11[k].abs() + f21[k].abs())
                + wgk[jtw2 - 1] * (f12[k].abs() + f22[k].abs());
        }

        let fsum1 = f11 + f21;
        let fsum2 = f12 + f22;

        resg += &(&fsum2 * wg[j - 1]);
        resk += &(fsum1 * wgk[jtw1 - 1]);
        resk += &(fsum2 * wgk[jtw2 - 1]);
    }

    if M / 2 != (M + 1) / 2 {
        let jtw1 = M;
        let absc = hlgth * xgk[jtw1 - 1];
        let f1 = f(centr - absc);
        let f2 = f(centr + absc);
        fv1.append(Axis(0), f1.view());
        fv2.append(Axis(0), f2.view());

        for k in 0..dim {
            resabs[k] += wgk[jtw1 - 1] * (f1[k].abs() + f2[k].abs());
        }

        resk += &((&f1 + &f2) * wgk[jtw1 - 1]);
    }

    let reskh = &resk * 0.5;

    let mut resasc = (&fc - &reskh).map(|x| x.abs() * wgk[M]);

    for j in 1..M + 1 {
        for k in 0..dim {
            resasc[k] += wgk[j - 1]
                * ((fv1[(j - 1) * dim + k] - reskh[k]).abs()
                    + (fv2[(j - 1) * dim + k] - reskh[k]).abs());
        }
    }

    let result = &resk * hlgth;

    resabs *= dhlgth;
    resasc *= dhlgth;

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
