use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::{criterion_group, criterion_main};
use quad::constants::FnVec;
use quad::qag::*;
use rgsl::*;
use std::sync::Arc;
use std::{thread, time};

fn qag_delay(c: &mut Criterion) {
    let mut group = c.benchmark_group("Qag");
    let rgsl_input = (0.0, 500.0, 1.0e-2, 0.0, 1000000, GaussKronrodRule::Gauss61);
    let range = [0, 1, 2, 3, 4, 5, 6];
    for z in range {
        group.bench_with_input(BenchmarkId::new("Rgsl_qag", z), &rgsl_input, |b, &inp| {
            let f = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10_i32.pow(z) as u64));
                x.cos()
            };
            b.iter(|| {
                let work = IntegrationWorkspace::new(1000000);
                work.expect("a")
                    .qag(f, inp.0, inp.1, inp.2, inp.3, inp.4, inp.5)
            });
        });

        let input = (0.0, 500.0, 1.0e-2, 0.0, 6, 1000000, false);
        group.bench_with_input(BenchmarkId::new("My_qag", z), &input, |b, &inp| {
            let f = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10_i32.pow(z) as u64));
                vec![x.cos()]
            };
            b.iter(|| {
                integrate(
                    &f,
                    inp.0,
                    inp.1,
                    inp.2,
                    inp.3,
                    inp.4,
                    inp.5,
                    [0.0; 0].to_vec(),
                    inp.6,
                )
            });
        });

        group.bench_with_input(BenchmarkId::new("My_qag_par", z), &input, |b, &inp| {
            let f = FnVec {
                components: Arc::new(|x: f64| {
                    thread::sleep(time::Duration::from_nanos(10_i32.pow(z) as u64));
                    vec![x.cos()]
                }),
            };
            b.iter(|| {
                integrate_par(
                    &f,
                    inp.0,
                    inp.1,
                    inp.2,
                    inp.3,
                    inp.4,
                    inp.5,
                    [0.0; 0].to_vec(),
                    inp.6,
                )
            });
        });
    }
    group.finish();
}

fn fn_lenght(c: &mut Criterion) {
    let mut group = c.benchmark_group("Fn_lenght");
    let rgsl_input = (0.0, 500.0, 1.0e-2, 0.0, 1000000, GaussKronrodRule::Gauss61);
    let range = [1, 2, 3, 4, 5];
    for z in range {
        group.bench_with_input(BenchmarkId::new("Rgsl_qag", z), &rgsl_input, |b, &inp| {
            let f = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10000));
                x.cos()
            };
            b.iter(|| {
                for _k in 0..z {
                    let work = IntegrationWorkspace::new(1000000);
                    let _res = work
                        .expect("a")
                        .qag(f, inp.0, inp.1, inp.2, inp.3, inp.4, inp.5);
                }
            });
        });

        let input = (0.0, 500.0, 1.0e-2, 0.0, 6, 1000000, false);
        group.bench_with_input(BenchmarkId::new("My_qag", z), &input, |b, &inp| {
            let g = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10000));
                x.cos()
            };
            let f = |x: f64| {
                let mut v = vec![];
                for i in 0..z {
                    v.push(g(x))
                }
                v
            };

            b.iter(|| {
                integrate(
                    &f,
                    inp.0,
                    inp.1,
                    inp.2,
                    inp.3,
                    inp.4,
                    inp.5,
                    [0.0; 0].to_vec(),
                    inp.6,
                )
            });
        });

        group.bench_with_input(BenchmarkId::new("My_qag_par", z), &input, |b, &inp| {
            let g = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10000));
                x.cos()
            };
            let f = FnVec {
                components: Arc::new(|x: f64| {
                    let mut v = vec![];
                    for i in 0..z {
                        v.push(g(x))
                    }
                    v
                }),
            };
            b.iter(|| {
                integrate_par(
                    &f,
                    inp.0,
                    inp.1,
                    inp.2,
                    inp.3,
                    inp.4,
                    inp.5,
                    [0.0; 0].to_vec(),
                    inp.6,
                )
            });
        });
    }
    group.finish();
}

fn number_of_interval_subdivision(c: &mut Criterion) {
    let mut group = c.benchmark_group("Fn_lenght");
    let rgsl_input = (0.0, 500.0, 1.0e-2, 0.0, 1000000, GaussKronrodRule::Gauss61);
    let range = [1.0, 100.0, 1000.0, 10000.0];
    for z in range {
        group.bench_with_input(BenchmarkId::new("Rgsl_qag", z), &rgsl_input, |b, &inp| {
            let f = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10));
                x.cos()
            };
            b.iter(|| {
                let work = IntegrationWorkspace::new(1000000);
                let _res = work
                    .expect("a")
                    .qag(f, inp.0, z, inp.2, inp.3, inp.4, inp.5);
            });
        });

        let input = (0.0, 500.0, 1.0e-2, 0.0, 6, 1000000, false);
        group.bench_with_input(BenchmarkId::new("My_qag", z), &input, |b, &inp| {
            let f = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10));
                vec![x.cos()]
            };

            b.iter(|| {
                integrate(
                    &f,
                    inp.0,
                    z,
                    inp.2,
                    inp.3,
                    inp.4,
                    inp.5,
                    [0.0; 0].to_vec(),
                    inp.6,
                )
            });
        });

        group.bench_with_input(BenchmarkId::new("My_qag_par", z), &input, |b, &inp| {
            let f = FnVec {
                components: Arc::new(|x: f64| {
                    thread::sleep(time::Duration::from_nanos(10));
                    vec![x.cos()]
                }),
            };
            b.iter(|| {
                integrate_par(
                    &f,
                    inp.0,
                    z,
                    inp.2,
                    inp.3,
                    inp.4,
                    inp.5,
                    [0.0; 0].to_vec(),
                    inp.6,
                )
            });
        });
    }
    group.finish();
}

fn key(c: &mut Criterion) {
    let mut group = c.benchmark_group("Qag");
    let rgsl_input = (0.0, 2500.0, 1.0e-2, 0.0, 1000000);
    let range = [1, 2, 3, 4, 5, 6];
    for z in range {
        group.bench_with_input(BenchmarkId::new("Rgsl_qag", z), &rgsl_input, |b, &inp| {
            let f = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10));
                x.cos()
            };
            b.iter(|| {
                let work = IntegrationWorkspace::new(1000000);
                let key = match z {
                    1 => GaussKronrodRule::Gauss15,
                    2 => GaussKronrodRule::Gauss21,
                    3 => GaussKronrodRule::Gauss31,
                    4 => GaussKronrodRule::Gauss41,
                    5 => GaussKronrodRule::Gauss51,
                    6 => GaussKronrodRule::Gauss61,
                    _ => GaussKronrodRule::Gauss15,
                };
                work.expect("a")
                    .qag(f, inp.0, inp.1, inp.2, inp.3, inp.4, key)
            });
        });

        let input = (0.0, 2500.0, 1.0e-2, 0.0, 6, 1000000, false);
        group.bench_with_input(BenchmarkId::new("My_qag", z), &input, |b, &inp| {
            let f = |x: f64| {
                thread::sleep(time::Duration::from_nanos(10));
                vec![x.cos()]
            };
            b.iter(|| {
                integrate(
                    &f,
                    inp.0,
                    inp.1,
                    inp.2,
                    inp.3,
                    z,
                    inp.5,
                    [0.0; 0].to_vec(),
                    inp.6,
                )
            });
        });

        group.bench_with_input(BenchmarkId::new("My_qag_par", z), &input, |b, &inp| {
            let f = FnVec {
                components: Arc::new(|x: f64| {
                    thread::sleep(time::Duration::from_nanos(10));
                    vec![x.cos()]
                }),
            };
            b.iter(|| {
                integrate_par(
                    &f,
                    inp.0,
                    inp.1,
                    inp.2,
                    inp.3,
                    z,
                    inp.5,
                    [0.0; 0].to_vec(),
                    inp.6,
                )
            });
        });
    }
    group.finish();
}

criterion_group!(benches1, qag_delay);
criterion_group!(benches2, fn_lenght);
criterion_group!(benches3, number_of_interval_subdivision);
criterion_group!(benches4, key);
criterion_main!(benches1, benches2, benches3, benches4);
