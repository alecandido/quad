use criterion::BenchmarkId;
use criterion::Criterion;
use criterion::{criterion_group, criterion_main};
use quad::*;
use rgsl::*;
use std::{thread, time};
use std::sync::Arc;
use quad::constants::FnVec;

fn bench_qag(c: &mut Criterion) {
    let mut group = c.benchmark_group("Qag");
    let rgsl_input = (0.0, 50000.0, 1.0e-2, 0.0, 1000000, GaussKronrodRule::Gauss61);
    group.bench_with_input(BenchmarkId::new("rgsl_qag", rgsl_input.0), &rgsl_input, |b, &inp| {
        let f = |x:f64| x.cos();
        b.iter(|| {
            let work = IntegrationWorkspace::new(1000000);
            work.expect("a").qag(f,inp.0,inp.1,inp.2,inp.3,inp.4,inp.5)}
        );

    });

    let input = (0.0, 50000.0, 1.0e-2, 0.0, 6, 1000000, false);
    group.bench_with_input(BenchmarkId::new("My_qag", "My_qag"), &input, |b, &inp| {
        let f = |x:f64| vec![x.cos()];
        b.iter(|| integrate(&f,inp.0,inp.1,inp.2,inp.3,inp.4,inp.5,[0.0;0].to_vec(),inp.6));
    });

    group.bench_with_input(BenchmarkId::new("My_qag_par", "My_qag_par"), &input, |b, &inp| {
        let f = FnVec{ components : Arc::new( |x:f64| vec![x.cos()])};
        b.iter(|| integrate_par(&f,inp.0,inp.1,inp.2,inp.3,inp.4,inp.5,[0.0;0].to_vec(),inp.6));
    });


    group.finish();
}


fn bench_qag_delay(c: &mut Criterion) {
    let mut group = c.benchmark_group("Qag");
    let rgsl_input = (0.0, 500.0, 1.0e-2, 0.0, 1000000, GaussKronrodRule::Gauss61);
    group.bench_with_input(BenchmarkId::new("Rgsl_qag", "Rgsl_qag"), &rgsl_input, |b, &inp| {
        let f = |x:f64| {
            thread::sleep(time::Duration::from_nanos(100));
            x.cos()
        };
        b.iter(|| {
            let work = IntegrationWorkspace::new(1000000);
            work.expect("a").qag(f,inp.0,inp.1,inp.2,inp.3,inp.4,inp.5)}
        );

    });

    let input = (0.0, 500.0, 1.0e-2, 0.0, 6, 1000000, false);
    group.bench_with_input(BenchmarkId::new("My_qag", "My_qag"), &input, |b, &inp| {
        let f = |x:f64| {
            thread::sleep(time::Duration::from_nanos(100));
            vec![x.cos()]
        };
        b.iter(|| integrate(&f,inp.0,inp.1,inp.2,inp.3,inp.4,inp.5,[0.0;0].to_vec(),inp.6));
    });

    group.bench_with_input(BenchmarkId::new("My_qag_par", "My_qag_par"), &input, |b, &inp| {
        let f = FnVec{ components : Arc::new( |x:f64| {
            thread::sleep(time::Duration::from_nanos(100));
            vec![x.cos()]
        })};
        b.iter(|| integrate_par(&f,inp.0,inp.1,inp.2,inp.3,inp.4,inp.5,[0.0;0].to_vec(),inp.6));
    });

    group.finish();
}

criterion_group!(benches, bench_qag_delay);
criterion_main!(benches);

