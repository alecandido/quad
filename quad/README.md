# Usage

The primary function is the method 'integrate' of the struct Qag.
In order to use it the [Qag](quad/src/qag.rs) struct need to be initialized.

E.g:

```
let qag = Qag {
    key: 1,
    limit: 50,
    points: vec![0.0; 0],
    number_of_thread: 8,
    more_info: false,
};
```

An FnVec function must then be provided, which is constructed as follows:

E.g:

```
let f = FnVec {
                components: Arc::new(|x: f64| array![x.sin(), x.cos()]),
            };
```

The qag method `integrate`, which returns a `Result<QagIntegrationResult, QagError>`,
can therefore be called as follows:

```
let (a, b, epsabs, epsrel) = (0.0, 1.0, 1.0e-2, 1.0e-2);
let res : Result<QagIntegrationResult, QagError> = qag.integrate(&f, a, b, epsabs, epsrel);
```

Where `a` and `b` are the interval extremes, and `epsabs` and `epsrel` are the required absolute
and relative errors.
