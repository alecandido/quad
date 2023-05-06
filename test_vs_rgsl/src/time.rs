use std::time::Instant;

pub fn time <X: std::fmt::Debug>(x : X) -> (){
    let start = Instant::now();
    let a = x;
    println!("{:?}",a);
    println!("{:?}",start.elapsed());
}