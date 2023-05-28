use std::cmp::Ordering;
use std::hash;
use std::sync::Arc;


#[derive(Clone)]
pub struct FnVecGen<const N:usize>{
    pub components : Arc< dyn Fn(f64)->[f64; N] + Send + Sync>
}

pub const EPMACH : f64 = f64::EPSILON;          // the largest relative spacing.
pub const UFLOW : f64 = f64::MIN_POSITIVE;      // the smallest positive magnitude.
//pub const OFLOW : f64 = f64::MAX;               // oflow is the largest positive magnitude.


pub fn norm_vec(v : &[f64]) -> f64{
    let mut norm = 0.0;
    for comp in v{
        norm += comp.powi(2);
    }
    norm = norm.sqrt();
    norm
}

pub fn res_update(v : &mut[f64], w: &[f64], z : &[f64], y : &[f64]){
    for k  in 0..v.len(){
        v[k] += w[k] + z[k] - y[k];
    }
}

pub fn add_res( v : &mut[f64], w : &[f64]){
    for k  in 0..v.len(){
        v[k] += w[k];
    }
}

#[derive(Debug,Clone)]
pub struct HeapItem {
    pub interval : (f64,f64),
    pub err : f64,
}

impl HeapItem {
    pub fn new( interval : (f64,f64) , err : f64) -> Self{
        Self{ interval,err}
    }
}

impl Eq for HeapItem{}

impl PartialEq for HeapItem{
    fn eq(&self, other: &Self) -> bool {
        self.err == other.err
    }
}

impl Ord for HeapItem {
    fn cmp(&self, other: &Self) -> Ordering {
        (self.err).partial_cmp(&other.err).unwrap()
    }
}

impl PartialOrd for HeapItem {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}



#[derive(Debug,Clone)]
pub struct Myf64{
    pub x : f64,
}
impl Myf64 {
    fn key(&self) -> u64 {
        self.x.to_bits()
    }
}

impl hash::Hash for Myf64 {
    fn hash<H>(&self, state: &mut H)
        where
            H: hash::Hasher,
    {
        self.key().hash(state)
    }
}

impl PartialEq for Myf64 {
    fn eq(&self, other: &Myf64) -> bool {
        self.key() == other.key()
    }
}

impl Eq for Myf64{}