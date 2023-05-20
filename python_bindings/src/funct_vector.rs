use std::sync::Arc;


pub struct FnVecGen<const n :usize>{
    pub components : Arc< dyn Fn(f64)->[f64;n] + Send + Sync>
}

