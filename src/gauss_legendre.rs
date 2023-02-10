use crate::*;
extern crate nalgebra as na;
use na::*;

pub trait Gauss{
    fn integrate(&self,funct : &impl Function) -> f64 ;
}



pub struct GaussLegendre {
    pub nodes: Vec<f64>,
    pub weights: Vec<f64>,
    pub scale_factor : f64,
}

impl GaussLegendre{
    pub fn new(number_of_points : usize) -> Self{
        let mut matrix = DMatrix::from_element(number_of_points,number_of_points,0.0);
        for i in 0..number_of_points{
            for j in 0..number_of_points{
                if j == i+1 {
                    matrix[(i,j)] = (j as f64) / ( 4.0 * (j as f64).powi(2) - 1.0).sqrt();
                }
                else if j as i32 == i as i32 -1 {
                    matrix[(i,j)] = i as f64 / ( 4.0 * (i as f64).powi(2) - 1.0).sqrt();
                }

            }
        }
        let eigen= SymmetricEigen::new(matrix);
        let weights : Vec<f64> = (eigen.eigenvectors.row(0).map(|x| x.powi(2)) * 2.0)
            .data
            .as_vec()
            .clone();
        let nodes : Vec<f64>  = eigen.eigenvalues.data.as_vec().clone();
        Self{
           nodes,
           weights,
           scale_factor : 0.5,
       }
    }
    fn argument_transformation(x: f64, a: f64, b: f64) -> f64 {
        0.5 * ((b - a) * x + (b + a))
    }

    fn scale_factor(a: f64, b: f64) -> f64 {
        0.5 * (b - a)
    }
}


impl Gauss for GaussLegendre{
    fn integrate(&self,funct : &impl Function) -> f64 {
        let mut tot = 0.0;
        for i in 0..self.weights.len(){
            tot += funct.evaluate(&(GaussLegendre::argument_transformation(self.nodes[i],0.0,1.0) ) ) * self.weights[i];
        }
        tot * GaussLegendre::scale_factor(0.0,1.0)
    }
}