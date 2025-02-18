#![allow(unused)]
#![allow(mixed_script_confusables)]
#![allow(non_snake_case)]

use std::{fs, io, thread};
use std::time::{Instant, Duration};
use std::str::FromStr;

fn main() -> io::Result<()> {
    println!("Runge-Kutta 4-Stage Numerical Integration");
    
    let (mut t_0, mut t_f, mut Δt, mut y_0): (f32, f32, f32, f32) = (0.0, 0.0, 0.0, 0.0); //  dep+ind var params
    let (mut s, mut n): (u8, u32) = (0, 0); //  stage number, step number
    let integrand = |λ_t: f32| -> f32 { //  closure that takes up to two arguments:: λ_t: f32, λ_y: f32
        rk_meth::integrating_functions::prod_func(λ_t)
    };
    let solution = |λ_t: &f32| -> f32 { //  (yet another) closure that takes up to two arguments:: λ_t: f32, λ_y: f32
        rk_meth::exact_solutions::prod_sol(λ_t)
    };
    let errata = |data: &mut Vec<Vec<f32>>| {
        for row in data.iter_mut() {
            let (t, y) = (row[0], row[1]);
            let (abs, rel) = mathematics::stats::error(&y, &solution(&t));
            
            row.push(solution(&t));
            row.push(abs);
            row.push(rel);
        }
    };
    
    //  code block for threading and user input(?)
    
    /*     
        thread::Builder::new().name("timekeeper".to_string());
        thread::Builder::new().name("error-stats".to_string());
    */

    /*
        {
            let mut param_check = |i: &u8| {
                let (aff, neg): (char, char) = ('y', 'n');
                let mut usr_inp = || -> String {
                    let mut input = String::with_capacity(8);
                    io::stdin()
                        .read_line(&mut input)
                        .expect("Failed to read line!");
                    
                    return input
                };
        
                match i {
                    0 => {  // stage number prompt, to be done later
                        todo!();
                    },
                    1 => {  // dep-var initial conditions & increment parameter
                        println!("Please specify the initial time(t_0):");
                        let dicks = u32::from_str(&usr_inp()).expect("t_0 not set!");
        
                        println!("Please specify the final time(t_f):");
                        t_f = usr_inp()
                            .parse::<f32>()
                            .expect("t_f not set!");
                        
                        println!("Would you like to integrate by fixing the\n(1) step number\nor\n(2) step size?");
                        let step_param: u8 = usr_inp()
                            .parse::<u8>()
                            .expect("Iteration parameter not correctly set!");
                        let fixed_step = |a: &f32, b: &f32, k: &u32| -> f32 {
                            return (b - a).abs() / (*k as f32)
                        };
        
                        match step_param {
                            1 => {
                                println!("Please specify the number of steps/iterations (n):");
                                n = usr_inp().parse::<u32>().expect("t_f not set!");
                                Δt = fixed_step(&t_0, &t_f, &n);
                            },
                            2 => {
                                println!("Please specify the step size (Δt):");
                                Δt = usr_inp().parse::<f32>().expect("t_f not set!");
                                Δt.abs();
                            },
                            _ => {
                                println!("bruh, wtf. this is being done with n = 1000, then :T");
                                n = 1000;
                                Δt = fixed_step(&t_0, &t_f, &n);
                            },
                        }
                    },
                    2 => {  // 
                        todo!();
                    },
                    3 => {  // ind-var boundary conditions
                        todo!();
                    },
                    _ => panic!("This literally isn't supposed to happen wtf"),
                }
            };

            param_check(&1);
        }
    */

    let dp_arr: (f32, f32, f32, f32) = (1.0, 10.0, 1.0e-2, 1.0); //  dummy-parameter array, to be axed
    (t_0, t_f, Δt, y_0) = (dp_arr.0, dp_arr.1, dp_arr.2, dp_arr.3); //  to be axed

    assert!(Δt.is_sign_positive()); //  so that nothing fucking breaks
    
    
    println!("Performing integration with the following parameters...
        t: [{}, {}]; Δt: {}, n: {} steps...",
        &t_0, &t_f, &Δt, ((&t_0 - &t_f) / &Δt).abs().ceil()
    );
    let now = Instant::now();
    let (result, mut data) = rk_meth::rk_4_std_aut(t_0, t_f, Δt, integrand); //  integration~
    let elapsed = now.elapsed();
    
    let y_int: f32 = solution(&t_f) - solution(&t_0);
    let (err, rel_err) = mathematics::stats::error(&result, &y_int);
    println!("sum: {}, ε: {} ({})‱\ntime elapsed: {} μs",
        &result, &err, &rel_err * 1.0e3, &elapsed.as_secs_f32() * 1.0e6
    );

    errata(&mut data);
    
    fs::write("errchk.csv", piping::export(&data))?;

    Ok(())
}

pub mod rk_meth {
    pub fn rk_4_std_aut<T: Fn(f32) -> f32>(t_a: f32, t_b: f32, Δt: f32, Δy_t: T) -> (f32, Vec<Vec<f32>>) {
        // initialsation of output, time, and function
        let (mut sum, mut t_i): (f32, f32) = (0.0, t_a);
        let mut integral: Vec<Vec<f32>> = Vec::with_capacity(
            capacity(t_a, t_b, Δt)
        ); // stores the total integration over iterations

        integral.push(vec![t_a, sum]);

        while t_i < t_b {
            let mut row: Vec<f32> = Vec::with_capacity(2);
            
            // slope generators
            let k_1 = Δy_t(t_i); // * Δt / 6.0
            let k_mid = Δy_t(t_i + Δt / 2.0) * 2.0;
            let k_4 = Δy_t(t_i + Δt);

            let mut k_vec: Vec<f32> = vec![k_1, k_mid, k_mid, k_4];
            k_vec.iter_mut().for_each(|mut u: &mut f32| *u = Δt * *u / 6.0);
            
            // update values
            sum += k_vec.iter().sum::<f32>();
            t_i += Δt;

            row.push(t_i);
            row.push(sum);

            integral.push(row);
        }

        return (sum, integral)
    }

    pub fn rk_4_std<T: Fn(f32, f32) -> f32>((t_0, y_0): (f32, f32), t_f: f32, Δt: f32, f_ty: T) -> (f32, Vec<Vec<f32>>) {
        let (mut t_i, mut y_i): (f32, f32) = (t_0, y_0);

        let mut integral: Vec<Vec<f32>> = Vec::with_capacity(
            capacity(t_0, t_f, Δt)
        ); // stores the total integration over iterations

        integral.push(vec![t_0, y_0]);

        while t_i < t_f {
            let mut row: Vec<f32> = Vec::with_capacity(2); // row initialisation
            
            let k_1 = f_ty(t_i, y_i);
            let k_2 = f_ty(Δt.mul_add(0.5, t_i), 0.5 * Δt.mul_add(k_1, y_i));
            let k_3 = f_ty(Δt.mul_add(0.5, t_i), 0.5 * Δt.mul_add(k_2, y_i));
            let k_4 = f_ty(t_i + Δt, Δt.mul_add(k_3, y_i));
            let mut k_vec = vec![k_1, 2.0 * k_2, 2.0 * k_3, k_4];

            k_vec.iter_mut().for_each(|u| *u = Δt * *u / 6.0);
            
            // update values
            t_i += Δt;
            y_i += k_vec.iter().sum::<f32>();

            row.push(t_i);
            row.push(y_i);
            
            integral.push(row);
        }

        return (y_i - y_0, integral)
    }

    fn capacity(a: f32, b: f32, Δt: f32) -> usize {
        assert!(2.0 * Δt < b - a);
        ((b - a).abs() / Δt).ceil() as usize
    }

    pub mod lin_alg {
        use std::ops::{Add, Sub, Mul};
        use std::ops::Neg;

        #[derive(Clone)]
        pub struct Vector {
            form: VecForm,
            dim: usize,
            elements: Vec<f32>,
        }

        #[derive(Clone)]
        enum VecForm {
            Row,
            Col,
        }

        pub fn inner_prod(a: &Vector, b: &Vector) -> f32 {
            assert_eq!(a.dim, b.dim);

            let mut u: Vec<f32> = a.elements.clone();
            
            u.iter_mut().zip(b.elements.clone()).for_each(|(v, w)| *v = *v * w);
            return u.iter().sum()
        }

        impl Vector {
            pub fn dim(&self) -> usize {
                return self.dim
            }

            pub fn transpose(&mut self) {
                match &self.form {
                    VecForm::Row => {
                        self.form = VecForm::Col;
                    },
                    VecForm::Col => {
                        self.form = VecForm::Row;
                    },
                }
            }

            pub fn mag(&self) -> f32 {
                let mut mag: Vec<f32> = self.elements.clone();
                
                mag.iter_mut().for_each(|mut u| *u = u.powi(2));

                return mag.iter().sum::<f32>().sqrt()
            }

            pub fn neg(&mut self) {
                self.elements.iter_mut().for_each(|u| {
                    u.neg();
                });
            }

            pub fn angle(&self, other: &Self) -> f32 {
                return (inner_prod(self, other) / (self.mag() * other.mag())).acos()
            }
        }

        impl Add for Vector {
            type Output = Self;

            fn add(mut self, other: Self) -> Self {
                assert_eq!(self.dim, other.dim);
                
                self.elements.iter_mut().zip(other.elements).for_each(|(u, v)| *u = *u + v);

                return self
            }
        }

        impl Sub for Vector {
            type Output = Self;

            fn sub(self, mut other: Self) -> Self {
                other.neg();
                
                return self + other
            }
        }

        // impl Mul for Vector {}
        
        #[derive(Clone)]
        pub struct Matrix {
            dim: (usize, usize),
            elements: Vec<f32>,
        }

        #[derive(Clone)]
        pub struct FlatMatrix {
            symmetric: bool,
            dim: usize,
            elements: Vec<f32>,
        }
    }
    
    pub mod integrating_functions {
        use std::ops::Neg;
        
        // autonomous functions
        pub fn quad_func (t: f32) -> f32 {
            return t.powi(2)
        }
    
        pub fn cube_func (t: f32) -> f32 {
            return t.powi(3)
        }
    
        pub fn prod_func (t: f32) -> f32 {
            return t.neg().exp() * t.powi(2)
        }
    
        // non-autonomous functions
        pub fn simple(t: f32, y: f32) -> f32 {
            return y * t.powi(-2)
        }

        fn well_behaved(t: f32, y: f32) -> f32 {
            return y.neg() + t.cos()
        }
    
        fn sing_discon(t: f32, y: f32) -> f32 {
            assert_ne!(t, 0.0);
            
            return 2.0 * y / t + t.powi(2) * t.neg().exp()
        }
    
        fn dual_discon(t: f32, y: f32) -> f32 {
            assert_ne!(t, -1.0);
            assert_ne!(t, 3.0);
            
            let q_t = t.powi(2) - 2.0 * t - 3.0;
            return y * q_t.recip()
        }
    }

    pub mod exact_solutions {
        use std::f32::consts::E;
        use std::ops::Neg;
        
        //  'simple/explicit/seperable(?)' integrals
        pub fn quad_sol(t: &f32) -> f32 {
            return t.powi(3) / 3.0
        }

        pub fn prod_sol(t: &f32) -> f32 {
            return t.neg().exp().neg() * (t.powi(2) + 2.0 * t + 2.0) // dyscalculia will be the death of me. lmao
        }

        pub fn cube_sol(t: &f32) -> f32 {
            return t.powi(4) * 0.25
        }

        //  linear ODEs
        pub fn simple_sol(t: &f32) -> f32 {
            return E.powi(2) * (2.0 * t.powi(-3).neg()).exp()
        }
    }
}

pub mod piping {
    // use std::ops::Index;

    pub mod version {
        const INIT_BUILD: char = 'α';
        
        struct Version {
            release: bool,
            build: Option<char>,
            version: u8,
            revision: u8,
            branch: Option<char>,
        }

        impl Version {
            fn tag(release: bool) -> Version {
                let dummy_version = Version {
                    release: false,
                    build: Some(INIT_BUILD),
                    version: 4,
                    revision: 20,
                    branch: Some('ω'),
                };

                return dummy_version
            }
        }
    }

    pub mod interface {
        use std::io;
        use std::fs;

        use std::fs::File;
        use std::io::Read;

        pub const AFF: char = 'y';
        pub const NEG: char = 'n';
        
        pub fn script(prompt_tree: File) -> io::Result<()> {
            unimplemented!();
            
            // let mut script: String = fs::read_to_string(dummy_file)?;
            let branch_prompt: &str = "!!!";




            return Ok(())
        }
    }

    const NEWLINE: &str = "\n";
    const SPCR: &str = ", ";

    pub fn export(data: &Vec<Vec<f32>>) -> String {
        let table_header: &str = "'n', 't(n)', 'y(t)', 'y_xct', 'ε_abs', 'ε_rel'\n";
        let mut output: String = String::new();
        let table = to_table(data);

        output += table_header;
    
        for row in table.iter() {
            row.iter().for_each(|chr| output.push(*chr));
            output += NEWLINE;
        }
    
        return output
    }

    fn to_table(data: &Vec<Vec<f32>>) -> Vec<Vec<char>> {
        let mut lines: Vec<Vec<char>> = vec![];
        let mut index = 0;

        for row in data.iter() {
            let c_0: String = index.to_string();
            let (c_1, c_2, c_3, c_4, c_5) = (
                row[0].to_string(),
                row[1].to_string(),
                row[2].to_string(),
                row[3].to_string(),
                row[4].to_string()
            );
            let output: String = c_0+&SPCR+&c_1+&SPCR+&c_2+&SPCR+&c_3+SPCR+&c_4+SPCR+&c_5;
            lines.push(output.chars().collect::<Vec<char>>());
            index += 1;
        }

        lines
    }

    fn meta_header<T: Copy> (parameters: &Vec<T>) -> String {
        let mut metadata = String::new();
        let com_tag: &str = "## ";
        let title: &str = "Runge-Kutta Numerical Integration";

        // let (app_info, parameters, stats): (_, _, _);
    
        return metadata
    }
}

pub mod mathematics {
    pub mod stats {
        pub fn mean(set: &Vec<f32>) -> f32 {
            return (set.iter().sum::<f32>() / set.len() as f32) as f32
        }

        pub fn variance(set: &Vec<f32>) -> f32 {
            let mut sum: f32 = 0.0;
            let mean = mean(set);

            set.iter().for_each(|x| sum += (x - mean).powi(2));

            return sum / (set.len() - 1) as f32
        }

        pub fn std_dev(set: &Vec<f32>) -> f32 {
            return variance(set).sqrt()
        }

        pub fn error(x: &f32, xct: &f32) -> (f32, f32) {
            let abs: f32 = (x - xct).abs();
            
            return (abs, abs / xct)
        }
    }

    mod linear_algebra {
        use std::iter::zip;
    
        fn levi_civita(x: &Vec<u8>) -> Option<bool> {
            let mut y = x.clone();
            let mut e = Vec::with_capacity(x.len());
    
            for i in [0,x.len()] {
                e.push(i);
            }
    
            // assert_eq!(x.len(), e.len()); // length check
            let e_check = zip(x.iter(), y.iter());
            let e_filter = e_check.filter(|(u, v)| u == v).collect::<Vec<_>>();
    
            if e_filter.len() == 0 {
                return None
            }
    
            let mut y_prime = y.clone();
    
            y_prime.rotate_right(1);
            
            let perm_pair = zip(y.iter(), y_prime.iter());
            let mut is_signed: bool = false;
            
            perm_pair.for_each(
                |(u, v)| if u < v { is_signed = !is_signed }
            ); 
    
            Some(is_signed)
        }

        fn tri_num(n: u8) -> u8 {
            return (n * (n + 1)) /2
        }
        
        enum Sign {
            NonZero(bool),
            Zero,
        }
        
        enum RealNum {
            Nat
                (Sign, u32),
            Int(
                (Sign, u32)
            ),
            Rat(
                (Sign, u16, u16)
            ),
            IrrRat, // lmao, fuck this one in particular
            Real(f32),
    
        }
    
        enum PlaneCoord {
            Carts,
            Polar,
        }
    
        struct CompNum {
            z: RealNum,
            form: PlaneCoord,
        }
    
        enum Element {
            BitMap(bool),
            Index(usize),
            Reals,
            Comp,
            Func,
        }
    
        enum Object {
    
        }
        
        struct Tensor {}
    }
}

//  testing, benchmarks, etc
#[cfg(test)]
mod tests{
    use super::*;
    use std::time::{Duration, Instant};

    // linear first-order ode homogeneous
    fn well_behaved_exact(c: f32, t: f32) -> f32 {
        
    }

    fn sing_discon_exact(c: f32, t: f32) -> f32 {
        assert_ne!(t, 0.0);
        
        unimplemented!();
    }

    fn dual_discon_exact(c: f32, t: f32) -> f32 {
        assert_ne!(t, -1.0);
        assert_ne!(t, 3.0);
        
        return (-0.25 * (t + 1.0).abs().ln() + 0.25 * (t - 3.0).abs().ln()).exp() * c
    }
}