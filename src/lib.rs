//! # rust-libcint
//!
//! The `rust-libcint` crate provides wrappers for libcint (C).
//!
//! In order to use the crate, the libcint should be installed with the outcome library `libcint.so` stored in a reachable path by users.
//!
//! Please visit <https://github.com/sunqm/libcint> for more details about the installation and the usage of libcint
//!
//! The `CINTR2CDATA` struct groups all necessary data for using `libcint`.
//! Various kinds of analytical Gaussian-type orbital (GTO) integrals provided by `libcint` are then wrapped as the methods defined on the `CINTR2CDATA` struct.
//!
//! Currently, only four kinds of integrals are available for both spheric and Cartesian GTOs, including 
//! 1) the one-electron kinetic integral (```CINTR2CDATA::cint_ijkin```),
//! 2) the one-electron nuclear attractive integral (```CINT2CDATA::cint_ijnuc```), 
//! 3) the one-electron overlap integral (```CINTR2CDATA::cint_ijovlp```), 
//! 4) the two-electron repulsive integral (```CINTR2CDATA::cint_ijkl```).
//! The other kinds of integrals are not yet ready in the current version.
//!
//! The integrals:```Vec<f64>``` produced by the CINTR2CDATA methods aformentioned are arranged in
//! the convention of `column-major` matrices according to the definition by `libcint`. 
//!
//! # Examples
//!
//!
//! ```
//! //=============================================================================
//! // Prepare `atm`, `bas` and `env` with the same data structures of those used by `libcint`.
//! // Refer to <https://github.com/sunqm/libcint/blob/master/doc/program_ref.pdf> 
//! // for the details of these data structures.
//! //=============================================================================
//! let mut atm: Vec<Vec<i32>> = vec![];
//! atm.push(vec![2,0,0,0,0,0]);
//! atm.push(vec![4,3,0,0,0,0]);
//! let mut natm = atm.len() as i32;
//! let mut bas: Vec<Vec<i32>> = vec![];
//! bas.push(vec![0,1,1,1,1,6,7,0]);
//! bas.push(vec![0,2,1,1,1,8,9,0]);
//! let mut nbas = bas.len() as i32;
//! let mut env: Vec<f64> = vec![0.0,0.0,0.0,0.7,0.0,0.0,1.0,1.0,0.5,1.0];
//! //=============================================================================
//! // Transfer `atm`, `bas`, and `env` to the raw pointers,
//! // and organize them by the data structure of `CINTR2CDATA`.
//! //=============================================================================
//! use rust_libcint::{CINTR2CDATA,CintType};
//! let mut cint_data = CINTR2CDATA::new();
//! cint_data.initial_r2c(&atm,natm,&bas,nbas,env);
//! //=============================================================================
//! //The 2-electron repulsive integrals (ERIs) for spheric Gaussian-type orbitals
//! //=============================================================================
//! // The GTO functions considered here are spheric.
//! // For Cartesian GTOs, replace `CintType::Spheric` by 
//! //  `CintType::Cartesian` on the following line:
//! cint_data.set_cint_type(&CintType::Spheric);
//! cint_data.cint2e_optimizer_rust();
//! let buf = cint_data.cint_ijkl(0,1,1,0);
//! let mut v1:f64=0.0;
//! &buf.into_iter().for_each(|i| {v1 += i.abs()});
//! println!("The reference data for cint2e ERIs: 0.5745411555937561; v1: {:18.16}; ",v1);
//! //=============================================================================
//! //The one-electron overlap integrals for spheric Gaussian-type orbitals
//! //=============================================================================
//! cint_data.cint_del_optimizer_rust();
//! // The GTO functions considered here are spheric
//! cint_data.set_cint_type(&CintType::Spheric);
//! cint_data.cint1e_ovlp_optimizer_rust();
//! let buf = cint_data.cint_ijovlp(0,1);
//! let mut v1:f64=0.0;
//! &buf.into_iter().for_each(|i| {v1 += i.abs()});
//! println!("The reference data for cint1e_ovlp: 0.7096366827378776; v1: {:18.16}; ",v1);
//! //=============================================================================
//! //The one-electron kinetic integrals for Cartesian Gaussian-type orbitals
//! //=============================================================================
//! cint_data.cint_del_optimizer_rust();
//! // The GTO functions considered here are Cartesian
//! cint_data.set_cint_type(&CintType::Cartesian);
//! cint_data.cint1e_kin_optimizer_rust();
//! let buf = cint_data.cint_ijkin(0,1);
//! let mut v1:f64=0.0;
//! &buf.into_iter().for_each(|i| {v1 += i.abs()});
//! println!("The reference data for cint1e_kin : 1.5780816190296618; v1: {:18.16}; ",v1);
//! //=============================================================================
//! //The one-electron nuclear attraction integrals for Cartesian Gaussian-type orbitals
//! //=============================================================================
//! cint_data.cint_del_optimizer_rust();
//! // The GTO functions considered here are Cartesian
//! cint_data.set_cint_type(&CintType::Cartesian);
//! cint_data.cint1e_nuc_optimizer_rust();
//! let buf = cint_data.cint_ijnuc(0,1);
//! let mut v1:f64=0.0;
//! &buf.into_iter().for_each(|i| {v1 += i.abs()});
//! println!("The reference data for cint1e_nuc : 4.0007622494430706; v1: {:18.16}; ",v1);
//! //=============================================================================
//! // Finally deallocate the memory by transferring the raw pointers back to RUST
//! // i.e. Vec::from_raw_parts();
//! //=============================================================================
//! cint_data.final_c2r();
//! ```

#![allow(unused)]
use std::os::raw::c_int;
use std::mem::ManuallyDrop;

mod cint;
use crate::cint::{CINTOpt,CINTdel_optimizer};

#[derive(Clone,Copy)]
pub enum CintType {
   Spheric,
   Cartesian,
   //Spinor,  // Not yet included
}

pub enum IJOPT {
    Ovlp,
    Kinetic,
    Nuclear,
}

pub enum IJIPOPT {
    IPOvlp,
    IPKin,
    IPNuc,
    IPRInv
}

#[derive(Clone)]
pub struct CINTR2CDATA {
    c_atm: (*const i32, usize, usize),
    c_bas: (*const i32, usize, usize),
    c_env: (*const f64, usize, usize),
    c_nbas: c_int,
    c_natm: c_int,
    c_opt: (*mut CINTOpt, usize, usize),
    c_ao_loc: (*const i32, usize,usize),
    cint_type: CintType,
}

impl CINTR2CDATA {
    /// create a new, empty CINTR2CDATA.
    pub fn new() -> CINTR2CDATA {
        CINTR2CDATA { 
            c_atm: (unsafe {std::ptr::null::<i32>() as *const i32}, 0,0),
            c_bas: (unsafe {std::ptr::null::<i32>() as *const i32}, 0,0),
            c_env: (unsafe {std::ptr::null::<f64>() as *const f64}, 0,0),
            c_ao_loc: (unsafe {std::ptr::null::<i32>() as *const i32}, 0,0),
            c_opt: (unsafe {std::ptr::null::<CINTOpt>() as *mut CINTOpt}, 0,0),
            c_nbas: 0 as c_int,
            c_natm: 0 as c_int,
            cint_type: CintType::Spheric,
            }
    }
    pub fn set_cint_type(&mut self, ctype: &CintType) {
        self.cint_type = *ctype;
    }
    //// 
    pub fn initial_r2c(&mut self, 
                    atm: &Vec<Vec<i32>>, natm:i32, 
                    bas: &Vec<Vec<i32>>, nbas:i32, 
                    env: &Vec<f64>) {
        unsafe {
            let r_atm = Vec::from_raw_parts(self.c_atm.0 as *mut i32,self.c_atm.1,self.c_atm.2);
            let r_bas = Vec::from_raw_parts(self.c_bas.0 as *mut i32,self.c_bas.1,self.c_bas.2);
            let r_env = Vec::from_raw_parts(self.c_env.0 as *mut f64,self.c_env.1,self.c_env.2);
            let r_ao_loc = Vec::from_raw_parts(self.c_ao_loc.0 as *mut i32,self.c_ao_loc.1,self.c_ao_loc.2);
        }

        let dim = bas.iter().map(|ibas| {match self.cint_type {
            CintType::Spheric => {(ibas[1]*2+1)*ibas[3]},
            CintType::Cartesian => {(ibas[1]+1)*(ibas[1]+2)/2*ibas[3]},
        }}).scan(0, |acc,x| {*acc += x; Some(*acc)});

        let mut ao_loc = vec![0; bas.len()+1];
        ao_loc[1..].iter_mut().zip(dim).for_each(|(i,j)| {*i=j});

        ao_loc.shrink_to_fit();
        let mut ao_loc = ManuallyDrop::new(ao_loc);
        //self.c_ao_loc = (ao_loc.as_mut_ptr(), ao_loc.len(), ao_loc.capacity());
        self.c_ao_loc = (ao_loc.as_ptr(), ao_loc.len(), ao_loc.capacity());

        let mut env_f= env.clone();
        env_f.shrink_to_fit();
        let mut env_f = ManuallyDrop::new(env_f);
        self.c_env = (env_f.as_ptr(), env_f.len(), env_f.capacity());

        let mut bas_f= bas.clone().into_iter().flatten().collect::<Vec<i32>>();
        bas_f.shrink_to_fit();
        let mut bas_f = ManuallyDrop::new(bas_f);
        self.c_bas = (bas_f.as_ptr(), bas_f.len(), bas_f.capacity());

        let mut atm_f = atm.clone().into_iter().flatten().collect::<Vec<i32>>();
        atm_f.shrink_to_fit();
        let mut atm_f = ManuallyDrop::new(atm_f);
        self.c_atm = (atm_f.as_ptr(), atm_f.len(), atm_f.capacity());

        self.c_natm = natm as c_int;
        self.c_nbas = nbas as c_int;

        self.c_opt = (unsafe {std::ptr::null::<CINTOpt>() as *mut CINTOpt}, 0,0);
    }
    pub fn final_c2r(&mut self) {
        ///```println!("Clean the unsafe data and transfer the ownership of the raw pointers in CINTR2CDATA to Rust");```
        unsafe {
            let r_atm = Vec::from_raw_parts(self.c_atm.0 as *mut i32,self.c_atm.1,self.c_atm.2);
            let r_bas = Vec::from_raw_parts(self.c_bas.0 as *mut i32,self.c_bas.1,self.c_bas.2);
            let r_env = Vec::from_raw_parts(self.c_env.0 as *mut f64,self.c_env.1,self.c_env.2);
        }
        self.cint_del_optimizer_rust();
        //r_env
    }
    pub fn cint_del_optimizer_rust(&mut self) {
        unsafe{
            CINTdel_optimizer(&mut self.c_opt.0);
        }
    }
    pub fn cint2c2e_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        //self.cint_init_2e_optimizer_rust();
        unsafe {
            cint::cint2c2e_optimizer(&mut self.c_opt.0, 
                                       self.c_atm.0, self.c_natm, 
                                       self.c_bas.0, self.c_nbas, 
                                       self.c_env.0);
        }
    }
    pub fn cint3c2e_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        //self.cint_init_2e_optimizer_rust();
        unsafe {
            cint::cint3c2e_optimizer(&mut self.c_opt.0, 
                                       self.c_atm.0, self.c_natm, 
                                       self.c_bas.0, self.c_nbas, 
                                       self.c_env.0);
        }
    }
    pub fn cint2e_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        //self.cint_init_2e_optimizer_rust();
        unsafe {
            cint::cint2e_optimizer(&mut self.c_opt.0, 
                                       self.c_atm.0, self.c_natm, 
                                       self.c_bas.0, self.c_nbas, 
                                       self.c_env.0);
        }
    }
    pub fn cint1e_ovlp_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        //self.cint_init_optimizer_rust();
        unsafe {
            cint::cint1e_ovlp_optimizer(&mut self.c_opt.0, 
                                       self.c_atm.0, self.c_natm, 
                                       self.c_bas.0, self.c_nbas, 
                                       self.c_env.0);
        }
    }
    pub fn cint1e_nuc_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        //self.cint_init_optimizer_rust();
        unsafe {
            cint::cint1e_nuc_optimizer(&mut self.c_opt.0, 
                                       self.c_atm.0, self.c_natm, 
                                       self.c_bas.0, self.c_nbas, 
                                       self.c_env.0);
        }
    }
    pub fn cint1e_kin_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        //self.cint_init_optimizer_rust();
        unsafe {
            cint::int1e_kin_optimizer(&mut self.c_opt.0, 
                                       self.c_atm.0, self.c_natm, 
                                       self.c_bas.0, self.c_nbas, 
                                       self.c_env.0);
        }
    }
    pub fn int1e_ipovlp_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        unsafe {
            cint::int1e_ipovlp_optimizer(&mut self.c_opt.0, 
                                        self.c_atm.0, self.c_natm, 
                                        self.c_bas.0, self.c_nbas, 
                                        self.c_env.0)
        }
    }
    pub fn int1e_ipkin_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        unsafe {
            cint::int1e_ipkin_optimizer(&mut self.c_opt.0, 
                                        self.c_atm.0, self.c_natm, 
                                        self.c_bas.0, self.c_nbas, 
                                        self.c_env.0)
        }
    }
    pub fn int1e_ipnuc_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        unsafe {
            cint::int1e_ipnuc_optimizer(&mut self.c_opt.0, 
                                        self.c_atm.0, self.c_natm, 
                                        self.c_bas.0, self.c_nbas, 
                                        self.c_env.0)
        }
    }
    pub fn int1e_iprinv_optimizer_rust(&mut self){
        self.cint_del_optimizer_rust();
        unsafe {
            cint::int1e_iprinv_optimizer(&mut self.c_opt.0, 
                                        self.c_atm.0, self.c_natm, 
                                        self.c_bas.0, self.c_nbas, 
                                        self.c_env.0)
        }
    }
    pub fn cint_cgto_rust(&self, index: i32) -> i32 {
        let mut dim: i32;
        unsafe {
            dim = match self.cint_type {
                CintType::Spheric  =>cint::CINTcgto_spheric(index as c_int, self.c_bas.0) as i32,
                CintType::Cartesian=>cint::CINTcgto_cart(index as c_int, self.c_bas.0) as i32,
            };
        }
        dim
    }
    pub fn cint_2c2e(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
        let mut buf: Vec<f64> = vec![0.0;(di*dj) as usize];
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => cint::cint2c2e_sph(c_buf, c_shls,
                                                    self.c_atm.0, self.c_natm,
                                                    self.c_bas.0,self.c_nbas,
                                                    self.c_env.0,
                                                    self.c_opt.0),
                CintType::Cartesian => cint::cint2c2e_cart(c_buf, c_shls,
                                                    self.c_atm.0, self.c_natm,
                                                    self.c_bas.0,self.c_nbas,
                                                    self.c_env.0,
                                                    self.c_opt.0),
            };
            //println!("debug 1 {}", &c_buf.read());
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }
    pub fn cint_3c2e(&mut self, i:i32,j:i32,k:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut dk = self.cint_cgto_rust(k);
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int,k as c_int];
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
        let mut buf: Vec<f64> = vec![0.0;(di*dj*dk) as usize];
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => cint::cint3c2e_sph(c_buf, c_shls,
                                                    self.c_atm.0, self.c_natm,
                                                    self.c_bas.0,self.c_nbas,
                                                    self.c_env.0,
                                                    self.c_opt.0),
                CintType::Cartesian => cint::cint3c2e_cart(c_buf, c_shls,
                                                    self.c_atm.0, self.c_natm,
                                                    self.c_bas.0,self.c_nbas,
                                                    self.c_env.0,
                                                    self.c_opt.0),
            };
            //println!("debug 1 {}", &c_buf.read());
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }
    pub fn cint_ijkl_by_shell(&mut self, i:i32,j:i32,k:i32,l:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut dk = self.cint_cgto_rust(k);
        let mut dl = self.cint_cgto_rust(l);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int,k as c_int,l as c_int];
        //shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = vec![0.0;(di*dj*dk*dl) as usize];
        //buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());

        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric => cint::cint2e_sph(c_buf, c_shls,
                                                    self.c_atm.0, self.c_natm,
                                                    self.c_bas.0,self.c_nbas,
                                                    self.c_env.0,
                                                    self.c_opt.0),
                CintType::Cartesian => cint::cint2e_cart(c_buf, c_shls,
                                                    self.c_atm.0, self.c_natm,
                                                    self.c_bas.0,self.c_nbas,
                                                    self.c_env.0,
                                                    self.c_opt.0),
            };
            //println!("debug 1 {}", &c_buf.read());
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }

    pub fn gto_norm(n:i32,a:f64) -> f64 {
        let mut r: f64 = 0.0_f64;
        unsafe {
            r = cint::CINTgto_norm(n as c_int, a);
        }
        r
    }
    pub fn cint_ij(&mut self, i:i32,j:i32,op_name: &String) -> Vec<f64> {
        // for 1e integrals: ovlp, kinetic, and nuclear
        let op_type = if op_name.to_lowercase() ==String::from("ovlp") {
            IJOPT::Ovlp
        } else if op_name.to_lowercase() ==String::from("kinetic") {
            IJOPT::Kinetic
        } else if op_name.to_lowercase() ==String::from("nuclear") {
            IJOPT::Nuclear
        } else {
            panic!("Error:: Unknown operator for GTO-ij integrals {}", op_name)
        };
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        //shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = [0.0f64].repeat((di*dj) as usize);
        //buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
            match op_type {
                IJOPT::Ovlp => {
                    match self.cint_type {
                        CintType::Spheric => cint::cint1e_ovlp_sph(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                        CintType::Cartesian => cint::cint1e_ovlp_cart(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                    }
                },
                IJOPT::Kinetic => {    
                    match self.cint_type {
                        CintType::Spheric => cint::cint1e_kin_sph(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                        CintType::Cartesian => cint::cint1e_kin_cart(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                    }
                },
                IJOPT::Nuclear => {    
                    match self.cint_type {
                        CintType::Spheric => cint::cint1e_nuc_sph(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                        CintType::Cartesian => cint::cint1e_nuc_cart(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                    }
                },
            };
            let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
        new_buf
    }

    pub fn cint_ip_ij(&mut self, i:i32,j:i32,op_name: &String) -> Vec<f64> {
        // for 1e integrals: ipovlp

        let op_type = if op_name.to_lowercase() == String::from("ipovlp") {
            IJIPOPT::IPOvlp
        } else if op_name.to_lowercase() == String::from("ipkin") {
            IJIPOPT::IPKin
        } else if op_name.to_lowercase() == String::from("ipnuc") {
            IJIPOPT::IPNuc
        } else if op_name.to_lowercase() == String::from("iprinv") {
            IJIPOPT::IPRInv
        } else {
            panic!("Error:: Unknown operator for GTO-ij-ip integrals {}", op_name)
        };
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        //shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = [0.0f64].repeat((3*di*dj) as usize);
        //buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
            match op_type {
                IJIPOPT::IPOvlp => {
                    match self.cint_type {
                        CintType::Spheric => cint::cint1e_ipovlp_sph(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                        CintType::Cartesian => cint::cint1e_ipovlp_cart(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                    }
                },
                IJIPOPT::IPKin => {
                    match self.cint_type {
                        CintType::Spheric => cint::cint1e_ipkin_sph(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                        CintType::Cartesian => cint::cint1e_ipkin_cart(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                    }
                },
                IJIPOPT::IPNuc => {
                    match self.cint_type {
                        CintType::Spheric => cint::cint1e_ipnuc_sph(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                        CintType::Cartesian => cint::cint1e_ipnuc_cart(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                    }
                },
                IJIPOPT::IPRInv => {
                    match self.cint_type {
                        CintType::Spheric => cint::cint1e_iprinv_sph(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                        CintType::Cartesian => cint::cint1e_iprinv_cart(
                                      c_buf, c_shls,
                                        self.c_atm.0, self.c_natm,
                                        self.c_bas.0,self.c_nbas,
                                        self.c_env.0,
                                        self.c_opt.0),
                    }
                }
            };
            //let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
            
            //println!("i={},j={},di={},dj={}", i,j,di,dj);
            //println!("new_buf={:?}", new_buf);
            //println!("new_buf_len={}", new_buf.len());
        }
        new_buf
    }
}

//pub fn cint2e_sph_rust(mut buf: Vec<f64>, mut shls: Vec<i32>, 
//                   c_atm: & *mut c_int, c_natm:c_int, 
//                   c_bas: & *mut c_int, c_nbas:c_int, 
//                   c_env: & *mut f64,
//                   c_opt: & *mut CINTOpt) -> Vec<f64> {
//    //
//    buf.shrink_to_fit();
//    let mut buf = ManuallyDrop::new(buf);
//    let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
//
//    shls.shrink_to_fit();
//    let mut shls = ManuallyDrop::new(shls);
//    let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
//    let mut new_buf:Vec<f64>;
//    unsafe {
//        cint::cint2e_sph(c_buf, c_shls, *c_atm, c_natm, *c_bas, c_nbas, *c_env, *c_opt);
//        //println!("debug 1 {}", &c_buf.read());
//        let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
//        new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
//        //println!("debug 2");
//        //vec![0.0,0.0]
//    }
//    new_buf
//}
//pub fn cint1e_ovlp_sph_rust(mut buf: Vec<f64>, mut shls: Vec<i32>, 
//                   c_atm: & *mut c_int, c_natm:c_int, 
//                   c_bas: & *mut c_int, c_nbas:c_int, 
//                   c_env: & *mut f64,
//                   c_opt: & *mut CINTOpt) -> Vec<f64> {
//    //
//    buf.shrink_to_fit();
//    let mut buf = ManuallyDrop::new(buf);
//    let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
//
//    shls.shrink_to_fit();
//    let mut shls = ManuallyDrop::new(shls);
//    let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
//    let mut new_buf:Vec<f64>;
//    unsafe {
//        cint::cint1e_ovlp_sph(c_buf, c_shls, *c_atm, c_natm, *c_bas, c_nbas, *c_env, *c_opt);
//        let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
//        new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
//    }
//    new_buf
//}
