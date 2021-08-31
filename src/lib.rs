//! # rust-libcint
//!
//! The `rust-libcint` crate provides wrappers for libcint (C).
//!
//! In order to use the crate correctly, the libcint should be installed with the outcome library `libcint.so` stored in a reachable path.
//!
//! Please visit <https://github.com/sunqm/libcint> for more details about the installation and the usage of libcint
//!
#![allow(unused)]
use std::os::raw::c_int;
//use std::os::raw::c_double;
//use std::ffi::c_void;
use std::ptr;
use std::mem::ManuallyDrop;

mod cint;
use crate::cint::{CINTOpt,CINTdel_optimizer};

pub enum CintType {
   Spheric,
   Cartesian,
   //Spinor,  // Not yet included
}

pub struct CINTR2CDATA {
    c_atm: (*mut i32, usize, usize),
    c_bas: (*mut i32, usize, usize),
    c_env: (*mut f64, usize, usize),
    c_nbas: c_int,
    c_natm: c_int,
    c_opt: (*mut CINTOpt, usize, usize),
    cint_type: CintType,
}

/// A new type that groups all necessary data to communicate with libcint
///
/// # Examples
///
/// ```
/// use rust_libcint::{CINT2CDATA,CintType};
/// let mut atm: Vec<Vec<i32>> = vec![];
/// atm.push(vec![2,0,0,0,0,0]);
/// let mut bas: Vec<Vec<i32>> = vec![];
/// bas.push(vec![0,1,1,1,1,3,4,0]);
/// bas.push(vec![0,2,1,1,1,5,6,0]);
/// let mut env: Vec<f64> = vec![0.0,0.0,0.0,2.0,1.0,1.5,0.5];
/// let mut cint_data = rust_libcint::CINTR2CDATA::new();
/// cint_data.initial_r2c(&atm,1,&bas,2,env);
/// ```
/// 2-electron analytic integrals for the shperic Gaussian-type orbitals
/// ```
/// cint_data.set_cint_type(CintType::Spheric);
/// cint_data.cint2e_optimizer_rust();
/// let buf = cint_data.cint_ijkl(0,1,1,0);
/// let mut v1:f64=0.0;
/// println!("{:?}",&buf);
/// &buf.into_iter().for_each(|i| {v1 += i.abs()});
/// println!("v1: {}",v1);
/// ```
/// the overlap integrals for the shperic Gaussian-type orbitals
/// ```
/// cint_data.cint_del_optimizer_rust();
/// cint_data.set_cint_type(CintType::Spheric);
/// cint_data.cint1e_ovlp_optimizer_rust();
/// let buf = cint_data.cint_ijovlp(0,1);
/// let mut v1:f64=0.0;
/// println!("{:?}",&buf);
/// &buf.into_iter().for_each(|i| {v1 += i.abs()});
/// println!("v1: {}",v1);
/// ```
/// the kinetic integrals for the cartesian Gaussian-type orbitals
/// ```
/// cint_data.cint_del_optimizer_rust();
/// cint_data.set_cint_type(CintType::Cartesian);
/// cint_data.cint1e_kin_optimizer_rust();
/// let buf = cint_data.cint_ijkin(0,1);
/// let mut v1:f64=0.0;
/// println!("{:?}",&buf);
/// &buf.into_iter().for_each(|i| {v1 += i.abs()});
/// println!("v1: {}",v1);
/// ```
impl CINTR2CDATA {
    /// create a new, empty CINTR2CDATA.
    pub fn new() -> CINTR2CDATA {
        CINTR2CDATA { 
            c_atm: (unsafe {std::ptr::null::<i32>() as *mut i32}, 0,0),
            c_bas: (unsafe {std::ptr::null::<i32>() as *mut i32}, 0,0),
            c_env: (unsafe {std::ptr::null::<f64>() as *mut f64}, 0,0),
            c_opt: (unsafe {std::ptr::null::<CINTOpt>() as *mut CINTOpt}, 0,0),
            c_nbas: 0 as c_int,
            c_natm: 0 as c_int,
            cint_type: CintType::Spheric,
            }
    }
    pub fn set_cint_type(&mut self, ctype: CintType) {
        self.cint_type = ctype;
    }
    //// 
    pub fn initial_r2c(&mut self, 
                    atm: &Vec<Vec<i32>>, natm:i32, 
                    bas: &Vec<Vec<i32>>, nbas:i32, 
                    mut env: Vec<f64>) {
        unsafe {
            let r_atm = Vec::from_raw_parts(self.c_atm.0,self.c_atm.1,self.c_atm.2);
            let r_bas = Vec::from_raw_parts(self.c_bas.0,self.c_bas.1,self.c_bas.2);
            let r_env = Vec::from_raw_parts(self.c_env.0,self.c_env.1,self.c_env.2);
        }
        env.shrink_to_fit();
        let mut env = ManuallyDrop::new(env);
        self.c_env = (env.as_mut_ptr(), env.len(), env.capacity());

        let mut bas_f:  Vec<i32> = vec![];
        &mut (0..bas.len()).for_each(|i| {
            bas_f.extend(bas[i].clone());
        });
        bas_f.shrink_to_fit();
        let mut bas_f = ManuallyDrop::new(bas_f);
        self.c_bas = (bas_f.as_mut_ptr(), bas_f.len(), bas_f.capacity());

        let mut atm_f:  Vec<i32> = vec![];
        &mut (0..atm.len()).for_each(|i| {
            atm_f.extend(atm[i].clone());
        });
        atm_f.shrink_to_fit();
        let mut atm_f = ManuallyDrop::new(atm_f);
        self.c_atm = (atm_f.as_mut_ptr(), atm_f.len(), atm_f.capacity());

        self.c_natm = natm as c_int;
        self.c_nbas = nbas as c_int;

        self.c_opt = (unsafe {std::ptr::null::<CINTOpt>() as *mut CINTOpt}, 0,0);
    }
    pub fn final_c2r(&mut self) {
        println!("Transfer the ownership of the raw pointers in CINTR2CDATA to Rust");
        unsafe {
            let r_atm = Vec::from_raw_parts(self.c_atm.0,self.c_atm.1,self.c_atm.2);
            let r_bas = Vec::from_raw_parts(self.c_bas.0,self.c_bas.1,self.c_bas.2);
            let r_env = Vec::from_raw_parts(self.c_env.0,self.c_env.1,self.c_env.2);
        }
        self.cint_del_optimizer_rust();
    }
    pub fn cint_del_optimizer_rust(&mut self) {
        unsafe{
            CINTdel_optimizer(&mut self.c_opt.0);
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
    pub fn cint_ijkl(&mut self, i:i32,j:i32,k:i32,l:i32) -> Vec<f64> {
        let mut di = self.cint_cgto_rust(i);
        let mut dj = self.cint_cgto_rust(j);
        let mut dk = self.cint_cgto_rust(k);
        let mut dl = self.cint_cgto_rust(l);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int,k as c_int,l as c_int];
        shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = [0.0f64].repeat((di*dj*dk*dl) as usize);
        buf.shrink_to_fit();
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
            let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
       new_buf
    }
    pub fn cint_ijovlp(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = [0.0f64].repeat((di*dj) as usize);
        buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
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
            };
            let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
        new_buf
    }
    pub fn cint_ijnuc(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = [0.0f64].repeat((di*dj) as usize);
        buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
            match self.cint_type {
                CintType::Spheric =>  cint::cint1e_nuc_sph(
                           c_buf, c_shls,
                             self.c_atm.0, self.c_natm,
                             self.c_bas.0,self.c_nbas,
                             self.c_env.0,
                             self.c_opt.0),
                CintType::Cartesian =>  cint::cint1e_nuc_cart(
                           c_buf, c_shls,
                             self.c_atm.0, self.c_natm,
                             self.c_bas.0,self.c_nbas,
                             self.c_env.0,
                             self.c_opt.0),
            };
            let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
        new_buf
    }

    pub fn cint_ijkin(&mut self, i:i32,j:i32) -> Vec<f64> {
        let mut di: i32 = self.cint_cgto_rust(i);
        let mut dj: i32 = self.cint_cgto_rust(j);
    
        let mut shls: Vec<c_int> = vec![i as c_int,j as c_int];
        shls.shrink_to_fit();
        let mut shls = ManuallyDrop::new(shls);
        let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    
        let mut buf: Vec<f64> = [0.0f64].repeat((di*dj) as usize);
        buf.shrink_to_fit();
        let mut buf = ManuallyDrop::new(buf);
        let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());
    
        let mut new_buf:Vec<f64>;
        unsafe {
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
            };
            let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
            new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        }
        new_buf
    }
}

pub fn cint2e_sph_rust(mut buf: Vec<f64>, mut shls: Vec<i32>, 
                   c_atm: & *mut c_int, c_natm:c_int, 
                   c_bas: & *mut c_int, c_nbas:c_int, 
                   c_env: & *mut f64,
                   c_opt: & *mut CINTOpt) -> Vec<f64> {
    //
    buf.shrink_to_fit();
    let mut buf = ManuallyDrop::new(buf);
    let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());

    shls.shrink_to_fit();
    let mut shls = ManuallyDrop::new(shls);
    let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    let mut new_buf:Vec<f64>;
    unsafe {
        cint::cint2e_sph(c_buf, c_shls, *c_atm, c_natm, *c_bas, c_nbas, *c_env, *c_opt);
        //println!("debug 1 {}", &c_buf.read());
        let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
        new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
        //println!("debug 2");
        //vec![0.0,0.0]
    }
    new_buf
}
pub fn cint1e_ovlp_sph_rust(mut buf: Vec<f64>, mut shls: Vec<i32>, 
                   c_atm: & *mut c_int, c_natm:c_int, 
                   c_bas: & *mut c_int, c_nbas:c_int, 
                   c_env: & *mut f64,
                   c_opt: & *mut CINTOpt) -> Vec<f64> {
    //
    buf.shrink_to_fit();
    let mut buf = ManuallyDrop::new(buf);
    let (c_buf, buf_len, buf_cap) = (buf.as_mut_ptr() as *mut f64, buf.len(), buf.capacity());

    shls.shrink_to_fit();
    let mut shls = ManuallyDrop::new(shls);
    let (c_shls,shls_len,shls_cap) = (shls.as_mut_ptr() as *mut c_int,shls.len(),shls.capacity());
    let mut new_buf:Vec<f64>;
    unsafe {
        cint::cint1e_ovlp_sph(c_buf, c_shls, *c_atm, c_natm, *c_bas, c_nbas, *c_env, *c_opt);
        let shls = Vec::from_raw_parts(c_shls, shls_len, shls_cap);
        new_buf = Vec::from_raw_parts(c_buf, buf_len, buf_cap);
    }
    new_buf
}
