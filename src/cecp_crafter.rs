use crate::cecp::*;
use crate::cecp_wrapper::*;
use crate::cint;
use crate::utilities::*;
use crate::{cecp::ECPOpt, CintType, CINTR2CDATA};
use core::ffi::c_int;
use itertools::Itertools;
use rayon::prelude::*;
use rayon::{current_thread_index, max_num_threads};
use std::ptr::{null, null_mut};
use std::sync::Mutex;

#[derive(Clone)]
pub struct ECPData {
    pub atm: Vec<[c_int; 6]>,
    pub bas: Vec<[c_int; 8]>,
    pub env: Vec<f64>,
    pub cint_type: CintType,
    c_opt: *mut ECPOpt,
}

unsafe impl Send for ECPData {}
unsafe impl Sync for ECPData {}

impl ECPData {
    pub fn new() -> ECPData {
        ECPData { atm: Vec::new(), bas: Vec::new(), env: Vec::new(), cint_type: CintType::Spheric, c_opt: null_mut() }
    }

    pub fn from_cint_data(cint_data: &CINTR2CDATA) -> ECPData {
        // initialize
        let mut ecp_data = ECPData {
            c_opt: null_mut(),
            atm: cint_data.atm.clone(),
            bas: cint_data.bas.clone(),
            env: cint_data.env.clone(),
            cint_type: cint_data.cint_type,
        };
        // concate ECP data
        ecp_data.bas.append(&mut cint_data.ecpbas.clone());
        ecp_data.env[AS_ECPBAS_OFFSET as usize] = cint_data.bas.len() as f64;
        ecp_data.env[AS_NECPBAS as usize] = cint_data.ecpbas.len() as f64;
        return ecp_data;
    }

    pub fn optimizer_destruct(&mut self) {
        unsafe {
            ECPdel_optimizer(&mut self.c_opt);
        }
    }

    pub fn optimizer<T>(&mut self)
    where
        T: ECPIntegrator,
    {
        self.optimizer_destruct();
        unsafe {
            T::optimizer(
                &mut self.c_opt,
                self.atm.as_ptr() as _,
                self.atm.len() as c_int,
                self.bas.as_ptr() as _,
                self.bas.len() as c_int,
                self.env.as_ptr(),
            );
        }
    }

    /* #region cgto size */

    /// Size (number of atomic orbitals) of spherical CGTO at certain basis
    /// index.
    ///
    /// Given `None` input, this function will return the maximum size of CGTO,
    /// which may be helpful for allocating output buffer size of
    /// [`Self::integral_block`].
    pub fn cgto_size_sph(&self, id_bas: i32) -> usize {
        use cint::{ANG_OF, BAS_SLOTS, NCTR_OF};
        let loc_ang = ANG_OF as usize;
        let loc_nctr = NCTR_OF as usize;
        let loc_bas = BAS_SLOTS as usize;
        let id_bas = id_bas as usize;
        let nctr = self.bas[id_bas][loc_nctr];
        let nang = self.bas[id_bas][loc_ang] * 2 + 1;
        return (nctr * nang) as usize;
    }

    /// Size (number of atomic orbitals) of cartesian CGTO at certain basis
    /// index.
    ///
    /// Given `None` input, this function will return the maximum size of CGTO,
    /// which may be helpful for allocating output buffer size of
    /// [`Self::integral_block`].
    pub fn cgto_size_cart(&self, id_bas: i32) -> usize {
        use cint::{ANG_OF, BAS_SLOTS, NCTR_OF};
        let loc_ang = ANG_OF as usize;
        let loc_nctr = NCTR_OF as usize;
        let loc_bas = BAS_SLOTS as usize;
        let id_bas = id_bas as usize;
        let nctr = self.bas[id_bas][loc_nctr];
        let l = self.bas[id_bas][loc_ang];
        let nang = (l + 1) * (l + 2) / 2;
        return (nctr * nang) as usize;
    }

    pub fn cgto_size(&self, id_bas: i32) -> usize {
        return match self.cint_type {
            CintType::Spheric => self.cgto_size_sph(id_bas),
            CintType::Cartesian => self.cgto_size_cart(id_bas),
            CintType::Spinor => panic!("Spinor GTOs are not yet supported"),
        };
    }

    /* #endregion */

    /* #region cgto loc */

    /// Location of atomic orbitals, for basis configuration of the current
    /// molecule.
    pub fn cgto_loc(&self) -> Vec<usize> {
        let nbas_without_ecp = self.env[AS_ECPBAS_OFFSET as usize] as i32;
        let size_vec = (0..nbas_without_ecp).map(|idx| self.cgto_size(idx)).collect::<Vec<usize>>();
        let mut loc = vec![0; size_vec.len() + 1];
        for idx in (0..size_vec.len() as usize) {
            loc[idx + 1] = loc[idx] + size_vec[idx];
        }
        return loc;
    }

    /// Location of atomic orbitals, for specified slice of shell.
    pub fn cgto_loc_slice(&self, shl_slice: &[i32; 2]) -> Vec<usize> {
        let loc = self.cgto_loc();
        return loc[shl_slice[0] as usize..shl_slice[1] as usize + 1].to_vec();
    }

    /// Location of atomic orbitals, for specified slices of shell.
    pub fn cgto_loc_slices(&self, shl_slices: &[[i32; 2]]) -> Vec<Vec<usize>> {
        return shl_slices.iter().map(|shl_slice| self.cgto_loc_slice(shl_slice)).collect();
    }

    /// Location of atomic orbitals, relative to the first AO index (start index
    /// to be 0), for specified slice of shell.
    pub fn cgto_loc_slice_relative(&self, shl_slice: &[i32; 2]) -> Vec<usize> {
        let loc_slice = self.cgto_loc_slice(shl_slice);
        return loc_slice.iter().map(|x| x - loc_slice[0]).collect();
    }

    /// Location of atomic orbitals, relative to the first AO index (start index
    /// to be 0), for specified slice of shell.
    pub fn cgto_loc_slices_relative(&self, shl_slices: &[[i32; 2]]) -> Vec<Vec<usize>> {
        return shl_slices.iter().map(|shl_slice| self.cgto_loc_slice_relative(shl_slice)).collect();
    }

    /* #endregion */

    /* #region shl_slices sanity check */

    pub fn check_shl_slices<T>(&self, shl_slices: &[[i32; 2]]) -> Result<(), String>
    where
        T: ECPIntegrator,
    {
        if self.env[AS_NECPBAS as usize] as i32 == 0 {
            return Err("No ECP data found in the molecule.".to_string());
        }
        for shl_slice in shl_slices {
            if shl_slice[1] < shl_slice[0] {
                return Err(format!("in {shl_slice:?}, second number is smaller than first number"));
            }
            if shl_slice[1] > self.bas.len() as i32 {
                return Err(format!("in {:?}, number of shells exceeds {:}", shl_slice, self.bas.len()));
            }
            if shl_slice[0] < 0 {
                return Err(format!("in {shl_slice:?}, shell number should not be negative"));
            }
        }
        return Ok(());
    }

    /* #endregion */

    /* #region cgto shape and buffer size */

    /// Shape of integral (in atomic orbital basis), for specified slices of
    /// shell.
    pub fn cgto_shape<T>(&self, shl_slices: &[[i32; 2]]) -> Vec<usize>
    where
        T: ECPIntegrator,
    {
        self.check_shl_slices::<T>(shl_slices).unwrap();
        let shape = shl_slices
            .iter()
            .map(|shl_slice| {
                let loc = self.cgto_loc_slice(shl_slice);
                loc.last().unwrap().clone() - loc.first().unwrap().clone() as usize
            })
            .collect::<Vec<usize>>();
        return shape;
    }

    /// Obtain cache size for integral.
    ///
    /// This function should be used with the shell slice one desired.
    ///
    /// If the shell slice is not known to you currently (molecule information
    /// has passed into `ECPData`), just pass empty `shls_slice = vec![]`,
    /// then it should give the maximum cache size for this molecule/intor.
    pub fn size_of_cache<T>(&mut self, shls_slice: &[[i32; 2]]) -> usize
    where
        T: ECPIntegrator,
    {
        let shls_min = shls_slice.iter().map(|x| x[0]).min().unwrap_or(0);
        let shls_max = shls_slice.iter().map(|x| x[1]).max().unwrap_or(self.bas.len() as i32);
        let cache_size = (shls_min..shls_max)
            .into_iter()
            .map(|shl| {
                let mut shls = [shl; 4];
                match self.cint_type {
                    CintType::Spheric => unsafe {
                        T::integral_sph(
                            null_mut(),
                            null(),
                            shls.as_mut_ptr(),
                            self.atm.as_ptr() as _,
                            self.atm.len() as c_int,
                            self.bas.as_ptr() as _,
                            self.bas.len() as c_int,
                            self.env.as_ptr(),
                            null(),
                            null_mut(),
                        ) as usize
                    },
                    CintType::Cartesian => unsafe {
                        T::integral_cart(
                            null_mut(),
                            null(),
                            shls.as_mut_ptr(),
                            self.atm.as_ptr() as _,
                            self.atm.len() as c_int,
                            self.bas.as_ptr() as _,
                            self.bas.len() as c_int,
                            self.env.as_ptr(),
                            null(),
                            null_mut(),
                        ) as usize
                    },
                    CintType::Spinor => panic!("Spinor GTOs are not yet supported"),
                }
            })
            .max()
            .unwrap();
        return cache_size;
    }

    // Obtain buffer size for integral.
    pub fn size_of_buffer<T>(&self, shl_slices: &[[i32; 2]]) -> usize
    where
        T: ECPIntegrator,
    {
        self.check_shl_slices::<T>(shl_slices).unwrap();
        T::n_comp()
            * shl_slices
                .iter()
                .map(|&[shl_0, shl_1]| (shl_0..shl_1).map(|shl| self.cgto_size(shl)).max().unwrap() as usize)
                .product::<usize>()
    }

    /* #endregion */

    /// Smallest unit of electron-integral function from libcint.
    ///
    /// This is not a safe wrapper function, though it is pure rust. Use with
    /// caution.
    ///
    /// * `out` - Output integral buffer, need to be allocated enough space
    ///   before calling this function.
    /// * `shls_slice` - shell indices, which size should be n_center; dimension
    ///   check should be performed before calling this function.
    /// * `cache` - Cache buffer, need to be allocated enough space before
    ///   calling this function; simply using `vec![]` should also works, which
    ///   lets libcint manages cache and efficiency decreases. See Also
    ///   [`Self::max_cache_size`] for guide of properly allocate cache.
    pub unsafe fn integral_block<T>(&self, out: &mut [f64], shls: &[i32], shape: &[i32], cache: &mut [f64])
    where
        T: ECPIntegrator,
    {
        let cache_ptr = match cache.len() {
            0 => null_mut(),
            _ => cache.as_mut_ptr(),
        };
        let shape_ptr = match shape.len() {
            0 => null_mut(),
            _ => shape.as_ptr(),
        };
        match self.cint_type {
            CintType::Spheric => unsafe {
                T::integral_sph(
                    out.as_mut_ptr(),
                    shape_ptr,
                    shls.as_ptr(),
                    self.atm.as_ptr() as _,
                    self.atm.len() as c_int,
                    self.bas.as_ptr() as _,
                    self.bas.len() as c_int,
                    self.env.as_ptr(),
                    self.c_opt,
                    cache_ptr,
                )
            },
            CintType::Cartesian => unsafe {
                T::integral_cart(
                    out.as_mut_ptr(),
                    shape_ptr,
                    shls.as_ptr(),
                    self.atm.as_ptr() as _,
                    self.atm.len() as c_int,
                    self.bas.as_ptr() as _,
                    self.bas.len() as c_int,
                    self.env.as_ptr(),
                    self.c_opt,
                    cache_ptr,
                )
            },
            CintType::Spinor => panic!("Spinor GTOs are not yet supported"),
        };
    }

    /// Main integral engine for s1 symmetry.
    ///
    /// This function a low-level API, which is not intended to be called by
    /// user. This function only works for f-contiguous integral (PySCF
    /// convention).
    pub fn integral_s1_inplace<T>(&mut self, out: &mut Vec<f64>, shl_slices: &[[i32; 2]])
    where
        T: ECPIntegrator,
    {
        /* #region 1. dimension definition and sanity check */

        self.check_shl_slices::<T>(shl_slices).unwrap();

        let n_comp = T::n_comp();
        let cgto_shape = self.cgto_shape::<T>(shl_slices);
        let cgto_shape_i32 = cgto_shape.iter().map(|&v| v as i32).collect::<Vec<i32>>();
        let index_shape = shl_slices.iter().map(|[shl_start, shl_stop]| (shl_stop - shl_start) as usize).collect_vec();
        let cgto_locs_rel = self.cgto_loc_slices_relative(shl_slices);

        /* #endregion */

        /* #region 2. preparation for integral engine */

        // optimizer (make integral faster)
        self.optimizer::<T>();

        // cache: thread-local
        let cache_size = self.size_of_cache::<T>(shl_slices);
        let thread_cache = (0..rayon::current_num_threads()).map(|n| Mutex::new(vec![0.; cache_size])).collect_vec();

        // out: enable mut vector by passing immut slice
        let out_const_slice = out.as_slice();

        // reverse iteration for f-contiguous
        let index_shape_rev = index_shape.into_iter().rev().collect_vec();
        let shl_slices_rev = shl_slices.iter().rev().collect_vec();
        let cgto_locs_rel_rev = cgto_locs_rel.iter().rev().collect_vec();
        let cgto_shape_rev = cgto_shape.iter().rev().collect_vec();

        /* #endregion */

        /* #region 3. parallel integral generation */

        // Following code of parallel is not fearless.
        // Variable `out` will be written in parallel, which should be considered
        // racing, when calling `self.integral_block::<T>`,
        // and racing would not actually happen if I am careful.

        (0..(index_shape_rev[0] * index_shape_rev[1])).into_par_iter().for_each(|idx_01| {
            let idx_0 = idx_01 / index_shape_rev[1];
            let idx_1 = idx_01 % index_shape_rev[1];
            let shl_0 = idx_0 as i32 + shl_slices_rev[0][0];
            let shl_1 = idx_1 as i32 + shl_slices_rev[1][0];
            let cgto_0 = cgto_locs_rel_rev[0][idx_0];
            let cgto_1 = cgto_locs_rel_rev[1][idx_1];

            let thread_index = current_thread_index().unwrap_or(0);
            let mut cache = thread_cache[thread_index].lock().unwrap();

            let shls = [shl_1, shl_0];
            let offset = cgto_1 + cgto_shape_rev[1] * cgto_0;

            unsafe {
                let out_with_offset = cast_mut_slice(&out_const_slice[offset..]);
                self.integral_block::<T>(out_with_offset, &shls, &cgto_shape_i32, &mut cache);
            }
        });
        /* #endregion */

        /* #region 4. cleanup */
        self.optimizer_destruct();
        /* #endregion */
    }

    pub fn integral_s1<T>(&mut self, shl_slices: Option<&[[i32; 2]]>) -> (Vec<f64>, Vec<usize>)
    where
        T: ECPIntegrator,
    {
        let shl_slices = match shl_slices {
            Some(shl_slices) => shl_slices,
            None => {
                let nbas = self.env[AS_ECPBAS_OFFSET as usize] as i32;
                &vec![[0, nbas]; 2]
            },
        };
        let mut out_shape = self.cgto_shape::<T>(&shl_slices);
        if T::n_comp() > 1 {
            out_shape.push(T::n_comp());
        }
        let out_size = out_shape.iter().product::<usize>();
        let mut out = Vec::<f64>::with_capacity(out_size);
        unsafe { out.set_len(out_size) };
        self.integral_s1_inplace::<T>(&mut out, &shl_slices);
        return (out, out_shape);
    }
}

impl CINTR2CDATA {
    pub fn integral_ecp_s1<T>(&mut self, shl_slices: Option<&[[i32; 2]]>) -> (Vec<f64>, Vec<usize>)
    where
        T: ECPIntegrator,
    {
        if self.ecpbas.len() == 0 {
            panic!("No ECP data found in the molecule.");
        }
        let mut ecp_data = ECPData::from_cint_data(&self);
        ecp_data.integral_s1::<T>(shl_slices)
    }
}
