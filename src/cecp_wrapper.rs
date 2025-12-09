use crate::cecp;
use crate::cecp::ECPOpt;
use std::os::raw::c_int;

pub trait ECPIntegrator {
    unsafe fn optimizer(opt: *mut *mut ECPOpt, atm: *const c_int, natm: c_int, bas: *const c_int, nbas: c_int, env: *const f64);
    unsafe fn integral_sph(
        out: *mut f64,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const ECPOpt,
        cache: *mut f64,
    ) -> c_int;
    unsafe fn integral_cart(
        out: *mut f64,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const ECPOpt,
        cache: *mut f64,
    ) -> c_int;
    fn n_comp() -> usize;
    fn integrator_type() -> &'static str;
    fn name() -> &'static str;
}

macro_rules! impl_ecpintegratorbase {
    (
        $integrator: ident,
        $optimizer: ident,
        $integral_sph: ident,
        $integral_cart: ident,
        $n_comp: expr,
        $integrator_type: literal,
        $name: literal
    ) => {
        #[allow(non_camel_case_types)]
        pub struct $integrator;
        impl ECPIntegrator for $integrator {
            unsafe fn optimizer(opt: *mut *mut ECPOpt, atm: *const c_int, natm: c_int, bas: *const c_int, nbas: c_int, env: *const f64) {
                cecp::$optimizer(opt, atm, natm, bas, nbas, env)
            }
            unsafe fn integral_sph(
                out: *mut f64,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const ECPOpt,
                cache: *mut f64,
            ) -> c_int {
                cecp::$integral_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
            }
            unsafe fn integral_cart(
                out: *mut f64,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const ECPOpt,
                cache: *mut f64,
            ) -> c_int {
                cecp::$integral_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
            }
            fn n_comp() -> usize {
                $n_comp as usize
            }
            fn integrator_type() -> &'static str {
                $integrator_type
            }
            fn name() -> &'static str {
                $name
            }
        }
    };
}

impl_ecpintegratorbase!(ECPscalar, ECPscalar_optimizer, ECPscalar_sph, ECPscalar_cart, 1, "ECP", "ECPIntegratorBase");
impl_ecpintegratorbase!(ECPscalar_ignuc, ECPscalar_ignuc_optimizer, ECPscalar_ignuc_sph, ECPscalar_ignuc_cart, 3, "ECP", "ECPscalar_ignuc");
impl_ecpintegratorbase!(ECPscalar_ipnuc, ECPscalar_ipnuc_optimizer, ECPscalar_ipnuc_sph, ECPscalar_ipnuc_cart, 3, "ECP", "ECPscalar_ipnuc");
impl_ecpintegratorbase!(
    ECPscalar_ipipnuc,
    ECPscalar_ipipnuc_optimizer,
    ECPscalar_ipipnuc_sph,
    ECPscalar_ipipnuc_cart,
    9,
    "ECP",
    "ECPscalar_ipipnuc"
);
impl_ecpintegratorbase!(
    ECPscalar_ipnucip,
    ECPscalar_ipnucip_optimizer,
    ECPscalar_ipnucip_sph,
    ECPscalar_ipnucip_cart,
    9,
    "ECP",
    "ECPscalar_ipnucip"
);
impl_ecpintegratorbase!(
    ECPscalar_iprinv,
    ECPscalar_iprinv_optimizer,
    ECPscalar_iprinv_sph,
    ECPscalar_iprinv_cart,
    3,
    "ECP",
    "ECPscalar_iprinv"
);
impl_ecpintegratorbase!(
    ECPscalar_ipiprinv,
    ECPscalar_ipiprinv_optimizer,
    ECPscalar_ipiprinv_sph,
    ECPscalar_ipiprinv_cart,
    9,
    "ECP",
    "ECPscalar_ipiprinv"
);
impl_ecpintegratorbase!(
    ECPscalar_iprinvip,
    ECPscalar_iprinvip_optimizer,
    ECPscalar_iprinvip_sph,
    ECPscalar_iprinvip_cart,
    9,
    "ECP",
    "ECPscalar_iprinvip"
);
