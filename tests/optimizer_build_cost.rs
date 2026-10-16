//! Micro-benchmark: cost of building libcint optimizers (`CInt::get_optimizer`)
//! for representative mol+aux concatenations.
//!
//! The unscreened `integrate*` drivers rebuild the optimizer on every call; this
//! benchmark quantifies the per-build cost that an optimizer cache would
//! amortize. Run with:
//!
//! ```text
//! cargo test --release --test optimizer_build_cost -- --nocapture
//! ```

use rest_libcint::prelude::*;
use rest_libcint::test_mol::{init_h2o_def2_jk, init_h2o_def2_tzvp, init_c10h22_def2_qzvp};
use rest_libcint::CInt;

fn total_primitives(cint: &CInt) -> usize {
    // bas: [ [c_int; 8] ], NPRIM_OF = index 1 per libcint bas layout
    cint.bas.iter().map(|b| b[1] as usize).sum()
}

fn bench(cint: &CInt, label: &str) {
    let intors = ["int3c2e", "int3c2e_ip1", "int3c2e_ip2", "int2e"];
    println!(
        "== {} : nbas(shells) = {}, nao = {}, tot_prim = {} (setij guard at 2048)",
        label,
        cint.nbas(),
        cint.ao_loc().last().unwrap(),
        total_primitives(cint),
    );
    for intor in intors {
        let integrator = CInt::get_integrator(intor);
        // warm-up + measure 5 builds
        let mut times = Vec::with_capacity(5);
        for _ in 0..5 {
            let t = std::time::Instant::now();
            let _opt = cint.get_optimizer(&*integrator);
            times.push(t.elapsed().as_secs_f64() * 1.0e3);
        }
        let _ = CInt::get_integrator(intor); // drop integrator after opts
        println!(
            "   {:<14} builds (ms): {:?}",
            intor,
            times.iter().map(|t| format!("{:.2}", t)).collect::<Vec<_>>()
        );
    }
}

#[test]
fn optimizer_build_cost() {
    let h2o = init_h2o_def2_tzvp();
    let h2o_jk = init_h2o_def2_jk();

    // (H2O)_10 def2-TZVP + universal-jkfit: merge 10 copies of each
    let mut mol10 = h2o.clone();
    for _ in 0..9 {
        mol10 = &mol10 + &h2o;
    }
    let mut aux10 = h2o_jk.clone();
    for _ in 0..9 {
        aux10 = &aux10 + &h2o_jk;
    }
    let merged10 = &mol10 + &aux10;

    // C10H22 def2-QZVP (large single molecule, no aux merge needed for cost probing)
    let c10 = init_c10h22_def2_qzvp();

    bench(&merged10, "(H2O)10 def2-TZVP + universal-jkfit (merged)");
    bench(&mol10, "(H2O)10 def2-TZVP (mol only)");
    bench(&c10, "C10H22 def2-QZVP");
}
