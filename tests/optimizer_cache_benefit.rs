//! Decisive A/B of the optimizer cache, measured **inside one process** so
//! machine drift cancels:
//!
//! - A ("old"): `CInt::integrate_cross` per batch — concatenates a fresh
//!   molecule every call, so its optimizer cache is always cold and libcint's
//!   O(nbas²) optimizer is rebuilt per call (the pre-cache behaviour).
//! - B ("new"): one persistent mol+aux concatenation, per-batch
//!   `merged.integrate` — the instance cache hits after the first batch.
//!
//! The batch structure mimics the analytical-gradient de-jk loop: the aux
//! shells are split into batches, each evaluating `int3c2e_ip1` (s1) and
//! `int3c2e_ip2` (s2ij). Run with:
//!
//! ```text
//! cargo test --release --test optimizer_cache_benefit -- --nocapture
//! ```

use rest_libcint::prelude::*;
use rest_libcint::test_mol::{init_h2o_def2_jk, init_h2o_def2_tzvp};

#[test]
fn optimizer_cache_benefit() {
    let h2o = init_h2o_def2_tzvp();
    let h2o_jk = init_h2o_def2_jk();

    // ~600-700 merged shells: (H2O)_15 def2-TZVP + universal-jkfit, close to
    // the 45-atom salt cluster's mol+aux shell count.
    let n = 15;
    let mut mol = h2o.clone();
    let mut aux = h2o_jk.clone();
    for _ in 1..n {
        mol = &mol + &h2o;
        aux = &aux + &h2o_jk;
    }
    let nbas_mol = mol.nbas();
    let nbas_aux = aux.nbas();
    let merged = &mol + &aux;
    println!(
        "mol shells = {}, aux shells = {}, merged shells = {}",
        nbas_mol, nbas_aux, merged.nbas()
    );

    // ~25 aux batches, as the memory-driven partition of a gradient would cut
    let nbatch = 25.min(nbas_aux);
    let per = (nbas_aux + nbatch - 1) / nbatch;
    let batches: Vec<[usize; 2]> = (0..nbatch)
        .map(|b| {
            let s = (b * per).min(nbas_aux);
            let e = ((b + 1) * per).min(nbas_aux);
            [s, e]
        })
        .collect();

    let run_old = |intor: &str, symm: &str, batches: &[[usize; 2]]| -> f64 {
        let t = std::time::Instant::now();
        for &[shl0, shl1] in batches {
            let slc = [[0, nbas_mol], [0, nbas_mol], [shl0, shl1]];
            let _ = CInt::integrate_cross(intor, [&mol, &mol, &aux], symm, slc);
        }
        t.elapsed().as_secs_f64()
    };
    let run_new = |intor: &str, symm: &str, batches: &[[usize; 2]]| -> f64 {
        let t = std::time::Instant::now();
        for &[shl0, shl1] in batches {
            let slc = [[0, nbas_mol], [0, nbas_mol], [nbas_mol + shl0, nbas_mol + shl1]];
            let _ = merged.integrate(intor, symm, slc);
        }
        t.elapsed().as_secs_f64()
    };

    // interleaved A/B, 3 rounds each, both intors
    let mut old_ip1 = Vec::new();
    let mut new_ip1 = Vec::new();
    let mut old_ip2 = Vec::new();
    let mut new_ip2 = Vec::new();
    for _ in 0..3 {
        old_ip1.push(run_old("int3c2e_ip1", "s1", &batches));
        new_ip1.push(run_new("int3c2e_ip1", "s1", &batches));
        old_ip2.push(run_old("int3c2e_ip2", "s2ij", &batches));
        new_ip2.push(run_new("int3c2e_ip2", "s2ij", &batches));
    }
    let med = |v: &Vec<f64>| {
        let mut s = v.clone();
        s.sort_by(|a, b| a.partial_cmp(b).unwrap());
        s[s.len() / 2]
    };
    let (o1, n1) = (med(&old_ip1), med(&new_ip1));
    let (o2, n2) = (med(&old_ip2), med(&new_ip2));
    println!(
        "ip1 旧(每批重建) {:.3} s vs 新(缓存) {:.3} s  -> 省 {:.1} ms/批",
        o1, n1, (o1 - n1) * 1.0e3 / nbatch as f64
    );
    println!(
        "ip2 旧(每批重建) {:.3} s vs 新(缓存) {:.3} s  -> 省 {:.1} ms/批",
        o2, n2, (o2 - n2) * 1.0e3 / nbatch as f64
    );
    println!(
        "合计每 force 当量 ({} 批 × 2 积分器): 旧 {:.2} s vs 新 {:.2} s, 省 {:.2} s",
        nbatch,
        o1 + o2,
        n1 + n2,
        (o1 + o2) - (n1 + n2)
    );
}
