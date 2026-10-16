//! Semantics of the instance-level optimizer cache: repeated `integrate`
//! calls reuse the cached optimizer, and any molecule mutation invalidates
//! it (correct results after `set_omega` / `set_geom` / `with_*`).

use rest_libcint::prelude::*;
use rest_libcint::test_mol::init_h2o_def2_tzvp;

fn sum_of(out: &[f64]) -> f64 {
    out.iter().sum()
}

#[test]
fn optimizer_cache_reuse_and_invalidation() {
    let mut cint = init_h2o_def2_tzvp();

    // warm the cache: two identical calls must agree exactly
    let (out1, _) = cint.integrate("int1e_ovlp", "s1", None).into();
    let (out2, _) = cint.integrate("int1e_ovlp", "s1", None).into();
    assert_eq!(sum_of(&out1), sum_of(&out2));

    // a mutation (omega) must invalidate: the SR two-electron sum differs from
    // the full one, and repeated SR calls must again be self-consistent
    let (full2e_1, _) = cint.integrate("int2e", "s8", None).into();
    let (full2e_2, _) = cint.integrate("int2e", "s8", None).into();
    assert_eq!(sum_of(&full2e_1), sum_of(&full2e_2));
    cint.set_omega(-1.0);
    let (sr1, _) = cint.integrate("int2e", "s8", None).into();
    let (sr2, _) = cint.integrate("int2e", "s8", None).into();
    assert_eq!(sum_of(&sr1), sum_of(&sr2));
    assert!(
        (sum_of(&sr1) - sum_of(&full2e_1)).abs() > 1.0,
        "omega change must alter the two-electron sum"
    );

    // geometry mutation invalidates as well; restore-verify via with_geom
    let coords: Vec<[f64; 3]> = cint.atom_coords().to_vec();
    let stretched: Vec<[f64; 3]> = coords.iter().map(|c| [c[0] * 1.1, c[1], c[2]]).collect();
    let sum_stretched = cint.with_geom(&stretched, |m| {
        let (o, _) = m.integrate("int1e_ovlp", "s1", None).into();
        sum_of(&o)
    });
    let sum_restored = {
        let (o, _) = cint.integrate("int1e_ovlp", "s1", None).into();
        sum_of(&o)
    };
    assert!(
        (sum_stretched - sum_restored).abs() > 1.0,
        "geometry change must alter the overlap sum"
    );
    assert!(
        (sum_restored - sum_of(&out1)).abs() < 1.0e-10,
        "restored geometry must give back the original overlap"
    );
}
