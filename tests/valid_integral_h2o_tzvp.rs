#[cfg(test)]
mod valid_integral_h2o_tzvp {
    use approx::*;
    use itertools::Itertools;
    use ndarray::prelude::*;
    use num_complex::*;
    use rest_libcint::prelude::*;
    use std::time::Instant;

    #[test]
    fn test_int3c2e_s1_full() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let (out, _) = cint_data.integral_s1::<int3c2e>(None);
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), 5372.255349662842, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), 549.9144320716556, max_relative = 1e-10);
    }

    #[test]
    fn test_int3c2e_s1_slice() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[3, 15], [6, 12], [2, 18]];
        let (out, _) = cint_data.integral_s1::<int3c2e>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), 2061.267953944021, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), 301.27148204542175, max_relative = 1e-10);
    }

    #[test]
    fn test_int3c2e_ip1_s1_slice() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[3, 15], [6, 12], [2, 18]];
        let (out, _) = cint_data.integral_s1::<int3c2e_ip1>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), -708.4938257710874, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), -206.84642420675885, max_relative = 1e-10);
    }

    #[test]
    fn test_int2c2e_ip1_s1_slice() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[3, 15], [6, 12]];
        let (out, _) = cint_data.integral_s1::<int2c2e_ip1>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), -133.92961573356592, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), 13.149333776387465, max_relative = 1e-10);
    }

    #[test]
    fn test_int2e_ip1ip2_s1_slice() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[3, 15], [6, 12], [2, 18], [7, 11]];
        let (out, _) = cint_data.integral_s1::<int2e_ip1ip2>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), 754.4717054716847, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), 153.0307989722583, max_relative = 1e-10);
    }

    #[test]
    fn test_int3c2e_s2ij_full() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let (out, _) = cint_data.integral_s2ij::<int3c2e>(None);
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), 3505.0701294680757, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), 354.77408197530985, max_relative = 1e-10);
    }

    #[test]
    fn test_int3c2e_ip2_s2ij_slice() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[0, 15], [0, 15], [6, 12]];
        let (out, _) = cint_data.integral_s2ij::<int3c2e_ip2>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), 484.24278790302156, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), -85.4382604552415, max_relative = 1e-10);
    }

    #[test]
    fn test_int2e_ip2_s2ij_slice() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[3, 15], [3, 15], [6, 12], [7, 11]];
        let (out, _) = cint_data.integral_s2ij::<int2e_ip2>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), -290.5233235501449, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), -178.07946594201977, max_relative = 1e-10);
    }

    #[test]
    fn test_int2c2e_s2ij_slice() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[3, 15], [3, 15]];
        let (out, _) = cint_data.integral_s2ij::<int2c2e>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), 929.2845254621801, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), 77.98393796393151, max_relative = 1e-10);
    }

    #[test]
    fn test_int3c2e_s2ij_cart_slice() {
        let mut cint_data = initialize();
        cint_data.set_cint_type(&CintType::Cartesian);
        let now = Instant::now();
        let (out, _) = cint_data.integral_s2ij::<int3c2e>(None);
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-1., 1., out.len());
        let out = Array::from_vec(out);
        assert_relative_eq!(out.sum(), 10467.283431217447, max_relative = 1e-10);
        assert_relative_eq!((out * scale).sum(), 245.7706210514359, max_relative = 1e-10);
    }

    #[test]
    fn test_int1e_ignuc_spinor_s1_full() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let (out, out_shape) = cint_data.integral_spinor_s1::<int1e_ignuc>(None);
        println!("out_shape: {:?}", out_shape);
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-2., 1., out.len());
        let out = Array::from_vec(out);
        let res = out.sum();
        assert_relative_eq!(res.re(), 0., epsilon = 1e-10);
        assert_relative_eq!(res.im(), -38.44576060597794, max_relative = 1e-10);
        let res = (out * scale).sum();
        assert_relative_eq!(res.re(), 27.346588674156333, max_relative = 1e-10);
        assert_relative_eq!(res.im(), -66.39923013085044, max_relative = 1e-10);
    }

    #[test]
    fn test_int2e_ip1ip2_spinor_s1_slice() {
        let mut cint_data = initialize();
        let now = Instant::now();
        let shl_slices = vec![[1, 7], [6, 12], [15, 19], [10, 14]];
        let (out, _) = cint_data.integral_spinor_s1::<int2e_ip1ip2>(Some(&shl_slices));
        println!("Elapsed: {:.3?}", now.elapsed());

        let scale = Array::linspace(-2., 1., out.len());
        let out = Array::from_vec(out);
        let res = out.sum();
        assert_relative_eq!(res.re(), -87.88756296443984, max_relative = 1e-10);
        assert_relative_eq!(res.im(), -0.6293582810846077, max_relative = 1e-10);
        let res = (out * scale).sum();
        assert_relative_eq!(res.re(), -30.068275921030978, max_relative = 1e-10);
        assert_relative_eq!(res.im(), -41.555332416720624, max_relative = 1e-10);
    }

    fn initialize() -> CINTR2CDATA {
        // mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5",
        // basis="def2-TZVP").build()
        let c_atm = vec![[8, 20, 1, 23, 0, 0], [1, 24, 1, 27, 0, 0], [1, 28, 1, 31, 0, 0]];
        let c_bas = vec![
            [0, 0, 6, 1, 0, 44, 50, 0],
            [0, 0, 2, 1, 0, 56, 58, 0],
            [0, 0, 1, 1, 0, 60, 61, 0],
            [0, 0, 1, 1, 0, 62, 63, 0],
            [0, 0, 1, 1, 0, 64, 65, 0],
            [0, 1, 4, 1, 0, 66, 70, 0],
            [0, 1, 1, 1, 0, 74, 75, 0],
            [0, 1, 1, 1, 0, 76, 77, 0],
            [0, 2, 1, 1, 0, 78, 79, 0],
            [0, 2, 1, 1, 0, 80, 81, 0],
            [0, 3, 1, 1, 0, 82, 83, 0],
            [1, 0, 3, 1, 0, 32, 35, 0],
            [1, 0, 1, 1, 0, 38, 39, 0],
            [1, 0, 1, 1, 0, 40, 41, 0],
            [1, 1, 1, 1, 0, 42, 43, 0],
            [2, 0, 3, 1, 0, 32, 35, 0],
            [2, 0, 1, 1, 0, 38, 39, 0],
            [2, 0, 1, 1, 0, 40, 41, 0],
            [2, 1, 1, 1, 0, 42, 43, 0],
        ];
        let c_env = vec![
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            1.7763425570911580e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            -4.4476065664656128e-01,
            0.0000000000000000e+00,
            1.7197618551510188e+00,
            0.0000000000000000e+00,
            3.4061340999999999e+01,
            5.1235746000000004e+00,
            1.1646626000000000e+00,
            9.0618446120248586e-01,
            1.6354784928239057e+00,
            2.4145128304249659e+00,
            3.2723041000000003e-01,
            1.0930883523645869e+00,
            1.0307241000000000e-01,
            4.5959135109675275e-01,
            8.0000000000000004e-01,
            2.2072263710762661e+00,
            2.7032382631000000e+04,
            4.0523871392000001e+03,
            9.2232722709999996e+02,
            2.6124070989000001e+02,
            8.5354641350999998e+01,
            3.1035035245000000e+01,
            3.0481181169845928e+00,
            5.6914576328642115e+00,
            9.7338835744432526e+00,
            1.5238733819733028e+01,
            2.0843228934131737e+01,
            2.2391049059992991e+01,
            1.2260860728000001e+01,
            4.9987076005000004e+00,
            1.0568131135849375e+01,
            3.3391469496791393e+00,
            1.1703108158000000e+00,
            2.8427648592056753e+00,
            4.6474740994000002e-01,
            1.4220922112658689e+00,
            1.8504536357000001e-01,
            7.1280983010446131e-01,
            6.3274954801000000e+01,
            1.4627049379000001e+01,
            4.4501223455999996e+00,
            1.5275799646999999e+00,
            6.2570323747894276e+00,
            6.9268656235998423e+00,
            6.0323599265415284e+00,
            3.5035168827833356e+00,
            5.2935117942999999e-01,
            1.3172379939563448e+00,
            1.7478421270000000e-01,
            3.2969483673949351e-01,
            2.3140000000000001e+00,
            1.1328313432935008e+01,
            6.4500000000000002e-01,
            1.2113199965714336e+00,
            1.4279999999999999e+00,
            4.3969226782656516e+00,
        ];

        let c_atm = c_atm.iter().map(|&v| v.to_vec()).collect_vec();
        let c_bas = c_bas.iter().map(|&v| v.to_vec()).collect_vec();
        let mut cint_data = CINTR2CDATA::new();
        cint_data.initial_r2c(&c_atm, c_atm.len() as i32, &c_bas, c_bas.len() as i32, &c_env);
        return cint_data;
    }
}
