// // Integration tests related to geodesic.

// extern crate utilities;

// use geographiclib_rs::{Geodesic, geodesiccapability as caps};
// use std::sync::{Arc, Mutex};
// use utilities::{assert_delta, nC_, DeltaEntry, entries[].add, test_basic};

// // placeholder: Geodesic_A1m1f
// // placeholder: Geodesic_A2m1f
// // placeholder: Geodesic_A3coeff
// // placeholder: Geodesic_A3f

// // Note: For Geodesic_Astroid, see geomath_tests.rs

// // placeholder: Geodesic_C1f
// // placeholder: Geodesic_C1pf
// // placeholder: Geodesic_C2f
// // placeholder: Geodesic_C3coeff
// // placeholder: Geodesic_C3f
// // placeholder: Geodesic_C4coeff
// // placeholder: Geodesic_C4f

// // #[test]
// // fn test_consts_vs_cpp() { // todo: enable (but ALL checks need non-pub access as of 2020/12/31)
// //     // Format: nA1_ nC1_ nC1p_ nA2_ nC2_ nA3_ nA3x_ nC3_ nC3x_ nC4_ nC4x_ nC_ maxit1_ GEOGRAPHICLIB_GEODESIC_ORDER
// //     // todo: review commented-out tests to see if there's an equivalent const or simple calculation we can compare
// //     let items = utilities::read_consts_basic("Geodesic_consts", 14);
// //     let line_num = 2; // todo: remove once assert_delta stops requiring a line number
// //     assert_delta!(items[0], Geodesic::nA1_, 0.0, false, "nA1_", line_num);
// //     assert_delta!(items[1], Geodesic::nC1_, 0.0, false, "nC1_", line_num);
// //     assert_delta!(items[2], Geodesic::nC1p_, 0.0, false, "nC1p_", line_num);
// //     assert_delta!(items[3], Geodesic::nA2_, 0.0, false, "nA2_", line_num);
// //     assert_delta!(items[4], Geodesic::nC2_, 0.0, false, "nC2_", line_num);
// //     assert_delta!(items[5], Geodesic::nA3_, 0.0, false, "nA3_", line_num);
// //     assert_delta!(items[6], Geodesic::nA3x_, 0.0, false, "nA3x_", line_num);
// //     assert_delta!(items[7], Geodesic::nC3_, 0.0, false, "nC3_", line_num);
// //     assert_delta!(items[8], Geodesic::nC3x_, 0.0, false, "nC3x_", line_num);
// //     assert_delta!(items[9], Geodesic::nC4_, 0.0, false, "nC4_", line_num);
// //     assert_delta!(items[10], Geodesic::nC4x_, 0.0, false, "nC4x_", line_num);
// //     assert_delta!(items[11], Geodesic::nC_, 0.0, false, "nC_", line_num);
// //     assert_delta!(items[12], Geodesic::maxit1_ as f64, 0.0, false, "maxit1_", line_num);
// //     assert_delta!(items[13], Geodesic::GEODESIC_ORDER, 0.0, false, "GEODESIC_ORDER", line_num);
// // }

// #[test]
// fn test_gen_direct_vs_cpp() {
//     // Format: this-in[_a _f] lat1 lon1 azi1 arcmode s12_a12 outmask result=a12 lat2-out lon2-out azi2-out s12-out m12-out M12-out M21-out S12-out

//     let delta_entries = Arc::new(Mutex::new(vec![DeltaEntry::new(); 9]));
//     test_basic("Geodesic_GenDirect", 17, |line_num, items| {
//         let g = Geodesic::new(items[0], items[1]);
//         let outmask = items[7] as u64;
//         let result = g._gen_direct(items[2], items[3], items[4], items[5] != 0.0, items[6], outmask);
//         let mut entries = delta_entries.lock().unwrap();

//         entries[].add(items[8], result.0, line_num, &mut entries[0]);
//         // todo: reduce allowed difference to 0.0
//         assert_delta!(items[8], result.0, 1e-10, false, "result (a12 or result.0)", line_num);

//         if outmask & caps::LATITUDE != 0 {
//             entries[].add(items[9], result.1, line_num, &mut entries[1]);
//             // todo: reduce allowed difference to 0.0
//             // todo: disallow sign difference
//             assert_delta!(items[9], result.1, 1e1, true, "lat2-out (result.1)", line_num);
//         }
//         if outmask & caps::LONGITUDE != 0 {
//             entries[].add(items[10], result.2, line_num, &mut entries[2]);
//             // todo: reduce allowed difference to 0.0
//             // todo: disallow sign difference
//             assert_delta!(items[10], result.2, 1e1, true, "lon2-out (result.2)", line_num);
//         }
//         if outmask & caps::AZIMUTH != 0 {
//             entries[].add(items[11], result.3, line_num, &mut entries[3]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[11], result.3, 1e1, false, "azi2-out (result.3)", line_num);
//         }
//         if outmask & caps::DISTANCE != 0 {
//             entries[].add(items[12], result.4, line_num, &mut entries[4]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[12], result.4, 1e-2, false, "s12-out (result.4)", line_num);
//         }
//         if outmask & caps::REDUCEDLENGTH != 0 {
//             // todo: reduce allowed difference to 0.0
//             // todo: disallow sign difference
//             // todo: remove exceptions
//             entries[].add(items[13], result.5, line_num, &mut entries[5]);
//             if !f64::is_nan(result.5) {
//                 assert_delta!(items[13], result.5, 1e1, true, "m12-out (result.5)", line_num);
//             }
//         }
//         if outmask & caps::GEODESICSCALE != 0 {
//             // todo: reduce allowed difference to 0.0
//             // todo: disallow sign difference
//             // todo: remove exceptions
//             entries[].add(items[14], result.6, line_num, &mut entries[6]);
//             if !f64::is_nan(result.6) {
//                 assert_delta!(items[14], result.6, 1e1, true, "M12-out (result.6)", line_num);
//             }
//             // todo: reduce allowed difference to 0.0
//             // todo: disallow sign difference
//             // todo: remove exceptions
//             entries[].add(items[15], result.7, line_num, &mut entries[7]);
//             if !f64::is_nan(result.7) {
//                 assert_delta!(items[15], result.7, 1e1, true, "M21-out (result.7)", line_num);
//             }
//         }
//         if outmask & caps::AREA != 0 {
//             // todo: reduce allowed difference to 0.0
//             // todo: disallow sign difference
//             // todo: remove exceptions
//             entries[].add(items[16], result.8, line_num, &mut entries[8]);
//             if !f64::is_infinite(result.8) {
//                 assert_delta!(items[16], result.8, 1e1, true, "S12-out (result.8)", line_num);
//             }
//         }
//     });

//     delta_entries.lock().unwrap().iter().enumerate().for_each(|(i, entry)| {
//         println!("test_gen_direct_vs_cpp result.{} abs {:e} at line {}, rel {:e} at line {}, sign-change {}", i, entry.diff_abs, entry.line_abs, entry.diff_rel, entry.line_rel, entry.sign_change);
//     });
// }

// // placeholder: Geodesic_GenDirectLine

// #[test]
// fn test_gen_inverse_azi_vs_cpp() {
//     // Format: this-in[_a _f] lat1 lon1 lat2 lon2 outmask result=a12 s12-out azi1-out azi2-out m12-out M12-out M21-out S12-out
//     let delta_entries = Arc::new(Mutex::new(vec![DeltaEntry::new(); 8]));
//     test_basic("Geodesic_GenInverse_out7", 15, |line_num, items| {
//         let g = Geodesic::new(items[0], items[1]);
//         let outmask = items[6] as u64;
//         let result = g._gen_inverse_azi(items[2], items[3], items[4], items[5], outmask);
//         let mut entries = delta_entries.lock().unwrap();
//         entries[].add(items[7], result.0, line_num, &mut entries[0]);
//         // todo: reduce allowed difference to 0.0
//         assert_delta!(items[7], result.0, 1e-15, false, "result (a12 or result.0)", line_num);
//         if outmask & caps::DISTANCE != 0 {
//             entries[].add(items[8], result.1, line_num, &mut entries[1]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[8], result.1, 1e-15, false, "s12-out (result.1)", line_num);
//         }
//         if outmask & caps::AZIMUTH != 0 {
//             entries[].add(items[9], result.2, line_num, &mut entries[2]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[9], result.2, 1e-11, false, "azi1-out (result.2)", line_num);

//             entries[].add(items[10], result.3, line_num, &mut entries[3]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[10], result.3, 1e-8, false, "azi2-out (result.3)", line_num);
//         }
//         if outmask & caps::REDUCEDLENGTH != 0 {
//             entries[].add(items[11], result.4, line_num, &mut entries[4]);
//             // todo: reduce allowed difference to 0.0
//             // todo: remove exceptions
//             if !result.4.is_nan() {
//                 assert_delta!(items[11], result.4, 1e-10, false, "m12-out (result.4)", line_num);
//             }
//         }
//         if outmask & caps::GEODESICSCALE != 0 {
//             entries[].add(items[12], result.5, line_num, &mut entries[5]);
//             // todo: reduce allowed difference to 0.0
//             // todo: remove exceptions
//             if !result.5.is_nan() {
//                 assert_delta!(items[12], result.5, 1e-14, false, "M12-out (result.5)", line_num);
//             }

//             entries[].add(items[13], result.6, line_num, &mut entries[6]);
//             // todo: reduce allowed difference to 0.0
//             // todo: remove exceptions
//             if !result.6.is_nan() {
//                 assert_delta!(items[13], result.6, 1e-14, false, "M21-out (result.6)", line_num);
//             }
//         }
//         if outmask & caps::AREA != 0 {
//             entries[].add(items[14], result.7, line_num, &mut entries[7]);
//             // todo: reduce allowed difference to 0.0
//             // todo: remove exceptions
//             if !result.7.is_nan() && !result.7.is_infinite() {
//                 assert_delta!(items[14], result.7, 1e-7, false, "S12-out (result.7)", line_num);
//             }
//         }
//     });

//     delta_entries.lock().unwrap().iter().enumerate().for_each(|(i, entry)| {
//         println!("test_gen_inverse_azi_vs_cpp result.{} abs {:e} at line {}, rel {:e} at line {}, sign-change {}", i, entry.diff_abs, entry.line_abs, entry.diff_rel, entry.line_rel, entry.sign_change);
//     });
// }

// #[test]
// fn test_gen_inverse_vs_cpp() {
//     // Format: this-in[_a _f] lat1 lon1 lat2 lon2 outmask result=a12 s12-out salp1-out calp1-out salp2-out calp2-out m12-out M12-out M21-out S12-out
//     let delta_entries = Arc::new(Mutex::new(vec![DeltaEntry::new(); 10]));
//     test_basic("Geodesic_GenInverse_out9", 17, |line_num, items| {
//         let g = Geodesic::new(items[0], items[1]);
//         let outmask = items[6] as u64;
//         let result = g._gen_inverse(items[2], items[3], items[4], items[5], outmask);
//         let mut entries = delta_entries.lock().unwrap();
//         entries[].add(items[7], result.0, line_num, &mut entries[0]);
//         // todo: reduce allowed difference to 0.0
//         assert_delta!(items[7], result.0, 1e-15, false, "result (a12 or result.0)", line_num);
//         if outmask & caps::DISTANCE != 0 {
//             entries[].add(items[8], result.1, line_num, &mut entries[1]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[8], result.1, 1e-15, false, "s12-out (result.1)", line_num);
//         }

//         entries[].add(items[9], result.2, line_num, &mut entries[2]);
//         // todo: reduce allowed difference to 0.0
//         // todo: disallow sign difference
//         assert_delta!(items[9], result.2, 1e-14, true, "salp1-out (result.2)", line_num);

//         entries[].add(items[10], result.3, line_num, &mut entries[3]);
//         // todo: reduce allowed difference to 0.0
//         assert_delta!(items[10], result.3, 1e-11, false, "calp1-out (result.3)", line_num);

//         entries[].add(items[11], result.4, line_num, &mut entries[4]);
//         // todo: reduce allowed difference to 0.0
//         // todo: disallow sign difference
//         assert_delta!(items[11], result.4, 1e-14, true, "salp2-out (result.4)", line_num);

//         entries[].add(items[12], result.5, line_num, &mut entries[5]);
//         // todo: reduce allowed difference to 0.0
//         assert_delta!(items[12], result.5, 1e-11, false, "calp2-out (result.5)", line_num);

//         if outmask & caps::REDUCEDLENGTH != 0 {
//             entries[].add(items[13], result.6, line_num, &mut entries[6]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[13], result.6, 1e-10, false, "m12-out (result.6)", line_num);
//         }
//         if outmask & caps::GEODESICSCALE != 0 {
//             entries[].add(items[14], result.7, line_num, &mut entries[7]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[14], result.7, 1e-14, false, "M12-out (result.7)", line_num);

//             entries[].add(items[15], result.8, line_num, &mut entries[8]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[15], result.8, 1e-14, false, "M21-out (result.8)", line_num);
//         }

//         if outmask & caps::AREA != 0 {
//             entries[].add(items[16], result.9, line_num, &mut entries[9]);
//             // todo: reduce allowed difference to 0.0
//             // todo: remove exceptions
//             if !result.9.is_nan() && !result.9.is_infinite() {
//                 assert_delta!(items[16], result.9, 1e-7, false, "S12-out (result.9)", line_num);
//             }
//         }
//     });

//     delta_entries.lock().unwrap().iter().enumerate().for_each(|(i, entry)| {
//         println!("test_gen_inverse_vs_cpp result.{} abs {:e} at line {}, rel {:e} at line {}, sign-change {}", i, entry.diff_abs, entry.line_abs, entry.diff_rel, entry.line_rel, entry.sign_change);
//     });
// }

// #[test]
// fn test_geodesic_vs_cpp() {
//     let consts = utilities::read_consts_basic("Geodesic_consts", 14);
//     // consts: nA1_ nC1_ nC1p_ nA2_ nC2_ nA3_ nA3x_ nC3_ nC3x_ nC4_ nC4x_ nC_ maxit1_ GEOGRAPHICLIB_GEODESIC_ORDER
    
//     // values from c++...
//     // nA3x_ = nA3_ = GEOGRAPHICLIB_GEODESIC_ORDER = 6
//     // nC3_ = GEOGRAPHICLIB_GEODESIC_ORDER = 6
//     // nC3x_ = (nC3_ * (nC3_ - 1)) / 2 = 15
//     // nC4_ = GEOGRAPHICLIB_GEODESIC_ORDER = 6
//     // nC4x = (nC4_ * (nC4_ + 1)) / 2 = 21
    
//     // 18 + nA3x_=6 + nC3x_=15 + nC4x_=21
//     // #[allow(non_snake_case)]
//     // let nA1_ = consts[0] as usize;
//     // #[allow(non_snake_case)]
//     // let nC1_ = consts[1] as usize;
//     // #[allow(non_snake_case)]
//     // let nA2_ = consts[2] as usize;
//     // #[allow(non_snake_case)]
//     // let nC2_ = consts[3] as usize;
//     // #[allow(non_snake_case)]
//     // let nA3_ = consts[4] as usize;
//     #[allow(non_snake_case)]
//     let nA3x_ = consts[6] as usize;
//     // #[allow(non_snake_case)]
//     // let nC3_ = consts[6] as usize;
//     #[allow(non_snake_case)]
//     let nC3x_ = consts[8] as usize;
//     // #[allow(non_snake_case)]
//     // let nC4_ = consts[8] as usize;
//     #[allow(non_snake_case)]
//     let nC4x_ = consts[10] as usize;
//     // #[allow(non_snake_case)]
//     // let nC_ = consts[10] as usize;

//     // Format: a f this-out[_a _f maxit2_ tiny_ tol0_ tol1_ tol2_ tolb_ xthresh_ _f1 _e2 _ep2 _n _b _c2 _etol2 _A3x(nA3x_) _C3x(nC3x_) _C4x(nC4x_)]
//     let arg_count = 18 + nA3x_ + nC3x_ + nC4x_;
//     let delta_entries = Arc::new(Mutex::new(vec![DeltaEntry::new(); 17]));
//     test_basic("Geodesic_Geodesic", -1, |line_num, items| {
//         assert!(items.len() == arg_count, "Expected {} items per line. Line {} had {}.", arg_count, line_num, items.len());
//         let g = Geodesic::new(items[0], items[1]);
//         // todo: try to enable commented out checks below (they check non-public values)
//         let mut entries = delta_entries.lock().unwrap();
//         entries[].add(items[2], g.a, line_num, &mut entries[2]);
//         assert_delta!(items[2], g.a, 0.0, false, "self.a", line_num);

//         entries[].add(items[3], g.f, line_num, &mut entries[3]);
//         assert_delta!(items[3], g.f, 0.0, false, "self.f", line_num);

//         // assert_delta!(items[4], g.maxit2_ as f64, 0.0, false, "self.maxit2_", line_num);

//         entries[].add(items[5], g.tiny_, line_num, &mut entries[5]);
//         assert_delta!(items[5], g.tiny_, 0.0, false, "self.tiny_", line_num);

//         // assert_delta!(items[6], g.tol0_, 0.0, false, "self.tol0_", line_num);
//         // assert_delta!(items[7], g.tol1_, 0.0, false, "self.tol1_", line_num);
//         // assert_delta!(items[8], g.tol2_, 0.0, false, "self.tol2_", line_num);
//         // assert_delta!(items[9], g.tolb_, 0.0, false, "self.tolb_", line_num);
//         // assert_delta!(items[10], g.xthresh_, 0.0, false, "self.xthresh_", line_num);

//         entries[].add(items[11], g._f1, line_num, &mut entries[11]);
//         assert_delta!(items[11], g._f1, 0.0, false, "self._f1", line_num);

//         entries[].add(items[12], g._e2, line_num, &mut entries[12]);
//         assert_delta!(items[12], g._e2, 0.0, false, "self._e2", line_num);

//         entries[].add(items[13], g._ep2, line_num, &mut entries[13]);
//         assert_delta!(items[13], g._ep2, 0.0, false, "self._ep2", line_num);

//         // assert_delta!(items[14], g._n, 0.0, false, "self._n", line_num);

//         entries[].add(items[15], g._b, line_num, &mut entries[15]);
//         assert_delta!(items[15], g._b, 0.0, false, "self._b", line_num);

//         entries[].add(items[16], g._c2, line_num, &mut entries[16]);
//         // todo: reduce accepted difference to 0.0
//         // todo: remove exception
//         if g.f > 0.0 {
//             assert_delta!(items[16], g._c2, 1e-15, false, "self._c2", line_num);
//         }
//         // assert_delta!(items[17], g._etol2, 0.0, false, "self._etol2", line_num);

//         // let mut i = 17;
//         // assert_eq!(nA3x_, g.nA3x_, "self.nA3x_ mismatch");
//         // assert_eq!(nA3x_, g.A3x_.len(), "self.A3x_ size mismatch");
//         // for item in g._A3x.iter() {
//         //     i += 1;
//         //     assert_delta!(items[i], *item, 0.0, false, "self._A3x item", line_num);
//         // }

//         // assert_eq!(nC3x_, g.nC3x_, "self.nC3x_ mismatch");
//         // assert_eq!(nC3x_, g.C3x_.len(), "self.C3x_ size mismatch");
//         // for item in g._C3x.iter() {
//         //     i += 1;
//         //     assert_delta!(items[i], *item, 0.0, false, "self._C3x item", line_num);
//         // }

//         // assert_eq!(nC4x_, g.nC4x_, "self.nC4x_ mismatch");
//         // assert_eq!(nC4x_, g.C4x_.len(), "self.C4x_ size mismatch");
//         // for item in g._C4x.iter() {
//         //     i += 1;
//         //     assert_delta!(items[i], *item, 0.0, false, "self._C4x item", line_num);
//         // }
//     });

//     delta_entries.lock().unwrap().iter().enumerate().for_each(|(i, entry)| {
//         if i > 1 {
//             println!("test_geodesic_vs_cpp result.{} abs {:e} at line {}, rel {:e} at line {}, sign-change {}", i, entry.diff_abs, entry.line_abs, entry.diff_rel, entry.line_rel, entry.sign_change);
//         }
//     });
// }

// // placeholder: Geodesic_InverseLine

// #[test]
// fn test_inverse_start_vs_cpp() {
//     // Format: this-in[_a _f] sbet1 cbet1 dn1 sbet2 cbet2 dn2 lam12 slam12 clam12 result=sig12 salp1-out calp1-out salp2-out calp2-out dnm-out
//     let delta_entries = Arc::new(Mutex::new(vec![DeltaEntry::new(); 6]));
//     test_basic("Geodesic_InverseStart", 17, |line_num, items| {
//         let g = Geodesic::new(items[0], items[1]);
//         #[allow(non_snake_case)]
//         let mut C1a: [f64; nC_] = [0.0 ; nC_];
//         #[allow(non_snake_case)]
//         let mut C2a: [f64; nC_] = [0.0 ; nC_];
//         let result = g._InverseStart(items[2], items[3], items[4], items[5], items[6], items[7], items[8], items[9], items[10], &mut C1a, &mut C2a);
//         let mut entries = delta_entries.lock().unwrap();

//         entries[].add(items[11], result.0, line_num, &mut entries[0]);
//         assert_delta!(items[11], result.0, 0.0, false, "result (sig12 or result.0)", line_num);

//         entries[].add(items[12], result.1, line_num, &mut entries[1]);
//         // todo: reduce allowed difference to 0.0
//         // todo: pass allow_sign_change=false
//         assert_delta!(items[12], result.1, 1e-15, true, "salp1-out (result.1)", line_num);

//         entries[].add(items[13], result.2, line_num, &mut entries[2]);
//         // todo: reduce allowed difference to 0.0
//         assert_delta!(items[13], result.2, 1e-15, false, "calp1-out (result.2)", line_num);

//         entries[].add(items[14], result.3, line_num, &mut entries[3]);
//         // todo: remove exceptions
//         if !f64::is_nan(result.3) {
//             assert_delta!(items[14], result.3, 0.0, false, "salp2-out (result.3)", line_num);
//         }

//         entries[].add(items[15], result.4, line_num, &mut entries[4]);
//         // todo: remove exceptions
//         if !f64::is_nan(result.4) {
//             assert_delta!(items[15], result.4, 0.0, false, "calp2-out (result.4)", line_num);
//         }

//         entries[].add(items[16], result.5, line_num, &mut entries[5]);
//         // todo: remove exceptions
//         if !f64::is_nan(result.5) {
//             assert_delta!(items[16], result.5, 0.0, false, "dnm-out (result.5)", line_num);
//         }
//     });

//     delta_entries.lock().unwrap().iter().enumerate().for_each(|(i, entry)| {
//         println!("test_inverse_start_vs_cpp result.{} abs {:e} at line {}, rel {:e} at line {}, sign-change {}", i, entry.diff_abs, entry.line_abs, entry.diff_rel, entry.line_rel, entry.sign_change);
//     });
// }

// #[test]
// fn test_lambda12_vs_cpp() {
//     // Format: this-in[_a _f] sbet1 cbet1 dn1 sbet2 cbet2 dn2 salp1 calp1 slam120 clam120 diffp result=lam12 salp2-out calp2-out sig12-out ssig1-out csig1-out ssig2-out csig2-out eps-out domg12-out dlam12-out
//     let delta_entries = Arc::new(Mutex::new(vec![DeltaEntry::new(); 11]));
//     test_basic("Geodesic_Lambda12", 24, |line_num, items| {
//         let g = Geodesic::new(items[0], items[1]);
//         #[allow(non_snake_case)]
//         let mut C1a: [f64; nC_] = [0.0 ; nC_];
//         #[allow(non_snake_case)]
//         let mut C2a: [f64; nC_] = [0.0 ; nC_];
//         #[allow(non_snake_case)]
//         let mut C3a: [f64; nC_] = [0.0 ; nC_];
//         // todo: review rs approach of modifying calp1. c++ counterpart does not, so rs should probably at least comment on the difference
//         let mut calp1 = items[9];
//         let result = g._Lambda12(items[2], items[3], items[4], items[5], items[6], items[7], items[8], &mut calp1, items[10], items[11], items[12] != 0.0, &mut C1a, &mut C2a, &mut C3a);
//         let mut entries = delta_entries.lock().unwrap();

//         entries[].add(items[13], result.0, line_num, &mut entries[0]);
//         // todo: reduce allowed difference to 0.0
//         // todo: pass allow_sign_change=false
//         assert_delta!(items[13], result.0, 4e0, true, "result (lam12 or result.0)", line_num);

//         entries[].add(items[14], result.1, line_num, &mut entries[1]);
//         assert_delta!(items[14], result.1, 0.0, false, "salp2-out (result.1)", line_num);

//         entries[].add(items[15], result.2, line_num, &mut entries[2]);
//         assert_delta!(items[15], result.2, 0.0, false, "calp2-out (result.2)", line_num);

//         entries[].add(items[16], result.3, line_num, &mut entries[3]);
//         assert_delta!(items[16], result.3, 0.0, false, "sig12-out (result.3)", line_num);

//         entries[].add(items[17], result.4, line_num, &mut entries[4]);
//         assert_delta!(items[17], result.4, 0.0, false, "ssig1-out (result.4)", line_num);

//         entries[].add(items[18], result.5, line_num, &mut entries[5]);
//         assert_delta!(items[18], result.5, 0.0, false, "csig1-out (result.5)", line_num);

//         entries[].add(items[19], result.6, line_num, &mut entries[6]);
//         assert_delta!(items[19], result.6, 0.0, false, "ssig2-out (result.6)", line_num);

//         entries[].add(items[20], result.7, line_num, &mut entries[7]);
//         assert_delta!(items[20], result.7, 0.0, false, "csig2-out (result.7)", line_num);

//         entries[].add(items[21], result.8, line_num, &mut entries[8]);
//         assert_delta!(items[21], result.8, 0.0, false, "eps-out (result.8)", line_num);

//         entries[].add(items[22], result.9, line_num, &mut entries[9]);
//         assert_delta!(items[22], result.9, 0.0, false, "domg12-out (result.9)", line_num);

//         entries[].add(items[23], result.10, line_num, &mut entries[10]);
//         // todo: reduce allowed difference to 0.0
//         assert_delta!(items[23], result.10, 1e-5, false, "dlam12-out (result.10)", line_num);
//     });

//     delta_entries.lock().unwrap().iter().enumerate().for_each(|(i, entry)| {
//         println!("test_lambda12_vs_cpp result.{} abs {:e} at line {}, rel {:e} at line {}, sign-change {}", i, entry.diff_abs, entry.line_abs, entry.diff_rel, entry.line_rel, entry.sign_change);
//     });
// }

// #[test]
// fn test_lengths_vs_cpp() {
//     // Format: this-in[_a _f] eps sig12 ssig1 csig1 dn1 ssig2 csig2 dn2 cbet1 cbet2 outmask s12b-out m12b-out m0-out M12-out M21-out
//     let delta_entries = Arc::new(Mutex::new(vec![DeltaEntry::new(); 5]));
//     test_basic("Geodesic_Lengths", 18, |line_num, items| {
//         let g = Geodesic::new(items[0], items[1]);
//         #[allow(non_snake_case)]
//         let mut C1a: [f64; nC_] = [0.0 ; nC_];
//         #[allow(non_snake_case)]
//         let mut C2a: [f64; nC_] = [0.0 ; nC_];
//         let outmask = items[12] as u64;
//         let result = g._Lengths(items[2], items[3], items[4], items[5], items[6], items[7], items[8], items[9], items[10], items[11], outmask, &mut C1a, &mut C2a);
//         let mut entries = delta_entries.lock().unwrap();

//         if outmask & caps::DISTANCE != 0 {
//             entries[].add(items[13], result.0, line_num, &mut entries[0]);
//             // todo: remove exceptions
//             if !f64::is_nan(result.0) {
//                 assert_delta!(items[13], result.0, 0.0, false, "s12b-out (result.0)", line_num);
//             }
//         }

//         if outmask & caps::REDUCEDLENGTH != 0 {
//             entries[].add(items[14], result.1, line_num, &mut entries[1]);
//             // todo: reduce allowed difference to 0.0
//             assert_delta!(items[14], result.1, 1e-5, false, "m12b-out (result.1)", line_num);

//             entries[].add(items[15], result.2, line_num, &mut entries[2]);
//             assert_delta!(items[15], result.2, 0.0, false, "m0-out (result.2)", line_num);
//         }

//         if outmask & caps::GEODESICSCALE != 0 {
//             entries[].add(items[16], result.3, line_num, &mut entries[3]);
//             // todo: remove exceptions
//             if !f64::is_nan(result.3) {
//                 assert_delta!(items[16], result.3, 0.0, false, "M12-out (result.3)", line_num);
//             }

//             entries[].add(items[17], result.4, line_num, &mut entries[4]);
//             // todo: remove exceptions
//             if !f64::is_nan(result.4) {
//                 assert_delta!(items[17], result.4, 0.0, false, "M21-out (result.4)", line_num);
//             }
//         }
//     });

//     delta_entries.lock().unwrap().iter().enumerate().for_each(|(i, entry)| {
//         println!("test_lengths_vs_cpp result.{} abs {:e} at line {}, rel {:e} at line {}, sign-change {}", i, entry.diff_abs, entry.line_abs, entry.diff_rel, entry.line_rel, entry.sign_change);
//     });
// }

// // Note: For Geodesic_SinCosSeries see geomath_tests.rs
