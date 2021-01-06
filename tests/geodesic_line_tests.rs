// todo: restore and complete if GeodesicLine is made public or if these become unit tests

// // Integration tests related to geodesic_line.

// extern crate utilities;

// use utilities::RecordDeltaEntry;
// use std::sync::Mutex;
// use std::sync::Arc;
// use utilities::record_delta;
// use geographiclib_rs::{GeodesicLine, geodesiccapability as caps};
// use utilities::{assert_delta, test_basic};

// // reminder: use "cargo test -- --nocapture" to let tests print to console.
// // reminder: if a test fails especially badly, review the C++ "instrumented-crude" logging for possible bugs there

// #[test]
// fn test_consts_vs_cpp() {
//     // Format: nC1_ nC1p_ nC2_ nC3_ nC4_
//     let items = utilities::read_consts_basic("GeodesicLine_consts", 5);
//     let line_num = 2; // todo: remove once assert_delta stops requiring a line number
//     assert_delta!(items[0], GeodesicLine::nC1_, 0.0, false, "nC1_", line_num);
//     assert_delta!(items[1], GeodesicLine::nC1p_, 0.0, false, "nC1p_", line_num);
//     assert_delta!(items[2], GeodesicLine::nC2_, 0.0, false, "nC2_", line_num);
//     assert_delta!(items[3], GeodesicLine::nC3_, 0.0, false, "nC3_", line_num);
//     assert_delta!(items[4], GeodesicLine::nC4_, 0.0, false, "nC4_", line_num);
// }


// // placeholder: GeodesicLine_GenPosition: this-in[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0] arcmode s12_a12 outmask result=a12 lat2-out lon2-out azi2-out s12-out m12-out M12-out M21-out S12-out
// // placeholder: GeodesicLine_GenSetDistance: this-in[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)] arcmode s13_a13 this-out[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0
// // placeholder: GeodesicLine_GeodesicLine_5arg: g[_a _f] lat1 lon1 azi1 caps this-out[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)]
// // placeholder: GeodesicLine_GeodesicLine_9arg: g[_a _f] lat1 lon1 azi1 salp1 calp1 caps arcmode s13_a13 this-out[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)]
// // placeholder: GeodesicLine_LineInit: this-in[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)] g[_a _f] lat1 lon1 azi1 salp1 calp1 caps this-out[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)]
// // placeholder: GeodesicLine_SetArc: this-in[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)] a13 this-out[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0]
// // placeholder: GeodesicLine_SetDistance: this-in[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0 tiny_ _b _c2 _f1 _k2 _salp1 _calp1 _ssig1 _csig1 _dn1 _stau1 _ctau1 _somg1 _comg1 _A1m1 _A2m1 _A3c _B11 _B21 _B31 _A4 _B41 _C1a(nC1_+1) _C1pa(nC1p_+1) _C2a(nC2_+1) _C3a(nC3_) _C4a(nC4_)] s13 this-out[_a _f _lat1 _lon1 _azi1 _a13 _s13 _caps _salp0 _calp0]
