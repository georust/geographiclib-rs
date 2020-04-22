use geographiclib_rs::{Geodesic, InverseGeodesic};
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() {
    let geod = Geodesic::wgs84();
    let file = File::open("input/GeodTest.dat").unwrap();
    let reader = BufReader::new(file);

    let inputs = reader.lines().map(|line| {
        let line = line.unwrap();
        let fields: Vec<f64> = line.split(" ").map(|s| s.parse::<f64>().unwrap()).collect();
        (fields[0], fields[1], fields[3], fields[4])
    });

    for (lat1, lon1, lat2, lon2) in inputs {
        #[allow(non_snake_case)]
        let (_s12, _azi1, _azi2, _m12, _M12, _M21, _S12, _a12) =
            geod.inverse(lat1, lon1, lat2, lon2);
    }
}
