use geographiclib_rs::{DirectGeodesic, Geodesic};
use std::fs::File;
use std::io::{BufRead, BufReader};

fn main() {
    let geod = Geodesic::wgs84();
    let file = File::open("input/GeodTest.dat").unwrap();
    let reader = BufReader::new(file);

    let inputs = reader.lines().map(|line| {
        let line = line.unwrap();
        let fields: Vec<f64> = line.split(" ").map(|s| s.parse::<f64>().unwrap()).collect();
        (fields[0], fields[1], fields[2], fields[6])
    });

    for (lat1, lon1, azi1, s12) in inputs {
        #[allow(non_snake_case)]
        let (_lat2, _lon2, _azi2, _m12, _M12, _M21, _S12, _a12) =
            geod.direct(lat1, lon1, azi1, s12);
    }
}
