use geographiclib_rs::{DirectGeodesic, Geodesic};
fn main() {
    let g = Geodesic::wgs84();
    let jfk_lat = 40.64;
    let jfk_lon = -73.78;
    let northeast_azimuth = 45.0;

    let (lat, lon, az) = g.direct(jfk_lat, jfk_lon, northeast_azimuth, 10e6);

    use approx::assert_relative_eq;
    assert_relative_eq!(lat, 32.621100463725796);
    assert_relative_eq!(lon, 49.05248709295982, epsilon = 1e-13);
    assert_relative_eq!(az, 140.4059858768007);

    println!("lat: {lat}, lon: {lon}, az: {az}");
}
