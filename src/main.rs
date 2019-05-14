mod geodesic;
mod geodesiccapability;
mod geomath;

fn main() {
    println!("{}", geomath::sq(12.0));
    println!("{}", geomath::cbrt(2.0));
    println!("{:?}", geomath::norm(2.0, 3.0));
    println!("{:?}", geomath::sum(2.0, 3.0));

    let coeff = vec![
        -3.0, 128.0, -2.0, -3.0, 64.0, -1.0, -3.0, -1.0, 16.0, 3.0, -1.0, -2.0, 8.0, 1.0, -1.0,
        2.0, 1.0, 1.0,
    ];
    println!(
        "{:?}",
        geomath::polyval(1, &coeff, 2, 0.0016792203863837047 / coeff[4] as f64)
    );
    println!("{:?}", geomath::ang_round(0.1));
    println!("{:?}", geomath::ang_normalize(0.1));
    println!("{:?}", geomath::lat_fix(-91.0));
    println!("{:?}", geomath::lat_fix(-80.0));
    println!("{:?}", geomath::ang_diff(-80.0, 43.0));
    println!("{:?}", geomath::sincosd(-80.0));
    println!("{:?}", geomath::atan2d(-80.0, 20.0));
    println!("{:?}", geomath::isfinite(20.0));
    println!("{:?}", geomath::isfinite(200000000000000.0));
    println!("{}", geodesiccapability::ALL);

    const WGS84_A: f64 = 6378137.0;
    const WGS84_F: f64 = 1.0 / 298.257223563;
    let geod = geodesic::Geodesic::new(WGS84_A, WGS84_F);
    println!("{:?}", geod);
}
