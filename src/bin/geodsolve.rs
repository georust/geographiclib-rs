use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use geographiclib_rs::{DirectGeodesic, Geodesic, InverseGeodesic};
use std::error::Error;

struct Runner {
    geod: Geodesic,
    is_full_output: bool,
    is_inverse: bool,
    input_filename: Option<String>,
}

fn main() -> Result<(), Box<dyn Error>> {
    let geod = Geodesic::wgs84();
    let mut is_full_output_arg: Option<bool> = None;
    let mut is_inverse_arg: Option<bool> = None;
    let mut input_file_arg: Option<String> = None;
    let mut next_arg_is_inputfile = false;

    for argument in env::args() {
        if argument == "--input-file" {
            next_arg_is_inputfile = true;
        } else if next_arg_is_inputfile {
            next_arg_is_inputfile = false;
            input_file_arg = Some(argument);
        } else if argument == "-i" {
            is_inverse_arg = Some(true);
        } else if argument == "-f" {
            is_full_output_arg = Some(true);
        }
    }

    let is_full_output = is_full_output_arg.unwrap_or(false);
    let is_inverse = is_inverse_arg.unwrap_or(false);
    Runner::new(geod, is_full_output, is_inverse, input_file_arg).run()
}

impl Runner {
    pub fn new(
        geod: Geodesic,
        is_full_output: bool,
        is_inverse: bool,
        input_filename: Option<String>,
    ) -> Self {
        Runner {
            geod,
            is_full_output,
            is_inverse,
            input_filename,
        }
    }

    pub fn run(&self) -> Result<(), Box<dyn Error>> {
        if let Some(input_filename) = &self.input_filename {
            let file = File::open(input_filename)?;
            let reader = BufReader::new(file);
            for line in reader.lines() {
                let line = line.unwrap();
                self.handle_line(line);
            }
        } else {
            for line in io::stdin().lock().lines() {
                let line = line.unwrap();
                self.handle_line(line);
            }
        }
        Ok(())
    }

    fn handle_line(&self, line: String) {
        let fields: Vec<f64> = line.split(' ').map(|s| s.parse::<f64>().unwrap()).collect();
        let output_fields = if self.is_inverse {
            self.compute_inverse(&fields)
        } else {
            self.compute_direct(&fields)
        };
        let output_strings: Vec<String> = output_fields.iter().map(|f| f.to_string()).collect();
        let output_line = output_strings.join(" ");
        println!("{}", output_line);
    }

    fn compute_direct(&self, fields: &Vec<f64>) -> Vec<f64> {
        assert_eq!(4, fields.len());
        let lat1 = fields[0];
        let lon1 = fields[1];
        let azi1 = fields[2];
        let s12 = fields[3];
        #[allow(non_snake_case)]
        let (
            computed_lat2,
            computed_lon2,
            computed_azi2,
            computed_m12,
            _computed_M12,
            _computed_M21,
            _computed_S12,
            computed_a12,
        ) = self.geod.direct(lat1, lon1, azi1, s12);

        if self.is_full_output {
            // TODO - we're currently omitting several fields, and only outputting what's
            //        necessary to pass the validation tool
            vec![
                lat1,
                lon1,
                azi1,
                computed_lat2,
                computed_lon2,
                computed_azi2,
                s12,
                computed_a12,
                computed_m12,
            ]
        } else {
            vec![computed_lat2, computed_lon2, computed_azi2]
        }
    }

    fn compute_inverse(&self, fields: &Vec<f64>) -> Vec<f64> {
        assert_eq!(4, fields.len());
        let input_lat1 = fields[0];
        let input_lon1 = fields[1];
        let input_lat2 = fields[2];
        let input_lon2 = fields[3];

        #[allow(non_snake_case)]
        let (s12, azi1, azi2, m12, _M12, _M21, _S12, a12) = self
            .geod
            .inverse(input_lat1, input_lon1, input_lat2, input_lon2);

        if self.is_full_output {
            // TODO - we're currently omitting several fields, and only outputting what's
            //        necessary to pass the validation tool
            vec![
                input_lat1, input_lon1, azi1, input_lat2, input_lon2, azi2, s12, a12, m12,
            ]
        } else {
            vec![azi1, azi2, s12]
        }
    }
}
