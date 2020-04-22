use std::env;
use std::fs::File;
use std::io::{self, BufRead, BufReader};

use geographiclib_rs::{capability, Geodesic};
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
        let fields: Vec<f64> = line.split(" ").map(|s| s.parse::<f64>().unwrap()).collect();
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
        let (
            computed_lat2,
            computed_lon2,
            computed_azi2,
            computed_m12,
            computed_M12,
            computed_M21,
            computed_S12,
            computed_a12,
        ) = self.geod.Direct(lat1, lon1, azi1, s12);
        let output_fields = if self.is_full_output {
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
        };
        output_fields
    }

    fn compute_inverse(&self, fields: &Vec<f64>) -> Vec<f64> {
        assert_eq!(4, fields.len());
        let lat1 = fields[0];
        let lon1 = fields[1];
        let lat2 = fields[2];
        let lon2 = fields[3];

        let outmask = capability::ALL;
        let result = self.geod.Inverse(lat1, lon1, lat2, lon2, outmask);

        let output_fields = if self.is_full_output {
            // TODO - we're currently omitting several fields, and only outputting what's
            //        necessary to pass the validation tool
            vec![
                lat1,
                lon1,
                result["azi1"],
                lat2,
                lon2,
                result["azi2"],
                result["s12"],
                result["a12"],
                result["m12"],
            ]
        } else {
            vec![result["azi1"], result["azi2"], result["s12"]]
        };
        output_fields
    }
}
