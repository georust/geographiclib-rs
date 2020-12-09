// #[macro_use]

use std::f64;
use std::fs::File;
use std::io::{self, BufRead};
use std::iter::Iterator;
use std::str;
use std::vec::Vec;


// todo: remove the *.dat files from residing directly in this repo, and instead pull them externally?
// todo: switch from truncated *.dat files to full *.dat files (after getting them out of this repo)

pub const NON_NUMS_NONE: &[usize; 0] = &[];


// Return the absolute and relative difference between two values, in that order.
// If both values are nan, consider the difference to be 0.
pub fn calc_delta(x: f64, y: f64) -> (f64, f64) {
  if x.is_nan() && y.is_nan() {
    return (0f64, 0f64);
  } else {
    let diff_abs = (x - y).abs();
    let diff_rel = 2.0 * diff_abs / (x.abs() + y.abs());
    return (diff_abs, diff_rel);
  }
}

// Check whether x and y are within delta of each other, considering both absolute and relative difference.
// Consider them close enough if absolute difference is within delta inclusive,
// or if relative difference is within delta exclusive.
// nan values are considered equal to each other, and outside delta from everything else.
#[macro_export]
macro_rules! assert_delta {
  ($x:expr, $y:expr, $delta:expr, $allow_sign_change:expr, $msg:expr) => {
      let delta2 = utilities::calc_delta($x, $y);
      // For the sign change check, allow (NAN vs NAN), but not (0.0 vs -0.0) or (NAN vs -NAN).
      assert!(
        $allow_sign_change || $x.is_sign_negative() == $y.is_sign_negative(),
        "assert_delta failed {}: {:e} vs {:e} sign difference disallowed", $msg, $x, $y
      );
      assert!(
          delta2.0 <= $delta || delta2.1 < $delta,
          "assert_delta failed {}: {:e} vs {:e} diff {:e} rel {:e} outside inclusive {:e}", $msg, $x, $y, delta2.0, delta2.1, $delta
      );
  }
}

// The parsed contents of a data file content line, not including operation type lines.
pub struct DataLine {
  // The full (except for trimming) original line.
  pub full_line: String,

  // Input values parsed as numbers. Uses nan for non-numeric values, and an empty vector if the section is missing or empty.
  pub inputs_num: Vec<f64>,

  // Output values parsed as numbers. Uses nan for non-numeric values, and an empty vector if the section is missing or empty.
  pub outputs_num: Vec<f64>,

  // Input values as strings. Uses an empty vector if the section is missing or empty.
  // Note that the values aren't truly raw, in that some escape sequences are translated.
  pub inputs_raw: Vec<String>,

  // Output values as strings. Uses an empty vector if the section is missing or empty.
  // Note that the values aren't truly raw, in that some escape sequences are translated.
  pub outputs_raw: Vec<String>,
  
  // Time in seconds reported for the line. If not present, nan is reported.
  pub time_sec: f64,

  // True if the line indicates that an error occurred for the set of inputs this line represents, when the line was created.
  // Typically, this means that the line represents inputs that the operation should reject.
  // Please note that this does NOT indicate that there was an error reading or parsing the line.
  pub error_flag: bool,
}

pub struct DataRequirements<'a> {
  pub op_type: &'a str,
  pub input_sizes: &'a [usize],
  pub output_sizes: &'a [usize],
  pub non_nums_in: &'a [usize],
  pub non_nums_out: &'a [usize],
}

impl DataRequirements<'_> {
  // Get an iterator over non-blank, non-comment lines of a data file.
  // Returns an Iterator over the trimmed lines of the file.
  pub fn read_data_lines<'a>(&'a self, file: &'a File) -> impl Iterator<Item=DataLine> + 'a {
    // Pass in the reader, rather than creating it here,
    // so that it can control the lifetime of the iterator.
    let reader = io::BufReader::new(file);
    let lines = reader.lines();
    let mapped = lines.map(|line| {
      line.expect("Failed to read data file line").trim().to_owned()
    });
    let filtered = mapped.filter(
      move |line| line.len() > 0 && !line.starts_with("#") && line.ne(self.op_type)
    );
    let parsed = filtered.map(
      move |line| self.parse(&line)
    );
    return parsed;
  }

  // Parse a data file content line (but not an operation line).
  // Expect to find non-numeric inputs at the 0-based indices indicated by non_nums_in (in ascending order).
  // Expect to find non-numeric outputs at the 0-based indices indicated by non_nums_out (in ascending order).
  // All other non-numeric inputs and outputs will be interpreted as errors.
  // Update parsed_line to reflect that values found.
  pub fn parse(&self, line: &str) -> DataLine {
    let mut parsed_line = DataLine {
      error_flag: false,
      full_line: line.to_string(),
      inputs_num: Vec::new(),
      outputs_num: Vec::new(),
      inputs_raw: Vec::new(),
      outputs_raw: Vec::new(),
      time_sec: f64::NAN,
    };
  
    let section_strs: Vec<&str> = line.split(';').map(|section_str| section_str.trim()).collect();
    let section_count = section_strs.len();
    if section_count > 0 {
      // First section: inputs
      self.read_values(section_strs[0], true, &mut parsed_line);
    }
    if section_count > 1 {
      // Second section: error flag
      let section = section_strs[1];
      if section.len() < 1 {
        // Do nothing
      } else if section.starts_with("ERROR") {
        // It's an error flag. Don't worry about whether ERROR is followed by anything else
        parsed_line.error_flag = true;
      } else {
        panic!("Invalid error section content. Expected empty or value starting with ERROR.");
      }
    }
    // If there's an error flag, don't even worry about recording remaining values, because their meaning is limited.
    if section_count > 2 && !parsed_line.error_flag {
      // Third section: outputs
      self.read_values(section_strs[2], false, &mut parsed_line);
    }
    if section_count > 3 && !parsed_line.error_flag {
      // Fourth section: time
      parsed_line.time_sec = section_strs[3].parse::<f64>().unwrap();
    }

    parsed_line
  }

  // Split section_str on spaces, and add the results to inputs_raw and (after parsing) inputs_num.
  // Indices included in non_nums are treated as nan without trying to parse.
  // non_nums is assumed to be in increasing order.
  pub fn read_values(&self, section_str: &str, read_inputs: bool, parsed_line: &mut DataLine) {
    // let items_raw = if read_inputs { &mut parsed_line.inputs_raw } else { &mut parsed_line.outputs_raw };
    // let items_num = if read_inputs { &mut parsed_line.inputs_num } else { &mut parsed_line.outputs_num };
    // let count_expected = if read_inputs { &self.input_sizes } else { &self.output_sizes };
    // let non_nums = if read_inputs {&self.non_nums_in} else {&self.non_nums_out};
    let (items_raw, items_num, count_expected, non_nums) = if read_inputs {
      (&mut parsed_line.inputs_raw, &mut parsed_line.inputs_num, &self.input_sizes, &self.non_nums_in)
    } else {
      (&mut parsed_line.outputs_raw, &mut parsed_line.outputs_num, &self.output_sizes, &self.non_nums_out) 
    };

    section_str.split(' ').filter(|item| item.len() > 0).for_each(|item| items_raw.push(item.to_owned()));
    let item_count = items_raw.len();
    if count_expected.len() > 0 && !count_expected.contains(&item_count) {
      // let count_combined: [String] = (**count_expected).iter().map(|count| count.to_string()).join("/");
      let mut count_combined = count_expected[0].to_string();
      for i in 1..count_expected.len() {
        count_combined.push('/');
        count_combined.push_str(&count_expected[i].to_string());
      }
      let item_type = if read_inputs { "in" } else { "out" };
      panic!("Expected {} {}puts but found {}", count_combined, item_type, item_count);
    }
    let mut non_num_meta_index = 0;
    for i in 0..item_count {
      let mod_item = unescape_all(&items_raw[i]);
      (*items_raw)[i] = mod_item;
      if non_num_meta_index < non_nums.len() && non_nums[non_num_meta_index] == i {
        // The value isn't expected to be numeric. Don't even try to parse it.
        items_num.push(f64::NAN);
        non_num_meta_index += 1;
      } else {
        let item = &*items_raw[i]; // convert from String to &str
        let as_num = match item {
          "nan" => f64::NAN,
          "-nan" => -f64::NAN,
          "inf" => f64::INFINITY,
          "-inf" => f64::NEG_INFINITY,
          _ => item.parse::<f64>().unwrap(),
        };
        items_num.push(as_num);
      }
    }
  }

}

fn unescape_all(s: &String) -> String {
  let chars1 = s.as_bytes();
  let mut chars2 = chars1.to_owned();
  let len = chars1.len();
  let mut i_from: usize = 0;
  let mut i_to: usize = 0;
  const BACKSLASH: u8 = '\\' as u8;
  const LOWER_S: u8 = 's' as u8;
  const SPACE: u8 = ' ' as u8;
  const UPPER_S: u8 = 'S' as u8;
  const SEMI: u8 = ';' as u8;
  const ZERO: u8 = '0' as u8;
  const NOTHING: u8 = '\0' as u8;
  while i_from < len {
    let c = chars1[i_from];
    i_from += 1;
    if c == BACKSLASH {
      if i_from >= len {
        panic!("Invalid escape sequence. Item ended with single backslash.");
      }
      let c_swap = match chars1[i_from] {
        BACKSLASH => BACKSLASH,
        LOWER_S => SPACE,
        UPPER_S => SEMI,
        ZERO => NOTHING,
        _ => panic!("Invalid escape sequence \\{}", &s[i_from..i_from+1])
      };
      if c_swap != NOTHING {
        chars2[i_to] = c_swap;
        i_to += 1;
      }
    } else {
      chars2[i_to] = c;
      i_to += 1;
    }
  }

  for i in i_to..len {
    chars2[i] = NOTHING;
  }
  let copy = String::from_utf8(chars2).expect("Failed to create new string");
  return copy;
}


// #[test]
// fn test_unescape_all() {
//   // Math_AngDiff_x_y_e
//   let before = "Math_AngDiff_x_y_e";
//   let after_expected = before;
//   let after_found = unescape_all(&before.to_string());
//   assert_eq!(after_expected, after_found);
// }