extern crate geographiclib_rs;

use crate::util;
use std::collections::{BTreeMap, HashMap};
use std::fmt::Display;


// A struct for taking difference values, splitting into buckets, and reporting back.
pub struct DeltaHistogram {
    pub num_nan: usize,
    pub num_inf: usize,
    pub num_zero: usize,
    pub exp_buckets: HashMap<isize, usize>,
}

impl DeltaHistogram {
    pub fn new() -> Self {
        DeltaHistogram {
            num_nan: 0,
            num_inf: 0,
            num_zero: 0,
            exp_buckets: HashMap::new(),
        }
    }

    // Add a new difference item to the dataset being tracked.
    pub fn add(&mut self, delta: f64) {
        if delta.is_nan() {
            self.num_nan += 1;
        } else if delta.is_infinite() {
            self.num_inf += 1;
        } else if delta == 0.0 {
            self.num_zero += 1;
        } else {
            let exp = delta.log10() as isize;
            let current: usize = match self.exp_buckets.get(&exp) {
                Some(val) => *val,
                _ => 0,
            };
            self.exp_buckets.insert(exp, current + 1);
        }
    }
}

impl Clone for DeltaHistogram {
    fn clone(&self) -> Self {
        DeltaHistogram {
            num_nan: self.num_nan,
            num_inf: self.num_inf,
            num_zero: self.num_zero,
            exp_buckets: self.exp_buckets.clone(),
        }
    }
}

// Round a value for use in DeltaHistogram display.
// Never round to 0 or 100. Only accept those values naturally.
fn to_percent(num_part: usize, num_all: usize) -> usize {
    let percent = 100f64 * num_part as f64 / num_all as f64;
    let rounded = if percent < 1.0 && num_part != 0 {
        1
    } else if percent > 99.0 && num_part != num_all {
        99
    } else {
        percent.round() as usize
    };
    rounded
}

impl Display for DeltaHistogram {
    // Display a summary, reduced down to a manageable number of buckets.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        const MAX_BUCKETS: usize = 7;
        let mut keys_asc: Vec<isize> = Vec::new();
        // histo_reduced map's keys are the original exponent.
        // Its values are (reduced_exponent_min, reduced_exponent_max, count).
        let mut histo_reduced: BTreeMap<isize, (isize, isize, usize)> = BTreeMap::new();
        let mut num_total = self.num_inf + self.num_nan + self.num_zero;
        self.exp_buckets.iter().for_each(|(&key, &val)| {
            keys_asc.push(key);
            histo_reduced.insert(key, (key, key, val));
            num_total += val;
        });
        keys_asc.sort();
        while histo_reduced.len() > MAX_BUCKETS {
            // Collapse the smallest bucket into its least-populated neighbor.
            let mut collapse_from = isize::MIN;
            let mut val_smallest = (collapse_from, collapse_from, usize::MAX);
            histo_reduced.iter().for_each(|(&key, &(_exp_min, _exp_max, count))| {
                if count < val_smallest.2 {
                    collapse_from = key;
                    val_smallest = (key, key, count);
                }
            });

            let index_smallest = keys_asc.iter().position(|&val| val == collapse_from).unwrap();
            // Note that we don't have to worry about the case of 2 or fewer buckets,
            // because we stop looping before that happens.
            let (collapse_to, val_to) = if index_smallest == 0 {
                let key_next = keys_asc[index_smallest + 1];
                let val_next = histo_reduced.get(&key_next).unwrap();
                (key_next, val_next)
            } else if index_smallest >= histo_reduced.len() - 1 {
                let key_prev = keys_asc[index_smallest - 1];
                let val_prev = histo_reduced.get(&key_prev).unwrap();
                (key_prev, val_prev)
            } else {
                // Favor collapsing into the smaller bucket, to reduce lopsided bucket sizes
                let key_prev = keys_asc[index_smallest - 1];
                let key_next = keys_asc[index_smallest + 1];
                let val_prev = histo_reduced.get(&key_prev).unwrap();
                let val_next = histo_reduced.get(&key_next).unwrap();
                if val_next.2 < val_prev.2 {
                    (key_next, val_next)
                } else {
                    (key_prev, val_prev)
                }
            };

            let val_sum = (isize::min(val_to.0, val_smallest.0), isize::max(val_to.1, val_smallest.1), val_to.2 + val_smallest.2);
            histo_reduced.insert(collapse_to, val_sum);
            histo_reduced.remove(&collapse_from);
            keys_asc.remove(index_smallest);
            assert_eq!(keys_asc.len(), histo_reduced.len(), "Size mismatch between key list and map");
        }

        write!(f, "count {}", num_total)?;

        if self.num_zero > 0 {
            let percent_zero = to_percent(self.num_zero, num_total); 
            write!(f, ", zero {}%", percent_zero)?;
        }

        // Convert counts to percentages
        histo_reduced.iter_mut().for_each(|(_key, (_exp_min, _exp_max, count))| {
            assert!(*count != 0, "Internal error: Bucket contains no items");
            *count = to_percent(*count, num_total);
        });
        for (key, (exp_min, exp_max, count)) in &histo_reduced {
            if exp_min == exp_max {
                write!(f, ", e{} {}%", key, count)?;
            } else {
                write!(f, ", e{} to e{} {}%", exp_min, exp_max, count)?;
            }
        }
        if self.num_inf > 0 {
            let percent_inf = to_percent(self.num_inf, num_total);
            write!(f, ", inf {}%", percent_inf)?;
        }
        if self.num_nan > 0 {
            let percent_nan = to_percent(self.num_nan, num_total);
            write!(f, ", nan {}%", percent_nan)?;
        }
        Ok(())
    }
}

// Summary of count of times a condition occurs for DeltaEntry,
// and information about a sample line where it occurred (which
// may be the first or worst, depending on context).
pub struct DeltaLineSummary {
    pub x: f64,
    pub y: f64,
    pub line_num: usize,
    pub count: usize,
}

impl Copy for DeltaLineSummary {
}

impl Clone for DeltaLineSummary {
    fn clone(&self) -> Self {
        DeltaLineSummary {
            x: self.x,
            y: self.y,
            line_num: self.line_num,
            count: self.count,
        }
    }
}

impl DeltaLineSummary {
    fn new() -> Self {
        DeltaLineSummary {
            x: f64::NAN,
            y: f64::NAN,
            line_num: 0,
            count: 0,
        }
    }

    // Update the summary based on a line.
    // If "worst" is true, update x, y, and line_num even if this isn't the first line added.
    fn add(&mut self, x: f64, y: f64, line_num: usize, worst: bool) {
        if worst || self.count == 0 {
            self.x = x;
            self.y = y;
            self.line_num = line_num;
        }
        self.count += 1;
    }
}

// Information about value differences.
// This is typically used to record the worst-case differences found
// among sets of comparable expected vs found values.
pub struct DeltaEntry {
    pub name: String,
    pub diff_abs: f64,
    pub diff_rel: f64,
    pub allow_diff: f64,
    pub allow_sign: bool,
    pub check_rel: bool,
    pub num_total: usize,
    pub num_diff_fail: usize,
    pub summary_diff: DeltaLineSummary,
    pub summary_sign_only: DeltaLineSummary,
    pub summary_sign_plus: DeltaLineSummary,
    pub histo: DeltaHistogram,
}

// An object for tracking a series of test results for a the same measurement type,
// recording how they compare to the expected value for the test case, and 
// reporting out those findings.
impl DeltaEntry {
    // todo: Consider adding support for cyclical values, like (-180, 180], where -180 and 180 should be treated as almost equal (but consider it a sign change)
    pub fn new() -> Self {
        DeltaEntry {
            name: "".to_string(),
            allow_diff: 0.0,
            allow_sign: false,
            check_rel: true,
            diff_abs: 0.0,
            diff_rel: 0.0,
            num_total: 0,
            num_diff_fail: 0,
            summary_diff: DeltaLineSummary::new(),
            summary_sign_only: DeltaLineSummary::new(),
            summary_sign_plus: DeltaLineSummary::new(),
            histo: DeltaHistogram::new(),
        }
    }

    // Create a vector of DeltaEntry based on a slice of tuples with the form:
    // (name, allow_diff, check_rel, allow_sign)
    pub fn new_vec(name_base: &str, infos: &[(&str, f64, bool, bool)]) -> Vec<Self> {
        infos.iter().map(|(name, allow_diff, check_rel, allow_sign)| {
            let mut name_full = name_base.to_string();
            name_full.push_str(name);
            DeltaEntry {
                name: name_full,
                allow_diff: *allow_diff,
                check_rel: *check_rel,
                allow_sign: *allow_sign,
                diff_abs: 0.0,
                diff_rel: 0.0,
                num_total: 0,
                num_diff_fail: 0,
                summary_diff: DeltaLineSummary::new(),
                summary_sign_only: DeltaLineSummary::new(),
                summary_sign_plus: DeltaLineSummary::new(),
                histo: DeltaHistogram::new(),
            }
        }).collect()
    }

    // Given x and y, calculate their absolute difference, relative difference,
    // and sign change status, then check whether any of those values is the
    // worst seen so far for comparable operations. If it is, record the line
    // number and the new worst difference.
    // Note that use of this function typically suggests that the values in question
    // are suspected to not be an ideal match for the expected values.
    pub fn add(&mut self, x: f64, y: f64, line_num: usize) {
        self.num_total += 1;
        let (diff_abs, diff_rel) = util::calc_delta(x, y);
        let diff_better = if self.check_rel && util::is_delta_worse(diff_abs, diff_rel) {
            diff_rel
        } else {
            diff_abs
        };
        let is_diff_worst = util::is_delta_worse(diff_better, self.diff_abs) && 
            (!self.check_rel || util::is_delta_worse(diff_better, self.diff_rel));
        // Funky negation on next line is intentional, to get desired nan behavior.
        if !(diff_abs == 0.0) {
            self.summary_diff.add(x, y, line_num, is_diff_worst);
            if is_diff_worst {
                self.diff_abs = diff_abs;
                self.diff_rel = diff_rel;
            }
            if !(diff_better <= self.allow_diff) {
                self.num_diff_fail += 1;
            }
        }
        // For the sign change check, allow (NAN vs NAN), but not (0.0 vs -0.0) or (NAN vs -NAN).
        if x.is_sign_negative() == y.is_sign_negative() {
        } else if diff_abs == 0.0 {
            self.summary_sign_only.add(x, y, line_num, false);
        } else {
            self.summary_sign_plus.add(x, y, line_num, false);
        }
        self.histo.add(diff_better);
    }

    pub fn assert(&self) {
        if self.check_rel {
            crate::assert_delta!(self.summary_diff.x, self.summary_diff.y, self.allow_diff, true, self.name, self.summary_diff.line_num);
        } else {
            crate::assert_delta_abs!(self.summary_diff.x, self.summary_diff.y, self.allow_diff, true, self.name, self.summary_diff.line_num);
        }
        if self.allow_sign {
        } else if self.summary_sign_plus.count > 0 {
            panic!("delta assert failed line {}: sign difference disallowed", self.summary_sign_plus.line_num);
        } else if self.summary_sign_only.count > 0 {
            panic!("delta assert failed line {}: zero vs negative zero sign difference disallowed", self.summary_sign_only.line_num);
        }
    }
}

impl Clone for DeltaEntry {
        fn clone(&self) -> Self {
            DeltaEntry {
                name: self.name.clone(),
                diff_abs: self.diff_abs,
                diff_rel: self.diff_rel,
                allow_diff: self.allow_diff,
                check_rel: self.check_rel,
                allow_sign: self.allow_sign,
                num_total: self.num_total,
                num_diff_fail: self.num_diff_fail,
                summary_diff: self.summary_diff.clone(),
                summary_sign_only: self.summary_sign_only.clone(),
                summary_sign_plus: self.summary_sign_plus.clone(),
                histo: self.histo.clone(),
            }
        }
}

impl Display for DeltaEntry {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::result::Result<(), std::fmt::Error> {
        assert!(self.num_diff_fail <= self.num_total);
        if self.num_total < 1 {
            write!(f, "{}: no diffs requested", self.name)?;
        } else if self.summary_diff.count < 1 {
            write!(f, "{}: no non-zero diffs", self.name)?;
        } else if self.check_rel {
            write!(
                f,
                "{}: worst line {} {}{:e} vs {}{:e} abs {:e}, rel {:e}, {}/{} failed tolerance {:e}",
                self.name,
                self.summary_diff.line_num,
                util::help_sign(self.summary_diff.x),
                self.summary_diff.x,
                util::help_sign(self.summary_diff.y),
                self.summary_diff.y,
                self.diff_abs,
                self.diff_rel,
                self.num_diff_fail,
                self.num_total,
                self.allow_diff,
            )?;
        } else {
            write!(
                f,
                "{}: worst line {} {}{:e} vs {}{:e} abs {:e}, {}/{} failed tolerance {:e}",
                self.name,
                self.summary_diff.line_num,
                util::help_sign(self.summary_diff.x),
                self.summary_diff.x,
                util::help_sign(self.summary_diff.y),
                self.summary_diff.y,
                self.diff_abs,
                self.num_diff_fail,
                self.num_total,
                self.allow_diff,
            )?;
        }
        write!(
            f,
            ", {}+{} sign diffs, histo {}",
            self.summary_sign_only.count,
            self.summary_sign_plus.count,
            self.histo,
        )?;
        Ok(())
    }
}
