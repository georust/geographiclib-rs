# geodesic-rs
A port of geographiclib in Rust

This is still an early stage work in progress, if you need something that works I
found [these bindings](https://github.com/savage13/geographiclib) to be usable and fast.

At the moment it's a 1:1 port of the python implementation (that strictly resembles the C++ one, but it's more readable for me).
I think it will stay like this until the port is complete and fully tested, so that any bug caused by errors in the port can be spotted easily.

Then it could be refactored into more idiomatic code, but I'm more concerned about performances.

The inverse algorithm is implemented, but (on my machine) it's 3 times slower than the rust bindings of the C implementation,
and if I want to be able to use it, this should be at least as fast as the bindings. Anyway I think there is room for improvements.

Same for the direct algorithm.

# Roadmap
- [x] Inverse
- [x] Direct

# Try it
If you really want to try it, run
```
cargo run --release
```
to run lots of test cases for ten thousands times for both the rust implementation and the rust bindings in the crate mentioned before, and print the execution times in milliseconds
![Screenshot](images/bench.png)
