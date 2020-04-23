load "ke_dbl.gnu"

plot f64_E(x,1.), '+' us ($1-sin($1)):1, f64_E_sgl(x,1.), i32_E(x,1.)

print K, kseq, i32_a
