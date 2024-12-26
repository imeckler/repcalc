First [install rust](https://www.rust-lang.org/tools/install) if you don't have it already.

To build:

```
cargo build --release
```

To run:

```
./target/release/repcalc
```

Running with `--help` will show the command line options

```
Usage: repcalc [OPTIONS] --precision <PRECISION>

Options:
  -z <x> <y>                       z parameter, x + i y
  -p, --precision <PRECISION>      Number of bits of precision for floating point arithmetic
      --word <WORD>                The word to calculate the value of, a string in {a,b,A,B}
  -r <p> <q>                       Obtain the word by locating the rational p/q in the Stern-Brocot tree
      --random-z                   Use a random value for z
      --random-word <RANDOM_WORD>  Use a uniform random (unreduced) word of the given length
  -h, --help                       Print help
  -V, --version                    Print version
```

So for example, to compute the value of the word aBabb when z is 1 + 2i, with 100 bits of precision, you would run
```
./target/release/repcalc --precision 100 -z 1 2 --word aBabb
```
and obtain
```

(10.554325208519245131314861609014 9.3708473002148348677819571766909e-1) (8.1693283935193764694614700059867e-1 -10.555507141472143962220978984526)
(19.683067160648062353053852999471 -20.944492858527856037779021016004) (-21.054325208519245131314861609216 -19.437084730021483486778195717869)
trace = (-10.500000000000000000000000000202 -18.500000000000000000000000000202)
```
First the rows of the matrix are displayed, and then the trace. A complex number x + i y is represented like (x y).

As an other example, instead of providing the word ababb, you could provide the rational 3/2 corresponding to it via the `-r` option:

```
./target/release/repcalc --precision 100 -z 1 2 -r 3 2
```
to obtain
```
(58.554325208519245131314861609342 16.937084730021483486778195717541) (17.017212478459341774321132216882 -58.455914868184428432016634451547)
(3.4827875215406582256788677830743 26.955914868184428432016634451244) (26.945674791480754868685138391314 -3.4370847300214834867781957178789)
trace = (85.500000000000000000000000000606 13.499999999999999999999999999659)
```
