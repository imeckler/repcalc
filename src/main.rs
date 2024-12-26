use std::cmp::Ordering;
use rand::{rngs::StdRng, SeedableRng, Rng};
use rug::{Complex, Rational, Integer};
use clap::{ArgGroup, ArgAction, Parser};

#[derive(Clone, Copy)]
struct M<A>([A; 4]);

type C = Complex;

fn rho_a(precision: u32, z: C) -> M<C> {
    let one: C = Complex::with_val(precision, 1);
    let c = (z.clone().square() - one.clone()).sqrt().recip();
    let cz = c.clone() * z;
    M([cz.clone(), c.clone(), c, cz])
}

fn rho_b(precision: u32, z: C) -> M<C> {
    let i : C = Complex::with_val(precision, (0, 1)); 
    let one: C = Complex::with_val(precision, 1);
    let y = (-z.clone()) / (z.clone().square() - one.clone()).sqrt();
    let c = (y.clone().square() - one).sqrt().recip();

    let cy = c.clone() * y;
    let ci = c.clone() * i;

    M([cy.clone(), ci.clone(), -ci, cy])
}

impl M<C> {
    fn det(&self) -> C {
        let [a, b, c, d] = &self.0;
        a.clone() * d.clone() - b.clone() * c.clone()
    }

    fn identity(precision: u32) -> Self {
        let one: C = Complex::with_val(precision, 1);
        let zero: C = Complex::with_val(precision, 0);
        M([one.clone(), zero.clone(), zero, one])
    }

    fn inv(self) -> Self {
        let det = &self.det();
        let [a, b, c, d] = self.0;
        M([det.clone()*d, det.clone()*(-b), det.clone()*(-c), det*a])
    }

    fn product(ms: Vec<Self>) -> Self {
        let mut res = ms[0].clone();
        for (i, m) in ms.into_iter().enumerate() {
            if i == 0 { continue }
            res = res.mul(m);
        }
        res
    }

    fn dominant_eigenvector(&self, precision: u32) -> (C, [C; 2]) {
        let two = Complex::with_val(precision, 2);
        let four = Complex::with_val(precision, 4);
        let [a, b, c, d] = &self.0;
        // sqrt(a^2 + 4*b*c - 2*a*d + d^2)
        // using the assumption that det = 1
        // let x = (a.clone().square() + four*b.clone()*c.clone() - two.clone()*a.clone()*d.clone() + d.clone().square()).sqrt();
        let x = ((a.clone() + d.clone()).square() - four.clone()).sqrt();
        // lambda^2 - (a + d) lambda + (ad - bc) = 0
        // lambda = ( (a+d) +/- sqrt((a + d)^2 - 4 (ad - bc)) ) / 2
        let lambda1 = (a.clone() + d.clone() - x.clone()) / two.clone();
        let lambda2 = (a.clone() + d.clone() + x.clone()) / two.clone();
        match lambda1.clone().cmp_abs(&lambda2).unwrap() {
            Ordering::Equal | Ordering::Greater =>
                (lambda1.clone(), [lambda1 - d, c.clone()]),
            Ordering::Less =>
                (lambda2.clone(), [b.clone(), lambda2 - a.clone()])
        }
    }

    fn mul(self, other: Self) -> Self {
        let [a1, b1, c1, d1] = self.0;
        let [a2, b2, c2, d2] = other.0;
        M([
            a1.clone()*a2.clone() + b1.clone()*c2.clone(),
            a1*b2.clone() + b1*d2.clone(),
            c1.clone()*a2 + d1.clone()*c2,
            c1*b2 + d1*d2,
        ])
    }

    fn is_eigenvector(&self, v: [C; 2]) -> bool {
        let [x, y] = v;
        let epsilon = Complex::with_val(x.prec(), 0.000001);
        let [a, b, c, d] = &self.0;
        let ux = a.clone() * x.clone() + b.clone() * y.clone();
        let uy = c.clone() * x.clone() + d.clone() * y.clone();

        let c = ux / x;
        // c * y should be close to uy
        (c * y - uy).cmp_abs(&epsilon).unwrap() == Ordering::Less
    }
}

fn parse_word(input: &str) -> Result<String, String> {
    // Check that every character is one of 'a', 'b', 'A', 'B'
    if input.chars().all(|c| matches!(c, 'a' | 'b' | 'A' | 'B')) {
        Ok(input.to_string())
    } else {
        Err("Value must contain only the letters 'a', 'b', 'A', 'B'.".to_string())
    }
}

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// z parameter, x + i y
    #[arg(
        short,
        // This enforces exactly 2 values
        num_args = 2,
        value_names = ["x", "y"]
    )]
    z: Option<Vec<f64>>,

    /// Number of bits of precision for floating point arithmetic
    #[arg(
        short,
        long,
    )]
    precision: u32,

    /// The word to calculate the value of, a string in {a,b,A,B}
    #[arg(long, value_parser = parse_word)]
    word: Option<String>,

    /// Obtain the word by locating the rational p/q in the Stern-Brocot tree
    #[arg(short, num_args = 2, value_names = ["p", "q"])]
    r: Option<Vec<u64>>,

    /// Use a random value for z
    #[arg(long, action = ArgAction::SetTrue)]
    random_z: bool,

    /// Use a uniform random (unreduced) word of the given length
    #[arg(long)]
    random_word: Option<usize>,
}

enum ExtendedRational {
    R(Rational),
    Infinity,
}

impl PartialEq for ExtendedRational {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (ExtendedRational::Infinity, ExtendedRational::Infinity) => true,
            (ExtendedRational::R(lhs), ExtendedRational::R(rhs)) => lhs == rhs,
            _ => false,
        }
    }
}
impl Eq for ExtendedRational {}
impl PartialOrd for ExtendedRational {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match (self, other) {
            (ExtendedRational::Infinity, ExtendedRational::Infinity) => Some(Ordering::Equal),
            (ExtendedRational::Infinity, ExtendedRational::R(_)) => Some(Ordering::Greater),
            (ExtendedRational::R(_), ExtendedRational::Infinity) => Some(Ordering::Less),
            (ExtendedRational::R(lhs), ExtendedRational::R(rhs)) => lhs.partial_cmp(rhs),
        }
    }
}
impl Ord for ExtendedRational {
    fn cmp(&self, other: &Self) -> Ordering {
        // Because partial_cmp always returns Some(...), it's safe to unwrap
        self.partial_cmp(other).unwrap()
    }
}

const ZERO: &'static Integer = & Integer::ZERO;

impl ExtendedRational {
    fn numer(&self) -> &Integer {
        match self {
            ExtendedRational::R(r) => r.numer(),
            ExtendedRational::Infinity => Integer::ONE,
        }
    }

    fn denom(&self) -> &Integer {
        match self {
            ExtendedRational::R(r) => r.denom(),
            ExtendedRational::Infinity => ZERO,
        }
    }

    fn mediant(&self, other: &Self) -> Self {
        let x = self.numer().clone() + other.numer().clone();
        let y = self.denom().clone() + other.denom().clone();
        if y.is_zero() {
            ExtendedRational::Infinity
        } else {
            ExtendedRational::R(Rational::from((x, y)))
        }
    }
}

fn stern_brocot_word(q: ExtendedRational, a: M<C>, b: M<C>) -> M<C> {
    match &q {
        ExtendedRational::Infinity => { return b },
        ExtendedRational::R(x) => {
            if x.eq(Rational::ONE) {
                return a;
            }
        }
    }

    let mut low = ExtendedRational::R(Rational::ZERO.clone());
    let mut low_m = a;
    let mut high = ExtendedRational::Infinity;
    let mut high_m = b;

    loop {
        let med = low.mediant(&high);
        if med < q {
            // q is in (med, high)
            low = med;
            low_m = low_m.mul(high_m.clone());
        } else if q < med {
            // q is in (low, med)
            high = med;
            high_m = low_m.clone().mul(high_m)
        } else {
            // finished
            return low_m.mul(high_m)
        }
    }
}

fn main() {
    let args = Args::parse();
    let precision = args.precision;
    let rng = &mut StdRng::from_seed([2u8; 32]);
    
    let z: C =
        if args.random_z {
            Complex::with_val(precision, (rng.gen::<f64>(), rng.gen::<f64>()))
        } else if let Some(z) = args.z {
            Complex::with_val(precision, (z[0], z[1]))
        } else {
            eprintln!("At least one of z, random-z must be provided.");
            std::process::exit(1)
        };

    let a = rho_a(precision, z.clone());
    let b = rho_b(precision, z);
    let a_inv = a.clone().inv();
    let b_inv = b.clone().inv();

    let res =
        if let Some(n) = args.random_word {
            M::product((0..n).map(|_| 
                match rng.gen_range(0usize..4) {
                    0 => a.clone(),
                    1 => b.clone(),
                    2 => a_inv.clone(),
                    3 => b_inv.clone(),
                    _ => panic!("impossible")
                }).collect())
        } else if let Some(word) = args.word {
            M::product(word.chars().map(|c|
                match c {
                    'a' => a.clone(),
                    'b' => b.clone(),
                    'A' => a_inv.clone(),
                    'B' => b_inv.clone(),
                    _ => panic!("impossible")
                }).collect())
        } else if let Some(r) = args.r {
            let p = r[0];
            let q = r[1];
            let x =
                if q == 0 {
                    ExtendedRational::Infinity
                } else {
                    ExtendedRational::R(Rational::from((p, q)))
                };
            stern_brocot_word(x, a, b)
        } else {
            eprintln!("At least one of --word, --random-word, -r must be provided.");
            std::process::exit(1);
        };

    let [x, y, z, w] = &res.0;
    println!("{} {}\n{} {}", x.clone(), y.clone(), z.clone(), w.clone());
    println!("trace = {}", x.clone() + w.clone());
    let (lambda, [vx, vy]) = res.dominant_eigenvector(precision);
    if !res.is_eigenvector([vx.clone(), vy.clone()]) {
        eprintln!("warning: output is not very close to an eigenvector, increase precision")
    }
    println!("dominant_eigenvalue = {}", lambda);
    println!("dominant_eigenvector = {} {}", vx, vy);
}
