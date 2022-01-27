use pest::iterators::{Pair, Pairs};
use pest::prec_climber::*;
use pest::Parser;

use crate::convert_chart::{convert, UnitType};

#[derive(Parser)]
#[grammar = "grammar.pest"]
struct Calculator;

lazy_static! {
    static ref PREC_CLIMBER: PrecClimber<Rule> = {
        use Assoc::*;
        use Rule::*;

        PrecClimber::new(vec![
            Operator::new(add, Left) | Operator::new(subtract, Left),
            Operator::new(multiply, Left) | Operator::new(divide, Left),
            Operator::new(modulus, Left),
            Operator::new(power, Right),
            Operator::new(percentOf, Left) | Operator::new(percentOn, Left),
            Operator::new(rightShift, Right) | Operator::new(leftShift, Right),
        ])
    };
}

fn eval(expression: Pairs<Rule>) -> f64 {
    PREC_CLIMBER.climb(
        expression,
        |pair: Pair<Rule>| match pair.as_rule() {
            Rule::convert => {
                let mut i = pair.into_inner();
                let value = i.next().unwrap().as_str().parse::<f64>().unwrap();
                // Try to figure out rule name for the conversion between units
                // weight = kilo to gram
                // length = kilometer to meter
                let si_unit_type = i
                    .clone()
                    .next()
                    .unwrap()
                    .into_inner()
                    .next()
                    .unwrap()
                    .as_rule();
                let from = i
                    .next()
                    .unwrap()
                    .into_inner()
                    .next()
                    .unwrap()
                    .into_inner()
                    .next()
                    .unwrap()
                    .as_rule();
                let to = i
                    .next()
                    .unwrap()
                    .into_inner()
                    .next()
                    .unwrap()
                    .into_inner()
                    .next()
                    .unwrap()
                    .as_rule();

                if let (Ok(from), Ok(to)) = (
                    format!("{:?}::{:?}", si_unit_type, from).parse::<UnitType>(),
                    format!("{:?}::{:?}", si_unit_type, to).parse::<UnitType>(),
                ) {
                    convert(value, from, to)
                } else {
                    f64::NAN
                }
            }
            Rule::function => {
                let mut i = pair.into_inner();
                let name = i.next().unwrap().as_str();
                let value = eval(i);
                apply_fun(name, value)
            }
            Rule::pi => std::f64::consts::PI,
            Rule::e => std::f64::consts::E,
            Rule::tau => std::f64::consts::TAU,
            Rule::num => pair.as_str().trim().parse::<f64>().unwrap(),
            Rule::expr => eval(pair.into_inner()),
            _ => f64::NAN,
        },
        |lhs: f64, op: Pair<Rule>, rhs: f64| match op.as_rule() {
            Rule::add => lhs + rhs,
            Rule::subtract => lhs - rhs,
            Rule::multiply => lhs * rhs,
            Rule::divide => lhs / rhs,
            Rule::power => lhs.powf(rhs),
            Rule::percentOf => percent_of(lhs, rhs),
            Rule::percentOn => percent_on(lhs, rhs),
            Rule::rightShift => (lhs as i64 >> rhs as i64) as f64,
            Rule::leftShift => ((lhs as i64) << rhs as i64) as f64,
            Rule::modulus => (lhs % rhs) as f64,
            _ => f64::NAN,
        },
    )
}

fn percent_on(a: f64, b: f64) -> f64 {
    a / 100_f64 * b + b
}

fn percent_of(a: f64, b: f64) -> f64 {
    a / 100_f64 * b
}

fn apply_fun(name: &str, arg: f64) -> f64 {
    match name {
        "sin" => arg.to_radians().sin(),
        "cos" => arg.to_radians().cos(),
        "tan" => arg.to_radians().tan(),
        "asin" => arg.asin(),
        "acos" => arg.cos(),
        "atan" => arg.atan(),
        "sinh" => arg.sinh(),
        "cosh" => arg.cosh(),
        "tanh" => arg.tanh(),
        "asinh" => arg.asinh(),
        "acosh" => arg.acosh(),
        "atanh" => arg.atanh(),
        "log" => arg.log10(),
        "sqrt" => arg.sqrt(),
        "cbrt" => arg.cbrt(),
        "round" => arg.round(),
        "ceil" => arg.ceil(),
        "floor" => arg.floor(),
        _ => f64::NAN,
    }
}

pub fn parse(input: &str) -> f64 {
    let parse_result = Calculator::parse(Rule::calculation, input);
    match parse_result {
        Ok(r) => eval(r),
        Err(_) => f64::NAN,
    }
}


















use newton_rootfinder as nrf;
use nrf::model::Model as nrf_m; // trait import


fn sq3(x: &nalgebra::DVector<f64>) -> nalgebra::DVector<f64> {
    let y = x * x * x * x * x;
    y
}

// Function to optimize: x**2 = 2
pub fn square2(x: &nalgebra::DVector<f64>) -> nalgebra::DVector<f64> {
    let mut y = x * sq3(x);
    y[0] -= 2.0;
   y
}

pub fn solve_rb() -> f64 {

    let problem_size = 1;

    // Parametrization of the iteratives variables
    let vec_iter_params = nrf::iteratives::default_vec_iteratives_fd(problem_size);
    let iter_params = nrf::iteratives::Iteratives::new(&vec_iter_params);

    // Parametrization of the residuals
    let stopping_residuals = vec![nrf::residuals::NormalizationMethod::Abs; problem_size];
    let update_methods = vec![nrf::residuals::NormalizationMethod::Abs; problem_size];
    let res_config = nrf::residuals::ResidualsConfig::new(&stopping_residuals, &update_methods);

    // Parametrization of the solver
    let init = nalgebra::DVector::from_vec(vec![1.0]);
    let resolution_method = nrf::solver::ResolutionMethod::NewtonRaphson;
    let damping = false;
    let mut rf = nrf::solver::default_with_guess(
        init,
        &iter_params,
        &res_config,
        resolution_method,
        damping,
    );

    // Adpatation of the function to solve to the Model trait.
    let mut user_model = nrf::model::UserModelFromFunction::new(problem_size, square2);

    rf.solve(&mut user_model).unwrap();

    // println!("{}", user_model.get_iteratives()[0]);
    user_model.get_iteratives()[0]
}