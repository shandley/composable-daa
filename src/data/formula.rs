//! Formula parsing for model specification.

use crate::error::{DaaError, Result};
use serde::{Deserialize, Serialize};

/// A term in a formula.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub enum Term {
    /// Intercept term (constant).
    Intercept,
    /// Main effect of a variable.
    Main(String),
    /// Interaction between two variables.
    Interaction(String, String),
}

impl Term {
    /// Get the variable names involved in this term.
    pub fn variables(&self) -> Vec<&str> {
        match self {
            Term::Intercept => vec![],
            Term::Main(v) => vec![v.as_str()],
            Term::Interaction(v1, v2) => vec![v1.as_str(), v2.as_str()],
        }
    }
}

impl std::fmt::Display for Term {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Term::Intercept => write!(f, "1"),
            Term::Main(v) => write!(f, "{}", v),
            Term::Interaction(v1, v2) => write!(f, "{}:{}", v1, v2),
        }
    }
}

/// A parsed formula specifying a linear model.
///
/// Supports R-style formula syntax:
/// - `~ group` - intercept + group
/// - `~ group + age` - intercept + group + age
/// - `~ group * age` - intercept + group + age + group:age
/// - `~ 0 + group` - no intercept, group
/// - `~ group:age` - intercept + interaction only
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Formula {
    /// Whether to include an intercept.
    pub intercept: bool,
    /// Terms in the formula (excluding intercept).
    pub terms: Vec<Term>,
    /// Original formula string.
    pub formula_str: String,
}

impl Formula {
    /// Parse a formula string.
    ///
    /// # Examples
    /// ```
    /// use composable_daa::data::Formula;
    /// let f = Formula::parse("~ group + age").unwrap();
    /// assert!(f.intercept);
    /// assert_eq!(f.terms.len(), 2);
    /// ```
    pub fn parse(formula: &str) -> Result<Self> {
        let formula_str = formula.to_string();
        let formula = formula.trim();

        // Must start with ~
        if !formula.starts_with('~') {
            return Err(DaaError::FormulaParse(
                "Formula must start with '~'".to_string(),
            ));
        }

        let rhs = formula[1..].trim();
        if rhs.is_empty() {
            return Err(DaaError::FormulaParse(
                "Formula right-hand side is empty".to_string(),
            ));
        }

        // Check for explicit no-intercept
        let (intercept, rhs) = if rhs == "0" || rhs == "-1" {
            // Just "0" or "-1" alone means no-intercept with no terms - invalid
            return Err(DaaError::FormulaParse(
                "Formula must have at least one term".to_string(),
            ));
        } else if rhs.starts_with("0 +") || rhs.starts_with("0+") {
            (false, rhs.trim_start_matches("0").trim_start_matches('+').trim())
        } else if rhs.starts_with("-1 +") || rhs.starts_with("-1+") {
            (false, rhs.trim_start_matches("-1").trim_start_matches('+').trim())
        } else {
            (true, rhs)
        };

        // Parse terms
        let mut terms = Vec::new();
        let term_strs: Vec<&str> = rhs.split('+').map(|s| s.trim()).collect();

        for term_str in term_strs {
            if term_str.is_empty() {
                continue;
            }

            // Check for * (expansion to main effects + interaction)
            if term_str.contains('*') {
                let parts: Vec<&str> = term_str.split('*').map(|s| s.trim()).collect();
                if parts.len() != 2 {
                    return Err(DaaError::FormulaParse(format!(
                        "Invalid interaction term: {}",
                        term_str
                    )));
                }
                let v1 = parts[0].to_string();
                let v2 = parts[1].to_string();

                // Add main effects if not already present
                let main1 = Term::Main(v1.clone());
                let main2 = Term::Main(v2.clone());
                if !terms.contains(&main1) {
                    terms.push(main1);
                }
                if !terms.contains(&main2) {
                    terms.push(main2);
                }
                terms.push(Term::Interaction(v1, v2));
            }
            // Check for : (interaction only)
            else if term_str.contains(':') {
                let parts: Vec<&str> = term_str.split(':').map(|s| s.trim()).collect();
                if parts.len() != 2 {
                    return Err(DaaError::FormulaParse(format!(
                        "Invalid interaction term: {}",
                        term_str
                    )));
                }
                terms.push(Term::Interaction(
                    parts[0].to_string(),
                    parts[1].to_string(),
                ));
            }
            // Main effect
            else {
                // Skip "1" as explicit intercept
                if term_str != "1" {
                    terms.push(Term::Main(term_str.to_string()));
                }
            }
        }

        if terms.is_empty() && !intercept {
            return Err(DaaError::FormulaParse(
                "Formula must have at least one term".to_string(),
            ));
        }

        Ok(Self {
            intercept,
            terms,
            formula_str,
        })
    }

    /// Get all variable names used in the formula.
    pub fn variables(&self) -> Vec<&str> {
        let mut vars: Vec<&str> = self
            .terms
            .iter()
            .flat_map(|t| t.variables())
            .collect();
        vars.sort();
        vars.dedup();
        vars
    }

    /// Check if a variable is used in the formula.
    pub fn uses_variable(&self, name: &str) -> bool {
        self.terms.iter().any(|t| t.variables().contains(&name))
    }
}

impl std::fmt::Display for Formula {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "~ ")?;
        if !self.intercept {
            write!(f, "0 + ")?;
        }
        let term_strs: Vec<String> = self.terms.iter().map(|t| t.to_string()).collect();
        write!(f, "{}", term_strs.join(" + "))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_simple() {
        let f = Formula::parse("~ group").unwrap();
        assert!(f.intercept);
        assert_eq!(f.terms.len(), 1);
        assert_eq!(f.terms[0], Term::Main("group".to_string()));
    }

    #[test]
    fn test_parse_multiple() {
        let f = Formula::parse("~ group + age").unwrap();
        assert!(f.intercept);
        assert_eq!(f.terms.len(), 2);
        assert_eq!(f.terms[0], Term::Main("group".to_string()));
        assert_eq!(f.terms[1], Term::Main("age".to_string()));
    }

    #[test]
    fn test_parse_no_intercept() {
        let f = Formula::parse("~ 0 + group").unwrap();
        assert!(!f.intercept);
        assert_eq!(f.terms.len(), 1);
    }

    #[test]
    fn test_parse_interaction() {
        let f = Formula::parse("~ group:age").unwrap();
        assert!(f.intercept);
        assert_eq!(f.terms.len(), 1);
        assert_eq!(
            f.terms[0],
            Term::Interaction("group".to_string(), "age".to_string())
        );
    }

    #[test]
    fn test_parse_star_expansion() {
        let f = Formula::parse("~ group * age").unwrap();
        assert!(f.intercept);
        assert_eq!(f.terms.len(), 3);
        assert!(f.terms.contains(&Term::Main("group".to_string())));
        assert!(f.terms.contains(&Term::Main("age".to_string())));
        assert!(f.terms.contains(&Term::Interaction(
            "group".to_string(),
            "age".to_string()
        )));
    }

    #[test]
    fn test_variables() {
        let f = Formula::parse("~ group + age + group:age").unwrap();
        let vars = f.variables();
        assert_eq!(vars, vec!["age", "group"]);
    }

    #[test]
    fn test_invalid_formula() {
        assert!(Formula::parse("group + age").is_err()); // missing ~
        assert!(Formula::parse("~").is_err()); // empty RHS
        assert!(Formula::parse("~ 0").is_err()); // no terms
    }
}
