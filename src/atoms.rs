/// Chemical element identified by atomic number.
///
/// Re-exports [`mendeleev::Element`] and provides convenience constructors
/// (`from_symbol`, `from_z`) used throughout the codebase.
pub use mendeleev::Element;

/// Look up an element by its standard symbol (e.g. "Si", "Fe", "C").
///
/// Returns `None` if the symbol does not match any known element.
pub fn from_symbol(s: &str) -> Option<Element> {
    Element::iter().find(|e| e.symbol() == s)
}

/// Look up an element by its atomic number Z (1-based).
///
/// Returns `None` if Z is outside the range 1..=118.
pub fn from_z(z: u32) -> Option<Element> {
    Element::list()
        .iter()
        .copied()
        .find(|e| e.atomic_number() == z)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_atomic_numbers() {
        assert_eq!(Element::H.atomic_number(), 1);
        assert_eq!(Element::Si.atomic_number(), 14);
        assert_eq!(Element::Rb.atomic_number(), 37);
        assert_eq!(Element::Ra.atomic_number(), 88);
        assert_eq!(Element::U.atomic_number(), 92);
    }

    #[test]
    fn test_symbols() {
        assert_eq!(Element::H.symbol(), "H");
        assert_eq!(Element::Fe.symbol(), "Fe");
        assert_eq!(Element::U.symbol(), "U");
    }

    #[test]
    fn test_from_symbol() {
        assert_eq!(from_symbol("Si"), Some(Element::Si));
        assert_eq!(from_symbol("U"), Some(Element::U));
        assert_eq!(from_symbol("Rb"), Some(Element::Rb));
        assert_eq!(from_symbol("Xx"), None);
    }

    #[test]
    fn test_from_z() {
        assert_eq!(from_z(1), Some(Element::H));
        assert_eq!(from_z(92), Some(Element::U));
        assert_eq!(from_z(0), None);
        assert_eq!(from_z(119), None);
    }

    #[test]
    fn test_roundtrip() {
        for z in 1..=92u32 {
            let elem = from_z(z).unwrap();
            assert_eq!(elem.atomic_number(), z);
            let sym = elem.symbol();
            assert_eq!(from_symbol(sym), Some(elem));
        }
    }
}
