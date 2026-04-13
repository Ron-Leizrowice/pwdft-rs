/// Chemical element identified by atomic number.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
#[repr(u32)]
pub enum Element {
    H = 1,
    He = 2,
    Li = 3,
    Be = 4,
    B = 5,
    C = 6,
    N = 7,
    O = 8,
    F = 9,
    Ne = 10,
    Na = 11,
    Mg = 12,
    Al = 13,
    Si = 14,
    P = 15,
    S = 16,
    Cl = 17,
    Ar = 18,
    K = 19,
    Ca = 20,
    Sc = 21,
    Ti = 22,
    V = 23,
    Cr = 24,
    Mn = 25,
    Fe = 26,
    Co = 27,
    Ni = 28,
    Cu = 29,
    Zn = 30,
    Ga = 31,
    Ge = 32,
    As = 33,
    Se = 34,
    Br = 35,
    Kr = 36,
    Rb = 37,
    Sr = 38,
    Y = 39,
    Zr = 40,
    Nb = 41,
    Mo = 42,
    Tc = 43,
    Ru = 44,
    Rh = 45,
    Pd = 46,
    Ag = 47,
    Cd = 48,
    In = 49,
    Sn = 50,
    Sb = 51,
    Te = 52,
    I = 53,
    Xe = 54,
    Cs = 55,
    Ba = 56,
    La = 57,
    Ce = 58,
    Pr = 59,
    Nd = 60,
    Pm = 61,
    Sm = 62,
    Eu = 63,
    Gd = 64,
    Tb = 65,
    Dy = 66,
    Ho = 67,
    Er = 68,
    Tm = 69,
    Yb = 70,
    Lu = 71,
    Hf = 72,
    Ta = 73,
    W = 74,
    Re = 75,
    Os = 76,
    Ir = 77,
    Pt = 78,
    Au = 79,
    Hg = 80,
    Tl = 81,
    Pb = 82,
    Bi = 83,
    Po = 84,
    At = 85,
    Rn = 86,
    Fr = 87,
    Ra = 88,
    Ac = 89,
    Th = 90,
    Pa = 91,
    U = 92,
}

impl Element {
    pub fn atomic_number(self) -> u32 {
        self as u32
    }

    pub fn symbol(self) -> &'static str {
        SYMBOLS[self as usize - 1]
    }

    pub fn from_symbol(s: &str) -> Option<Self> {
        SYMBOLS
            .iter()
            .position(|&sym| sym == s)
            .map(|i| Self::from_z(i as u32 + 1).unwrap())
    }

    pub fn from_z(z: u32) -> Option<Self> {
        if z >= 1 && z <= 92 {
            // SAFETY: Element is repr(u32) with contiguous values 1..=92
            Some(unsafe { std::mem::transmute(z) })
        } else {
            None
        }
    }
}

const SYMBOLS: [&str; 92] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne",
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca",
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",
];

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
        assert_eq!(Element::from_symbol("Si"), Some(Element::Si));
        assert_eq!(Element::from_symbol("U"), Some(Element::U));
        assert_eq!(Element::from_symbol("Rb"), Some(Element::Rb));
        assert_eq!(Element::from_symbol("Xx"), None);
    }

    #[test]
    fn test_from_z() {
        assert_eq!(Element::from_z(1), Some(Element::H));
        assert_eq!(Element::from_z(92), Some(Element::U));
        assert_eq!(Element::from_z(0), None);
        assert_eq!(Element::from_z(93), None);
    }

    #[test]
    fn test_roundtrip() {
        for z in 1..=92u32 {
            let elem = Element::from_z(z).unwrap();
            assert_eq!(elem.atomic_number(), z);
            let sym = elem.symbol();
            assert_eq!(Element::from_symbol(sym), Some(elem));
        }
    }
}
