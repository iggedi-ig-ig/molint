macro_rules! define_atoms {
    (
        $(#[derive($($derive:ident),*)])?
        $vis:vis enum $enum_name:ident {
            $(
                #[doc = $doc:expr]
                $name:ident = $ordinal:expr
            ),*
        }
    ) => {
        $(#[derive($($derive),*)])?
        #[repr(usize)]
        $vis enum $enum_name {
            $(
                #[doc = $doc]
                $name = $ordinal
            ),*
        }

        #[allow(unused)]
        impl $enum_name {
            $vis fn from_ordinal(ordinal: usize) -> Option<Self> {
                Some(match ordinal {
                    $(
                        $ordinal => Self::$name,
                    )*
                    _ => return None,
                })
            }

            $vis fn from_symbol(symbol: &str) -> Option<Self> {
                Some(match symbol {
                    $(
                        s if s == stringify!($name) => Self::$name,
                    )*
                    _ => return None,
                })
            }

            $vis fn symbol(&self) -> &str {
                match self {
                    $(
                        Self::$name => stringify!($name)
                    ),*
                }
            }

            fn ordinal_str(&self) -> &str {
                match self {
                    $(
                        Self::$name => stringify!($ordinal)
                    ),*
                }
            }

            fn from_ordinal_str(ordinal_str: &str) -> Option<Self> {
                Some(match ordinal_str.parse::<usize>() {
                    $(
                        Ok(i) if i == $ordinal => Self::$name,
                    )*
                    _ => return None,
                })
            }
        }

        impl serde::Serialize for $enum_name {
            fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
            where
                S: serde::Serializer,
            {
                serializer.serialize_str(self.ordinal_str())
            }
        }

        impl<'de> serde::Deserialize<'de> for $enum_name {
            fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
            where
                D: serde::Deserializer<'de>,
            {
                let s = String::deserialize(deserializer)?;
                Self::from_ordinal_str(&s).ok_or_else(|| serde::de::Error::custom(format!("unknown atom ordinal {}", s)))
            }
        }
    };
}

define_atoms! {
    #[derive(Copy, Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
    pub enum ElementType {
        #[doc = "Hydrogen"]
        H = 1,
        #[doc = "Helium"]
        He = 2,
        #[doc = "Lithium"]
        Li = 3,
        #[doc = "Beryllium"]
        Be = 4,
        #[doc = "Boron"]
        B = 5,
        #[doc = "Carbon"]
        C = 6,
        #[doc = "Nitrogen"]
        N = 7,
        #[doc = "Oxygen"]
        O = 8,
        #[doc = "Fluorine"]
        F = 9,
        #[doc = "Neon"]
        Ne = 10,
        #[doc = "Sodium"]
        Na = 11,
        #[doc = "Magnesium"]
        Mg = 12,
        #[doc = "Aluminum"]
        Al = 13,
        #[doc = "Silicon"]
        Si = 14,
        #[doc = "Phosphorus"]
        P = 15,
        #[doc = "Sulfur"]
        S = 16,
        #[doc = "Chlorine"]
        Cl = 17,
        #[doc = "Argon"]
        Ar = 18,
        #[doc = "Potassium"]
        K = 19,
        #[doc = "Calcium"]
        Ca = 20,
        #[doc = "Scandium"]
        Sc = 21,
        #[doc = "Titanium"]
        Ti = 22,
        #[doc = "Vanadium"]
        V = 23,
        #[doc = "Chromium"]
        Cr = 24,
        #[doc = "Manganese"]
        Mn = 25,
        #[doc = "Iron"]
        Fe = 26,
        #[doc = "Cobalt"]
        Co = 27,
        #[doc = "Nickel"]
        Ni = 28,
        #[doc = "Copper"]
        Cu = 29,
        #[doc = "Zinc"]
        Zn = 30,
        #[doc = "Gallium"]
        Ga = 31,
        #[doc = "Germanium"]
        Ge = 32,
        #[doc = "Arsenic"]
        As = 33,
        #[doc = "Selenium"]
        Se = 34,
        #[doc = "Bromine"]
        Br = 35,
        #[doc = "Krypton"]
        Kr = 36,
        #[doc = "Rubidium"]
        Rb = 37,
        #[doc = "Strontium"]
        Sr = 38,
        #[doc = "Yttrium"]
        Y = 39,
        #[doc = "Zirconium"]
        Zr = 40,
        #[doc = "Niobium"]
        Nb = 41,
        #[doc = "Molybdenum"]
        Mo = 42,
        #[doc = "Technetium"]
        Tc = 43,
        #[doc = "Ruthenium"]
        Ru = 44,
        #[doc = "Rhodium"]
        Rh = 45,
        #[doc = "Palladium"]
        Pd = 46,
        #[doc = "Silver"]
        Ag = 47,
        #[doc = "Cadmium"]
        Cd = 48,
        #[doc = "Indium"]
        In = 49,
        #[doc = "Tin"]
        Sn = 50,
        #[doc = "Antimony"]
        Sb = 51,
        #[doc = "Tellurium"]
        Te = 52,
        #[doc = "Iodine"]
        I = 53,
        #[doc = "Xenon"]
        Xe = 54,
        #[doc = "Cesium"]
        Cs = 55,
        #[doc = "Barium"]
        Ba = 56,
        #[doc = "Lanthanum"]
        La = 57,
        #[doc = "Cerium"]
        Ce = 58,
        #[doc = "Praseodymium"]
        Pr = 59,
        #[doc = "Neodymium"]
        Nd = 60,
        #[doc = "Promethium"]
        Pm = 61,
        #[doc = "Samarium"]
        Sm = 62,
        #[doc = "Europium"]
        Eu = 63,
        #[doc = "Gadolinium"]
        Gd = 64,
        #[doc = "Terbium"]
        Tb = 65,
        #[doc = "Dysprosium"]
        Dy = 66,
        #[doc = "Holmium"]
        Ho = 67,
        #[doc = "Erbium"]
        Er = 68,
        #[doc = "Thulium"]
        Tm = 69,
        #[doc = "Ytterbium"]
        Yb = 70,
        #[doc = "Lutetium"]
        Lu = 71,
        #[doc = "Hafnium"]
        Hf = 72,
        #[doc = "Tantalum"]
        Ta = 73,
        #[doc = "Tungsten"]
        W = 74,
        #[doc = "Rhenium"]
        Re = 75,
        #[doc = "Osmium"]
        Os = 76,
        #[doc = "Iridium"]
        Ir = 77,
        #[doc = "Platinum"]
        Pt = 78,
        #[doc = "Gold"]
        Au = 79,
        #[doc = "Mercury"]
        Hg = 80,
        #[doc = "Thallium"]
        Tl = 81,
        #[doc = "Lead"]
        Pb = 82,
        #[doc = "Bismuth"]
        Bi = 83,
        #[doc = "Polonium"]
        Po = 84,
        #[doc = "Astatine"]
        At = 85,
        #[doc = "Radon"]
        Rn = 86,
        #[doc = "Francium"]
        Fr = 87,
        #[doc = "Radium"]
        Ra = 88,
        #[doc = "Actinium"]
        Ac = 89,
        #[doc = "Thorium"]
        Th = 90,
        #[doc = "Protactinium"]
        Pa = 91,
        #[doc = "Uranium"]
        U = 92,
        #[doc = "Neptunium"]
        Np = 93,
        #[doc = "Plutonium"]
        Pu = 94,
        #[doc = "Americium"]
        Am = 95,
        #[doc = "Curium"]
        Cm = 96,
        #[doc = "Berkelium"]
        Bk = 97,
        #[doc = "Californium"]
        Cf = 98,
        #[doc = "Einsteinium"]
        Es = 99,
        #[doc = "Fermium"]
        Fm = 100,
        #[doc = "Mendelevium"]
        Md = 101,
        #[doc = "Nobelium"]
        No = 102,
        #[doc = "Lawrencium"]
        Lr = 103,
        #[doc = "Rutherfordium"]
        Rf = 104,
        #[doc = "Dubnium"]
        Db = 105,
        #[doc = "Seaborgium"]
        Sg = 106,
        #[doc = "Bohrium"]
        Bh = 107,
        #[doc = "Hassium"]
        Hs = 108,
        #[doc = "Meitnerium"]
        Mt = 109,
        #[doc = "Darmstadtium"]
        Ds = 110,
        #[doc = "Roentgenium"]
        Rg = 111,
        #[doc = "Copernicium"]
        Cn = 112,
        #[doc = "Nihonium"]
        Nh = 113,
        #[doc = "Flerovium"]
        Fl = 114,
        #[doc = "Moscovium"]
        Mc = 115,
        #[doc = "Livermorium"]
        Lv = 116,
        #[doc = "Tennessine"]
        Ts = 117,
        #[doc = "Oganesson"]
        Og = 118
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let hydrogen: Result<ElementType, _> = serde_json::from_str(r#""1""#);
        let helium: Result<ElementType, _> = serde_json::from_str(r#""2""#);
        let uranium: Result<ElementType, _> = serde_json::from_str(r#""92""#);

        assert!(matches!(hydrogen, Ok(ElementType::H)));
        assert!(matches!(helium, Ok(ElementType::He)));
        assert!(matches!(uranium, Ok(ElementType::U)));
    }
}
