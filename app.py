import pandas as pd
import streamlit as st

st.set_page_config(page_title="Chemical Property Explorer", page_icon="🧪", layout="wide")

# ---------------------------------------------------------------
# Data (approximate literature values at ~20-25 °C, 1 atm)
# Add more compounds by adding another dict entry.
# ---------------------------------------------------------------
COMPOUNDS = {
    "Benzene": {
        "Formula": "C6H6",
        "Class": "Aromatic hydrocarbon",
        "Appearance": "Colorless liquid, sweet odor",
        "Molar mass (g/mol)": 78.11,
        "Density (g/cm³)": 0.8765,
        "Melting point (°C)": 5.5,
        "Boiling point (°C)": 80.1,
        "Flash point (°C)": -11,
        "Autoignition (°C)": 498,
        "Water solubility (g/L)": 1.8,
        "Dipole moment (D)": 0.0,
        "Refractive index": 1.5011,
        "Viscosity (mPa·s)": 0.604,
        "log P": 2.13,
        "Vapor pressure (kPa)": 12.7,
        "Crystal / structure": "Planar hexagonal ring, D6h symmetry",
        "Key chemistry": "Electrophilic aromatic substitution (nitration, halogenation, Friedel-Crafts)",
        "Uses": "Precursor to styrene, phenol, cyclohexane; solvent",
        "Hazards": "Flammable, carcinogenic (Group 1), toxic to bone marrow",
    },
    "Toluene": {
        "Formula": "C7H8",
        "Class": "Aromatic hydrocarbon",
        "Appearance": "Colorless liquid, paint-thinner odor",
        "Molar mass (g/mol)": 92.14,
        "Density (g/cm³)": 0.8669,
        "Melting point (°C)": -95,
        "Boiling point (°C)": 110.6,
        "Flash point (°C)": 4,
        "Autoignition (°C)": 480,
        "Water solubility (g/L)": 0.52,
        "Dipole moment (D)": 0.36,
        "Refractive index": 1.4961,
        "Viscosity (mPa·s)": 0.560,
        "log P": 2.73,
        "Vapor pressure (kPa)": 3.8,
        "Crystal / structure": "Benzene ring with one methyl group",
        "Key chemistry": "Methyl group activates ring (ortho/para directing); side-chain oxidation to benzoic acid",
        "Uses": "Solvent, gasoline additive, precursor to benzene, TDI",
        "Hazards": "Flammable, CNS depressant, reproductive toxicity",
    },
    "Aniline": {
        "Formula": "C6H7N",
        "Class": "Aromatic amine",
        "Appearance": "Colorless to brown oily liquid",
        "Molar mass (g/mol)": 93.13,
        "Density (g/cm³)": 1.0217,
        "Melting point (°C)": -6.3,
        "Boiling point (°C)": 184.1,
        "Flash point (°C)": 70,
        "Autoignition (°C)": 615,
        "Water solubility (g/L)": 36,
        "Dipole moment (D)": 1.13,
        "Refractive index": 1.5863,
        "Viscosity (mPa·s)": 3.85,
        "log P": 0.90,
        "Vapor pressure (kPa)": 0.04,
        "Crystal / structure": "Benzene ring with an -NH2 group",
        "Key chemistry": "Weak base (pKaH ≈ 4.6); diazotization to diazonium salts; strongly activated ring",
        "Uses": "Dyes, polyurethane (MDI) precursor, rubber chemicals, pharmaceuticals",
        "Hazards": "Toxic (methemoglobinemia), suspected carcinogen, environmental hazard",
    },
    "Wurtzite Boron Nitride (w-BN)": {
        "Formula": "BN",
        "Class": "Inorganic ceramic / superhard material",
        "Appearance": "White to off-white solid",
        "Molar mass (g/mol)": 24.82,
        "Density (g/cm³)": 3.48,
        "Melting point (°C)": None,  # no clean melting point at 1 atm
        "Boiling point (°C)": None,
        "Flash point (°C)": None,
        "Autoignition (°C)": None,
        "Water solubility (g/L)": 0.0,
        "Dipole moment (D)": None,
        "Refractive index": None,
        "Viscosity (mPa·s)": None,
        "log P": None,
        "Vapor pressure (kPa)": None,
        "Crystal / structure": "Hexagonal, space group P6₃mc; a ≈ 2.55 Å, c ≈ 4.2 Å; sp³ bonded",
        "Key chemistry": "Metastable high-pressure phase; reverts toward hexagonal BN at high temperature",
        "Uses": "Superhard coatings and abrasives (research), wide-gap electronics",
        "Hazards": "Low toxicity; avoid inhaling fine dust",
        "Hardness (GPa, Vickers, est.)": "~46 (reported/computed)",
        "Band gap (eV)": "Wide, ~4.5-6 (reported range)",
        "Synthesis": "Shock compression or high pressure/high temperature from h-BN",
    },
}

NUMERIC = [
    "Molar mass (g/mol)", "Density (g/cm³)", "Melting point (°C)", "Boiling point (°C)",
    "Flash point (°C)", "Autoignition (°C)", "Water solubility (g/L)",
    "Dipole moment (D)", "Refractive index", "Viscosity (mPa·s)", "log P", "Vapor pressure (kPa)",
]

# ---------------------------------------------------------------
# UI
# ---------------------------------------------------------------
st.title("🧪 Chemical Property Explorer")
st.caption("Select a compound to view its physical and chemical properties. Values are approximate literature data.")

tab1, tab2 = st.tabs(["Single compound", "Compare compounds"])

with tab1:
    name = st.selectbox("Choose a compound", list(COMPOUNDS))
    data = COMPOUNDS[name]

    st.subheader(f"{name}  ·  {data['Formula']}")
    st.write(f"**Class:** {data['Class']}  |  **Appearance:** {data['Appearance']}")

    c1, c2, c3, c4 = st.columns(4)
    c1.metric("Molar mass", f"{data['Molar mass (g/mol)']} g/mol")
    c2.metric("Density", f"{data['Density (g/cm³)']} g/cm³")
    mp, bp = data["Melting point (°C)"], data["Boiling point (°C)"]
    c3.metric("Melting point", f"{mp} °C" if mp is not None else "n/a")
    c4.metric("Boiling point", f"{bp} °C" if bp is not None else "n/a")

    left, right = st.columns(2)
    with left:
        st.markdown("#### Physical properties")
        phys = {k: v for k, v in data.items() if k in NUMERIC and v is not None}
        st.table(pd.DataFrame([(k, str(v)) for k, v in phys.items()], columns=["Property", "Value"]))
    with right:
        st.markdown("#### Chemical / other properties")
        for key in ["Crystal / structure", "Key chemistry", "Uses", "Hazards",
                    "Hardness (GPa, Vickers, est.)", "Band gap (eV)", "Synthesis"]:
            if key in data:
                st.markdown(f"**{key}:** {data[key]}")

with tab2:
    picks = st.multiselect("Select compounds", list(COMPOUNDS), default=list(COMPOUNDS)[:3])
    if picks:
        df = pd.DataFrame({p: {k: COMPOUNDS[p][k] for k in NUMERIC} for p in picks}).astype(float)
        st.dataframe(df, use_container_width=True)

        prop = st.selectbox("Plot property", NUMERIC, index=3)
        series = df.loc[prop].dropna()
        if series.empty:
            st.info("No data available for this property in the selected compounds.")
        else:
            st.bar_chart(series)
    else:
        st.info("Pick at least one compound.")

st.divider()
st.caption("Built with Streamlit. Verify data against an SDS or the CRC Handbook before real use.")
