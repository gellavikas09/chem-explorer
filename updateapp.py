import numpy as np
import pandas as pd
import py3Dmol
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, Crippen, Lipinski, rdDepictor, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

st.set_page_config(page_title="Chemical Property Explorer", page_icon="🧪", layout="wide")

# ---------------------------------------------------------------
# Example data (approximate literature values, ~20-25 °C, 1 atm)
# ---------------------------------------------------------------
EXAMPLES = {
    "Benzene": {
        "smiles": "c1ccccc1",
        "props": {
            "Class": "Aromatic hydrocarbon",
            "Appearance": "Colorless liquid, sweet odor",
            "Density (g/cm³)": 0.8765, "Melting point (°C)": 5.5, "Boiling point (°C)": 80.1,
            "Flash point (°C)": -11, "Autoignition (°C)": 498, "Water solubility (g/L)": 1.8,
            "Dipole moment (D)": 0.0, "Refractive index": 1.5011, "Viscosity (mPa·s)": 0.604,
            "Vapor pressure (kPa)": 12.7,
            "Key chemistry": "Electrophilic aromatic substitution (nitration, halogenation, Friedel-Crafts)",
            "Uses": "Precursor to styrene, phenol, cyclohexane; solvent",
            "Hazards": "Flammable, carcinogenic (Group 1), bone-marrow toxin",
        },
    },
    "Toluene": {
        "smiles": "Cc1ccccc1",
        "props": {
            "Class": "Aromatic hydrocarbon",
            "Appearance": "Colorless liquid, paint-thinner odor",
            "Density (g/cm³)": 0.8669, "Melting point (°C)": -95, "Boiling point (°C)": 110.6,
            "Flash point (°C)": 4, "Autoignition (°C)": 480, "Water solubility (g/L)": 0.52,
            "Dipole moment (D)": 0.36, "Refractive index": 1.4961, "Viscosity (mPa·s)": 0.560,
            "Vapor pressure (kPa)": 3.8,
            "Key chemistry": "Methyl activates ring (ortho/para); side-chain oxidation to benzoic acid",
            "Uses": "Solvent, gasoline additive, precursor to benzene and TDI",
            "Hazards": "Flammable, CNS depressant, reproductive toxicity",
        },
    },
    "Aniline": {
        "smiles": "Nc1ccccc1",
        "props": {
            "Class": "Aromatic amine",
            "Appearance": "Colorless to brown oily liquid",
            "Density (g/cm³)": 1.0217, "Melting point (°C)": -6.3, "Boiling point (°C)": 184.1,
            "Flash point (°C)": 70, "Autoignition (°C)": 615, "Water solubility (g/L)": 36,
            "Dipole moment (D)": 1.13, "Refractive index": 1.5863, "Viscosity (mPa·s)": 3.85,
            "Vapor pressure (kPa)": 0.04,
            "Key chemistry": "Weak base (pKaH ≈ 4.6); diazotization; strongly activated ring",
            "Uses": "Dyes, MDI/polyurethane precursor, rubber chemicals, pharmaceuticals",
            "Hazards": "Toxic (methemoglobinemia), suspected carcinogen, aquatic toxin",
        },
    },
    "Wurtzite Boron Nitride (w-BN)": {
        "smiles": None,  # crystal, not a molecule
        "props": {
            "Class": "Inorganic superhard ceramic",
            "Appearance": "White to off-white solid",
            "Formula": "BN", "Molar mass (g/mol)": 24.82, "Density (g/cm³)": 3.48,
            "Crystal system": "Hexagonal, space group P6₃mc",
            "Lattice (Å)": "a ≈ 2.55, c ≈ 4.17",
            "Hardness (GPa, Vickers, est.)": "~46 (reported/computed)",
            "Band gap (eV)": "Wide, ~4.5-6 (reported range)",
            "Water solubility": "Insoluble",
            "Key chemistry": "Metastable high-pressure phase; reverts toward h-BN at high temperature",
            "Synthesis": "Shock compression or high pressure/high temperature from h-BN",
            "Uses": "Superhard coatings, abrasives (research), wide-gap electronics",
            "Hazards": "Low toxicity; avoid inhaling fine dust",
        },
    },
}

UPLOAD_FORMATS = {"mol": "sdf", "sdf": "sdf", "pdb": "pdb", "xyz": "xyz", "cif": "cif", "mol2": "mol2"}


# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------
def wurtzite_bn_xyz(nx=3, ny=3, nz=2, a=2.55, c=4.17, u=0.375):
    """Build a wurtzite BN supercell as an XYZ string."""
    a1 = np.array([a, 0.0, 0.0])
    a2 = np.array([-a / 2, a * np.sqrt(3) / 2, 0.0])
    a3 = np.array([0.0, 0.0, c])
    basis = [
        ("B", (1 / 3, 2 / 3, 0.0)), ("B", (2 / 3, 1 / 3, 0.5)),
        ("N", (1 / 3, 2 / 3, u)), ("N", (2 / 3, 1 / 3, 0.5 + u)),
    ]
    atoms = []
    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                for el, (x, y, z) in basis:
                    p = (i + x) * a1 + (j + y) * a2 + (k + z) * a3
                    atoms.append(f"{el} {p[0]:.4f} {p[1]:.4f} {p[2]:.4f}")
    return f"{len(atoms)}\nWurtzite BN supercell\n" + "\n".join(atoms)


def mol_from_smiles_3d(smiles):
    """Return (2D mol, 3D molblock) from a SMILES string, or (None, None)."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    mh = Chem.AddHs(mol)
    if AllChem.EmbedMolecule(mh, randomSeed=42) != 0:
        return mol, None
    try:
        AllChem.MMFFOptimizeMolecule(mh)
    except Exception:
        AllChem.UFFOptimizeMolecule(mh)
    return mol, Chem.MolToMolBlock(mh)


def parse_uploaded_mol(text, ext):
    """Try to parse an uploaded file with RDKit (for 2D drawing and descriptors)."""
    try:
        if ext in ("mol", "sdf"):
            return Chem.MolFromMolBlock(text, removeHs=False)
        if ext == "pdb":
            return Chem.MolFromPDBBlock(text, removeHs=False)
        if ext == "mol2":
            return Chem.MolFromMol2Block(text, removeHs=False)
        if ext == "xyz":
            mol = Chem.MolFromXYZBlock(text)
            if mol is not None:
                try:
                    from rdkit.Chem import rdDetermineBonds
                    rdDetermineBonds.DetermineBonds(mol, charge=0)
                except Exception:
                    pass
            return mol
    except Exception:
        return None
    return None


def draw_2d(mol, size=(450, 350)):
    m = Chem.Mol(mol)
    m = Chem.RemoveHs(m)
    rdDepictor.Compute2DCoords(m)
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    drawer.drawOptions().addStereoAnnotation = True
    drawer.DrawMolecule(m)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def show_3d(data, fmt, style, bg, spin, cell=False, height=480):
    view = py3Dmol.view(width=620, height=height)
    view.addModel(data, fmt)
    styles = {
        "Ball & stick": {"stick": {"radius": 0.12}, "sphere": {"scale": 0.28}},
        "Stick": {"stick": {"radius": 0.15}},
        "Spacefill": {"sphere": {}},
        "Line": {"line": {}},
    }
    view.setStyle(styles[style])
    if cell:
        try:
            view.addUnitCell()
        except Exception:
            pass
    view.setBackgroundColor(bg)
    view.zoomTo()
    if spin:
        view.spin(True)
    components.html(view._make_html(), height=height + 10)


def rdkit_properties(mol):
    m = Chem.Mol(mol)
    props = {
        "Formula": rdMolDescriptors.CalcMolFormula(m),
        "Molar mass (g/mol)": round(Descriptors.MolWt(m), 3),
        "Exact mass": round(Descriptors.ExactMolWt(m), 4),
        "Heavy atoms": m.GetNumHeavyAtoms(),
        "Rings": rdMolDescriptors.CalcNumRings(m),
        "Aromatic rings": rdMolDescriptors.CalcNumAromaticRings(m),
        "logP (Crippen)": round(Crippen.MolLogP(m), 3),
        "TPSA (Å²)": round(rdMolDescriptors.CalcTPSA(m), 2),
        "H-bond donors": Lipinski.NumHDonors(m),
        "H-bond acceptors": Lipinski.NumHAcceptors(m),
        "Rotatable bonds": Lipinski.NumRotatableBonds(m),
        "Canonical SMILES": Chem.MolToSmiles(Chem.RemoveHs(m)),
    }
    try:
        props["InChI"] = Chem.MolToInchi(m)
        props["InChIKey"] = Chem.MolToInchiKey(m)
    except Exception:
        pass
    return props


def show_table(d):
    st.table(pd.DataFrame([(k, str(v)) for k, v in d.items()], columns=["Property", "Value"]))


# ---------------------------------------------------------------
# Sidebar
# ---------------------------------------------------------------
st.title("🧪 Chemical Property Explorer")
st.caption("View 2D and 3D structures, enter SMILES, or upload a structure file.")

with st.sidebar:
    st.header("Input")
    mode = st.radio("Choose input", ["Examples", "SMILES", "Upload file"])
    st.header("3D viewer")
    style = st.selectbox("Style", ["Ball & stick", "Stick", "Spacefill", "Line"])
    bg = st.selectbox("Background", ["white", "black", "#f0f2f6"])
    spin = st.checkbox("Spin", value=False)

mol = None          # RDKit mol for 2D / descriptors
model_data = None   # data for 3D viewer
model_fmt = None
is_crystal = False
extra_props = {}
title = ""

# ---------------------------------------------------------------
# Resolve input
# ---------------------------------------------------------------
if mode == "Examples":
    title = st.selectbox("Compound", list(EXAMPLES))
    entry = EXAMPLES[title]
    extra_props = entry["props"]
    if entry["smiles"]:
        mol, model_data = mol_from_smiles_3d(entry["smiles"])
        model_fmt = "sdf"
        st.write(f"**SMILES:** `{entry['smiles']}`")
    else:
        model_data, model_fmt, is_crystal = wurtzite_bn_xyz(), "xyz", True

elif mode == "SMILES":
    smi = st.text_input("Enter a SMILES string", value="CC(=O)Oc1ccccc1C(=O)O")
    st.caption("Examples: `c1ccccc1` benzene · `CCO` ethanol · `CC(=O)Oc1ccccc1C(=O)O` aspirin")
    title = smi
    if smi.strip():
        mol, model_data = mol_from_smiles_3d(smi.strip())
        model_fmt = "sdf"
        if mol is None:
            st.error("Invalid SMILES. Please check the string.")
        elif model_data is None:
            st.warning("Could not generate a 3D structure for this molecule.")

else:
    up = st.file_uploader("Upload a structure file", type=list(UPLOAD_FORMATS))
    if up is not None:
        ext = up.name.rsplit(".", 1)[-1].lower()
        text = up.getvalue().decode("utf-8", errors="ignore")
        title = up.name
        model_data, model_fmt = text, UPLOAD_FORMATS[ext]
        is_crystal = ext == "cif"
        mol = parse_uploaded_mol(text, ext)
        if mol is None and not is_crystal:
            st.info("3D view is available, but RDKit could not read this file for 2D/descriptors.")
    else:
        st.info("Upload a .mol, .sdf, .pdb, .xyz, .mol2 or .cif file to begin.")

# ---------------------------------------------------------------
# Display
# ---------------------------------------------------------------
if model_data:
    st.subheader(title)
    col2d, col3d = st.columns(2)

    with col2d:
        st.markdown("#### 2D structure")
        if mol is not None:
            try:
                st.image(draw_2d(mol), use_container_width=True)
            except Exception as e:
                st.warning(f"Could not draw 2D structure: {e}")
        else:
            st.info("2D drawing is not available for crystals or unreadable files.")

    with col3d:
        st.markdown("#### 3D structure")
        show_3d(model_data, model_fmt, style, bg, spin, cell=is_crystal)

    tab_a, tab_b, tab_c = st.tabs(["Computed (RDKit)", "Experimental / reference", "Download"])

    with tab_a:
        if mol is not None:
            show_table(rdkit_properties(mol))
        else:
            st.info("No computed descriptors for this input.")

    with tab_b:
        if extra_props:
            show_table(extra_props)
        else:
            st.info("Reference data is only included for the built-in examples.")

    with tab_c:
        ext_map = {"sdf": "sdf", "xyz": "xyz", "pdb": "pdb", "cif": "cif", "mol2": "mol2"}
        st.download_button("Download 3D structure", model_data,
                           file_name=f"structure.{ext_map.get(model_fmt, 'txt')}")
        if mol is not None:
            st.download_button("Download 2D drawing (SVG)", draw_2d(mol), file_name="structure.svg")

st.divider()
st.caption("Built with Streamlit, RDKit and py3Dmol. Reference values are approximate; verify with an SDS or the CRC Handbook.")
