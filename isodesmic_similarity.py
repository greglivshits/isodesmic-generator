import sys
import os
import numpy as np
from scipy.optimize import linprog
from rdkit import Chem
from rdkit.Chem import rdFMCS
from rdkit.Chem import DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont
import io
from collections import defaultdict

# --- Configuration ---
INPUT_FILE = "molecules.txt"
OUTPUT_FILE = "isodesmic_structural_similarity.txt"
IMAGE_DIR = "reaction_images_similarity"

if not os.path.exists(IMAGE_DIR):
    os.makedirs(IMAGE_DIR)

# --- Reference Library ---
REF_LIBRARY_SMILES = {
    # --- Aliphatics ---
    'methane': 'C',
    'ethane': 'CC',
    'ethene': 'C=C',
    'ethyne': 'C#C',
    
    # --- Oxygenated ---
    'water': 'O',
    'methanol': 'CO',
    'formaldehyde': 'C=O',
    'dimethyl_ether': 'COC',
    'dioxolane': 'C1COCO1',
    
    # --- Nitrogen ---
    'ammonia': 'N',
    'methylamine': 'CN',
    'dimethylamine': 'CNC', 
    'hydrogen_cyanide': 'C#N',
    'guanidine': 'NC(=N)N',
    'aminoethene': 'C=CN',
    'diaminoethene': 'C=C(N)N',
    
    # --- Nitro / Nitrate / Nitramines ---
    'nitromethane': 'C[N+](=O)[O-]', 
    'ethyl_nitrate': 'CCO[N+](=O)[O-]', 
    'dimethylnitroamine': 'CN(C)[N+](=O)[O-]', 
    'nitroethene': 'C=C[N+]([O-])=O', 
    
    # --- Aromatics ---
    'benzene': 'c1ccccc1',
    'toluene': 'Cc1ccccc1', 
    'phenol': 'Oc1ccccc1', 
    'nitrobenzene': 'c1ccccc1[N+](=O)[O-]', 
    'aniline': 'Nc1ccccc1',
    
    # --- Heterocycles & Benzofused Rings ---
    'benzofuroxan': 'c1ccc2c(c1)[n+]([O-])on2',
    'benzofurazan': 'c1ccc2nonc2c1',
    'furazan': 'c1nonc1',
    'pyridine': 'c1ccncc1',
    
    # --- Nitroso Compounds ---
    'nitrosomethane': 'CN=O',
    'nitrosobenzene': 'O=Nc1ccccc1'
}

# --- Scoring Functions ---

def calculate_similarity_score(target_mol, reactants, products, lib_mols, radius=2):
    """
    Calculates the Average Tanimoto Similarity of UNIQUE participants.
    """
    mfgen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=2048)
    fp_target = mfgen.GetFingerprint(target_mol)
    
    unique_species = set()
    all_participants = reactants + products
    
    for name, coeff in all_participants:
        mol = lib_mols.get(name)
        if not mol: continue
        if mol.HasSubstructMatch(target_mol) and target_mol.HasSubstructMatch(mol):
            continue
        unique_species.add(name)
    
    if not unique_species: return 0.0

    total_sim = 0.0
    for name in unique_species:
        mol = lib_mols[name]
        fp_mol = mfgen.GetFingerprint(mol)
        sim = DataStructs.TanimotoSimilarity(fp_target, fp_mol)
        total_sim += sim
        
    return total_sim / len(unique_species)

def calculate_mcs_score(target_mol, products, lib_mols):
    weighted_score = 0.0
    total_coeffs = 0.0
    params = rdFMCS.MCSParameters()
    params.AtomCompareParameters.RingMatchesRingOnly = False 
    params.BondCompareParameters.RingMatchesRingOnly = False

    for name, coeff in products:
        mol_prod = lib_mols.get(name)
        if not mol_prod: continue
        res = rdFMCS.FindMCS([target_mol, mol_prod], params)
        if mol_prod.GetNumAtoms() > 0:
            ratio = res.numAtoms / mol_prod.GetNumAtoms()
        else:
            ratio = 0.0
        weighted_score += ratio * coeff
        total_coeffs += coeff
        
    if total_coeffs < 1e-6: return 0.0
    return weighted_score / total_coeffs

def get_bond_signature(mol_in, kekulize=False):
    mol = Chem.Mol(mol_in)
    try: Chem.SanitizeMol(mol)
    except: pass
    if kekulize:
        try: Chem.Kekulize(mol, clearAromaticFlags=True)
        except: pass 
    mol = Chem.AddHs(mol)
    features = defaultdict(int)
    for atom in mol.GetAtoms():
        features[atom.GetSymbol()] += 1
    for bond in mol.GetBonds():
        a1, a2 = bond.GetBeginAtom(), bond.GetEndAtom()
        if kekulize:
            btype = 'DOUBLE' if str(bond.GetBondType()) == 'AROMATIC' else str(bond.GetBondType())
            a1_tag, a2_tag = "al", "al"
        else:
            is_ar = a1.GetIsAromatic() and a2.GetIsAromatic()
            btype = "AROMATIC" if is_ar else str(bond.GetBondType())
            a1_tag = "ar" if a1.GetIsAromatic() else "al"
            a2_tag = "ar" if a2.GetIsAromatic() else "al"
        
        desc = sorted([f"{a1.GetSymbol()}({a1_tag})", f"{a2.GetSymbol()}({a2_tag})"])
        features[f"{desc[0]}-{desc[1]}-{btype}"] += 1
    return features

# --- Solver Logic ---

def solve_reaction_multimode(target_mol, library_mols, kekulize=False, num_attempts=500, n_saves=40):
    target_feats = get_bond_signature(target_mol, kekulize=kekulize)
    lib_names = list(library_mols.keys())
    all_keys = set(target_feats.keys())
    lib_feat_maps = []
    for name in lib_names:
        f = get_bond_signature(library_mols[name], kekulize=kekulize)
        lib_feat_maps.append(f)
        all_keys.update(f.keys())
        
    sorted_keys = sorted(list(all_keys))
    n_feats, n_mols = len(sorted_keys), len(lib_names)
    
    A = np.zeros((n_feats, n_mols))
    b = np.zeros(n_feats)
    for i, k in enumerate(sorted_keys):
        b[i] = target_feats.get(k, 0)
        for j, f in enumerate(lib_feat_maps):
            A[i, j] = f.get(k, 0)
            
    A_eq = np.hstack([A, -A])
    bounds = [(0, None) for _ in range(2 * n_mols)]
    unique_solutions = {} 
    
    for _ in range(num_attempts):
        cost_vec = np.random.uniform(0.9, 1.1, size=n_mols)
        n_drop = np.random.randint(0, 8) 
        if n_drop > 0:
            cost_vec[np.random.choice(n_mols, n_drop, replace=False)] = 10000.0 
        
        res = linprog(np.concatenate([cost_vec, cost_vec]), A_eq=A_eq, b_eq=b, bounds=bounds, method='highs')
        
        if res.success:
            u, v = res.x[:n_mols], res.x[n_mols:]
            reactants, products, coeffs_map = [], [], {}
            
            for i in range(n_mols):
                net = u[i] - v[i]
                if abs(net) < 1e-3: continue
                val = round(float(net), 2)
                name = lib_names[i]
                coeffs_map[name] = val
                if val > 0: products.append((name, val))
                else: reactants.append((name, abs(val)))
            
            if not products: continue
            signature = "|".join([f"{n}:{coeffs_map[n]}" for n in sorted(coeffs_map.keys())])
            
            if signature not in unique_solutions:
                cost_val = sum([abs(v) for v in coeffs_map.values()])
                num_comp = len(reactants) + len(products)
                mcs_score = calculate_mcs_score(target_mol, products, library_mols)
                sim_score = calculate_similarity_score(target_mol, reactants, products, library_mols)
                
                # --- NEW MASTER SCORE FORMULA ---
                # 1. Base Score = (10 - N). Gives priority to simple reactions (N=2 -> 8.0, N=3 -> 7.0)
                # 2. Add Sim Score (0.0 - 1.0). Ranks by similarity WITHIN the complexity tier.
                # 3. Subtract tiny cost. Breaks ties.
                iso_score = (10.0 - num_comp) + sim_score - (cost_val * 0.0001)
                
                unique_solutions[signature] = (reactants, products, cost_val, mcs_score, num_comp, sim_score, iso_score)

    if not unique_solutions: return []

    # Sort DESCENDING by the new IsoScore
    all_sols = list(unique_solutions.values())
    all_sols.sort(key=lambda x: -x[6]) 
    
    return all_sols[:n_saves]

# --- Reporting & Drawing ---

def format_reaction_text(target_name, reactants, products):
    lhs = [target_name] + [f"{c:.2f} * {n}" for n, c in reactants]
    rhs = [f"{c:.2f} * {n}" for n, c in products]
    return f"{' + '.join(lhs)} -> {' + '.join(rhs)}"

def get_balance_report(target_mol, reactants, products, library_mols, kekulize=False):
    lhs_totals, rhs_totals = defaultdict(float), defaultdict(float)
    
    t_feats = get_bond_signature(target_mol, kekulize=kekulize)
    for k, v in t_feats.items(): lhs_totals[k] += v
    for n, c in reactants:
        for k, v in get_bond_signature(library_mols[n], kekulize=kekulize).items(): lhs_totals[k] += v * c
    for n, c in products:
        for k, v in get_bond_signature(library_mols[n], kekulize=kekulize).items(): rhs_totals[k] += v * c
            
    lines = [f"\n    --- Balance Verification Table {'(Ring Break)' if kekulize else '(Strict)'} ---"]
    lines.append(f"    {'Feature':<40} | {'Left':<8} | {'Right':<8} | {'Diff'}")
    lines.append("    " + "-" * 70)
    
    is_perfect = True
    for k in sorted(set(lhs_totals) | set(rhs_totals)):
        l, r = lhs_totals.get(k, 0.0), rhs_totals.get(k, 0.0)
        if abs(l - r) > 1e-3: is_perfect = False
        if abs(l) > 1e-4 or abs(r) > 1e-4:
            lines.append(f"    {k:<40} | {l:8.2f} | {r:8.2f} | {l-r:6.3f}")
            
    if not is_perfect: lines.append("    [WARNING] Balance mismatch detected!")
    return "\n".join(lines)

def get_best_font(font_size, bold=False):
    try: return ImageFont.truetype("arialbd.ttf" if bold else "arial.ttf", font_size)
    except: return ImageFont.load_default()

def create_text_image(text, size=(120, 300), font_size=100, bold=True):
    img = Image.new("RGB", size, (255, 255, 255))
    draw = ImageDraw.Draw(img)
    font = get_best_font(font_size, bold)
    bbox = draw.textbbox((0, 0), text, font=font)
    draw.text(((size[0]-bbox[2]+bbox[0])/2, (size[1]-bbox[3]+bbox[1])/2), text, fill="black", font=font)
    return img

def mol_to_pil(mol, legend="", size=(300, 300)):
    if not mol: return Image.new("RGB", size, (255, 255, 255))
    d2d = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    d2d.drawOptions().legendFontSize = 24
    d2d.drawOptions().padding = 0.05
    d2d.DrawMolecule(mol, legend=legend)
    d2d.FinishDrawing()
    return Image.open(io.BytesIO(d2d.GetDrawingText()))

def create_coeff_mol_block(coeff, mol, name, mol_size=(300, 300)):
    mol_img = mol_to_pil(mol, legend=name, size=mol_size)
    coeff_img = create_text_image(f"{coeff:.2f}", size=(160, mol_size[1]), font_size=40, bold=True)
    combined = Image.new("RGB", (coeff_img.width + mol_img.width, max(coeff_img.height, mol_img.height)), (255, 255, 255))
    combined.paste(coeff_img, (0, 0))
    combined.paste(mol_img, (coeff_img.width, 0))
    return combined

def draw_scheme(target_name, target_mol, reactants, products, lib_mols, var_idx, cost, mcs, n_sp, sim, iso_score):
    imgs = [mol_to_pil(target_mol, legend=f"Target: {target_name}")]
    for n, c in reactants:
        imgs.append(create_text_image("+", size=(100,300), font_size=40))
        imgs.append(create_coeff_mol_block(c, lib_mols[n], n))
    imgs.append(create_text_image("---->", size=(200,300), font_size=40))
    first=True
    for n, c in products:
        if not first: imgs.append(create_text_image("+", size=(100,300), font_size=40))
        imgs.append(create_coeff_mol_block(c, lib_mols[n], n))
        first=False
    
    w, h = sum(i.width for i in imgs), 300
    final = Image.new("RGB", (w, h + 50), (255, 255, 255))
    draw = ImageDraw.Draw(final)
    # Updated Info Header with IsoScore
    info_text = f"Var {var_idx+1} | Score: {iso_score:.2f} | Sim: {sim:.2f} | N={n_sp} | Cost: {cost:.2f}"
    draw.text((10, 10), info_text, fill="blue", font=get_best_font(30, bold=True))
    
    x=0
    for i in imgs:
        final.paste(i, (x, (h-i.height)//2 + 50))
        x+=i.width
    
    safe_name = "".join([c for c in target_name if c.isalnum()]).rstrip()
    final.save(os.path.join(IMAGE_DIR, f"{safe_name}_var_{var_idx+1}.png"))

def main():
    lib_mols = {}
    for name, smi in REF_LIBRARY_SMILES.items():
        m = Chem.MolFromSmiles(smi)
        if m: 
            try: Chem.SanitizeMol(m)
            except: pass
            lib_mols[name] = m
    
    results_text = []
    if os.path.exists(INPUT_FILE):
        with open(INPUT_FILE, 'r') as f: lines = f.readlines()
        print(f"Processing {len(lines)} molecules...")
    else:
        print("Demo mode: RDX")
        lines = ["RDX C1(N(CN(CN1[N+](=O)[O-])[N+](=O)[O-])[N+](=O)[O-])"]

    print("-" * 60)
    for line in lines:
        parts = line.strip().split()
        if len(parts) < 2: continue
        mol_name, smiles = parts[0], parts[1]
        target_mol = Chem.MolFromSmiles(smiles)
        if not target_mol: continue
        try: Chem.SanitizeMol(target_mol)
        except: pass

        print(f"Analyzing {mol_name}...")
        
        # 1. Strict
        solutions = solve_reaction_multimode(target_mol, lib_mols, kekulize=False)
        mode = "Strict"
        # 2. Ring Break
        if not solutions:
            print("  > Trying Ring Decomposition...")
            solutions = solve_reaction_multimode(target_mol, lib_mols, kekulize=True)
            mode = "Kekulize"

        if not solutions:
            print(f"{mol_name}: FAILED")
            continue

        print(f"  > Found {len(solutions)} variants (Sorted by IsoScore) ({mode})")
        results_text.append(f"=== {mol_name} ({mode}) ===")
        
        for idx, sol in enumerate(solutions):
            reactants, products, cost, mcs_score, num_comp, sim_score, iso_score = sol
            rxn_str = format_reaction_text(mol_name, reactants, products)
  
            info_str = f"  [Var {idx+1}] Score={iso_score:.2f} | Sim={sim_score:.3f} | MCS={mcs_score:.2f} | N={num_comp} | Cost={cost:.2f}: {rxn_str}"
            print(info_str)
            
            if idx == 0:
                print(get_balance_report(target_mol, reactants, products, lib_mols, kekulize=(mode=="Kekulize")))

            results_text.append(f"  Variant {idx+1} (Score {iso_score:.2f}, Sim {sim_score:.3f}, MCS {mcs_score:.2f}, N={num_comp}, Cost {cost:.2f}):\n  {rxn_str}\n")
            draw_scheme(mol_name, target_mol, reactants, products, lib_mols, idx, cost, mcs_score, num_comp, sim_score, iso_score)
        
        results_text.append("-" * 40 + "\n")

    with open(OUTPUT_FILE, 'w') as f: f.write("\n".join(results_text))
    print(f"Done. Detailed results saved to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()