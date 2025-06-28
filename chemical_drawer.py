from flask import Flask, render_template, request, jsonify
import requests
import json
import os
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import base64
from io import BytesIO
import re

app = Flask(__name__)

def get_compound_info(compound_name):
    """Get compound information from PubChem API"""
    try:
        # Search for compound by name
        search_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{compound_name}/JSON"
        response = requests.get(search_url)
        
        if response.status_code == 200:
            data = response.json()
            cid = data['PC_Compounds'][0]['id']['id']['cid']
            
            # Get SMILES string
            smiles_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/IsomericSMILES/JSON"
            smiles_response = requests.get(smiles_url)
            
            if smiles_response.status_code == 200:
                smiles_data = smiles_response.json()
                smiles = smiles_data['PropertyTable']['Properties'][0]['IsomericSMILES']
                return {
                    'cid': cid,
                    'smiles': smiles,
                    'name': compound_name,
                    'success': True
                }
        
        return {'success': False, 'error': 'Compound not found'}
    except Exception as e:
        return {'success': False, 'error': str(e)}

def draw_molecule(smiles, compound_name):
    """Draw molecule structure using RDKit"""
    try:
        # Create RDKit molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        
        # Generate 2D coordinates
        AllChem.Compute2DCoords(mol)
        
        # Draw the molecule
        img = Draw.MolToImage(mol, size=(400, 400))
        
        # Convert to base64 for web display
        buffer = BytesIO()
        img.save(buffer, format='PNG')
        img_str = base64.b64encode(buffer.getvalue()).decode()
        
        return img_str
    except Exception as e:
        print(f"Error drawing molecule: {e}")
        return None

def get_compound_properties(cid):
    """Get additional compound properties"""
    try:
        properties_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/MolecularWeight,MolecularFormula,IUPACName/JSON"
        response = requests.get(properties_url)
        
        if response.status_code == 200:
            data = response.json()
            props = data['PropertyTable']['Properties'][0]
            return {
                'molecular_weight': props.get('MolecularWeight', 'N/A'),
                'molecular_formula': props.get('MolecularFormula', 'N/A'),
                'iupac_name': props.get('IUPACName', 'N/A')
            }
        return None
    except:
        return None

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/draw_compound', methods=['POST'])
def draw_compound():
    compound_name = request.form.get('compound_name', '').strip()
    
    if not compound_name:
        return jsonify({'success': False, 'error': 'Please enter a compound name'})
    
    # Get compound information
    compound_info = get_compound_info(compound_name)
    
    if not compound_info['success']:
        return jsonify(compound_info)
    
    # Draw the molecule
    img_base64 = draw_molecule(compound_info['smiles'], compound_info['name'])
    
    if img_base64 is None:
        return jsonify({'success': False, 'error': 'Could not generate structure'})
    
    # Get additional properties
    properties = get_compound_properties(compound_info['cid'])
    
    return jsonify({
        'success': True,
        'image': img_base64,
        'compound_name': compound_info['name'],
        'smiles': compound_info['smiles'],
        'cid': compound_info['cid'],
        'properties': properties
    })

@app.route('/suggest_compounds')
def suggest_compounds():
    query = request.args.get('q', '').lower()
    
    # Common compounds for suggestions
    common_compounds = [
        'water', 'ethanol', 'methane', 'benzene', 'glucose', 'caffeine',
        'aspirin', 'acetaminophen', 'ibuprofen', 'vitamin c', 'sucrose',
        'lactose', 'fructose', 'glycerol', 'acetic acid', 'formic acid',
        'propane', 'butane', 'hexane', 'octane', 'toluene', 'phenol',
        'aniline', 'pyridine', 'thiophene', 'furan', 'pyrrole', 'imidazole'
    ]
    
    suggestions = [comp for comp in common_compounds if query in comp]
    return jsonify(suggestions[:10])

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000) 