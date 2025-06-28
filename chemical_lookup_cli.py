import requests
import json
from PIL import Image
import io
import os
import sys

class ChemicalCompoundLookupCLI:
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    
    def search_compound(self, compound_name):
        """Search for a compound by name and return information"""
        print(f"Searching for: {compound_name}")
        
        try:
            # First, search for the compound by name
            search_url = f"{self.base_url}/compound/name/{compound_name}/cids/JSON"
            response = requests.get(search_url, timeout=10)
            
            if response.status_code != 200:
                # Try fuzzy search
                print("Exact match not found, trying fuzzy search...")
                search_url = f"{self.base_url}/compound/name/{compound_name}/cids/JSON?name_type=word"
                response = requests.get(search_url, timeout=10)
                
                if response.status_code != 200:
                    print("Error: Compound not found. Please try a different name.")
                    return None
            
            data = response.json()
            
            # Try both old and new API response formats
            cids = []
            if 'InformationList' in data:
                cids = data.get('InformationList', {}).get('Information', [{}])[0].get('CID', [])
            elif 'IdentifierList' in data:
                cids = data.get('IdentifierList', {}).get('CID', [])
            
            if not cids:
                print("Error: No compounds found with that name.")
                return None
            
            # Get the first (best match) CID
            cid = cids[0]
            print(f"Found compound with CID: {cid}")
            
            # Get compound information
            info_url = f"{self.base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
            info_response = requests.get(info_url, timeout=10)
            
            if info_response.status_code != 200:
                print("Error: Could not retrieve compound information.")
                return None
            
            info_data = info_response.json()
            properties = info_data.get('PropertyTable', {}).get('Properties', [{}])[0]
            
            # Get compound image
            image_url = f"{self.base_url}/compound/cid/{cid}/PNG"
            image_response = requests.get(image_url, timeout=10)
            
            return {
                'name': compound_name,
                'cid': cid,
                'properties': properties,
                'image_response': image_response
            }
            
        except requests.RequestException as e:
            print(f"Network error: {str(e)}")
            return None
        except Exception as e:
            print(f"Error: {str(e)}")
            return None
    
    def display_results(self, result):
        """Display the search results"""
        if not result:
            return
        
        print("\n" + "="*50)
        print("CHEMICAL COMPOUND INFORMATION")
        print("="*50)
        print(f"Compound: {result['name']}")
        print(f"CID: {result['cid']}")
        print(f"Molecular Formula: {result['properties'].get('MolecularFormula', 'N/A')}")
        print(f"Molecular Weight: {result['properties'].get('MolecularWeight', 'N/A')} g/mol")
        print(f"IUPAC Name: {result['properties'].get('IUPACName', 'N/A')}")
        print(f"SMILES: {result['properties'].get('CanonicalSMILES', 'N/A')}")
        
        # Save image if available
        if result['image_response'].status_code == 200:
            try:
                image_data = result['image_response'].content
                image = Image.open(io.BytesIO(image_data))
                
                # Save image with compound name
                filename = f"{result['name'].replace(' ', '_')}_{result['cid']}.png"
                image.save(filename)
                print(f"\nImage saved as: {filename}")
                
                # Display image info
                print(f"Image size: {image.size[0]} x {image.size[1]} pixels")
                
            except Exception as e:
                print(f"Could not save image: {str(e)}")
        else:
            print("\nImage not available")
        
        print("="*50)
    
    def interactive_mode(self):
        """Run in interactive mode"""
        print("Chemical Compound Lookup Tool")
        print("Enter compound names to search (type 'quit' to exit)")
        print("-" * 40)
        
        while True:
            compound_name = input("\nEnter compound name: ").strip()
            
            if compound_name.lower() in ['quit', 'exit', 'q']:
                print("Goodbye!")
                break
            
            if not compound_name:
                print("Please enter a compound name.")
                continue
            
            result = self.search_compound(compound_name)
            self.display_results(result)

def main():
    """Main function"""
    if len(sys.argv) > 1:
        # Command line mode
        compound_name = ' '.join(sys.argv[1:])
        lookup = ChemicalCompoundLookupCLI()
        result = lookup.search_compound(compound_name)
        lookup.display_results(result)
    else:
        # Interactive mode
        lookup = ChemicalCompoundLookupCLI()
        lookup.interactive_mode()

if __name__ == "__main__":
    main() 