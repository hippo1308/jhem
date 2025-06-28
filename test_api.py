import requests
import json

def test_pubchem_api():
    base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    compound_name = "aspirin"
    
    print(f"Testing PubChem API for: {compound_name}")
    
    # Test the search URL
    search_url = f"{base_url}/compound/name/{compound_name}/cids/JSON"
    print(f"Search URL: {search_url}")
    
    try:
        response = requests.get(search_url, timeout=10)
        print(f"Status code: {response.status_code}")
        
        if response.status_code == 200:
            data = response.json()
            print(f"Response data: {json.dumps(data, indent=2)}")
            
            # Try to extract CIDs
            cids = data.get('InformationList', {}).get('Information', [{}])[0].get('CID', [])
            print(f"Found CIDs: {cids}")
            
            if cids:
                cid = cids[0]
                print(f"Using CID: {cid}")
                
                # Test getting compound info
                info_url = f"{base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
                info_response = requests.get(info_url, timeout=10)
                print(f"Info status code: {info_response.status_code}")
                
                if info_response.status_code == 200:
                    info_data = info_response.json()
                    print(f"Info data: {json.dumps(info_data, indent=2)}")
                    
                    # Test getting image
                    image_url = f"{base_url}/compound/cid/{cid}/PNG"
                    image_response = requests.get(image_url, timeout=10)
                    print(f"Image status code: {image_response.status_code}")
                    print(f"Image content length: {len(image_response.content)}")
                else:
                    print(f"Info request failed: {info_response.text}")
            else:
                print("No CIDs found in response")
        else:
            print(f"Search request failed: {response.text}")
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_pubchem_api() 