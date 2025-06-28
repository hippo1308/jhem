import requests
from PIL import Image
import io

def test_image_download():
    print("Testing image download...")
    
    # Test with aspirin CID
    cid = 2244
    image_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
    
    try:
        response = requests.get(image_url, timeout=10)
        print(f"Image response status: {response.status_code}")
        print(f"Image content length: {len(response.content)}")
        
        if response.status_code == 200:
            # Try to open the image
            image = Image.open(io.BytesIO(response.content))
            print(f"Image opened successfully: {image.size}")
            
            # Try to save the image
            filename = "test_aspirin.png"
            image.save(filename)
            print(f"Image saved as: {filename}")
            
            # Check if file exists
            import os
            if os.path.exists(filename):
                print(f"File exists: {filename}")
                print(f"File size: {os.path.getsize(filename)} bytes")
            else:
                print("File was not created")
                
        else:
            print(f"Failed to get image: {response.text}")
            
    except Exception as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    test_image_download() 