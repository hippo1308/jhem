import requests
import json
from PIL import Image, ImageTk
import io
import tkinter as tk
from tkinter import ttk, messagebox
import threading
import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import webbrowser

class ChemicalCompoundLookup:
    def __init__(self):
        self.base_url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
        self.setup_gui()
    
    def setup_gui(self):
        """Setup the graphical user interface"""
        self.root = tk.Tk()
        self.root.title("Chemical Compound Lookup")
        self.root.geometry("600x500")
        self.root.configure(bg='#f0f0f0')
        
        self.notebook = ttk.Notebook(self.root)
        self.notebook.pack(fill='both', expand=True)

        # --- Compound Name Tab ---
        self.compound_frame = ttk.Frame(self.notebook, padding="20")
        self.notebook.add(self.compound_frame, text="Name to Structure")
        self.setup_name_tab(self.compound_frame)

        # --- SMILES to IUPAC Tab ---
        self.smiles_frame = ttk.Frame(self.notebook, padding="20")
        self.notebook.add(self.smiles_frame, text="SMILES to IUPAC")
        self.setup_smiles_tab(self.smiles_frame)
    
    def setup_name_tab(self, main_frame):
        # Title
        title_label = ttk.Label(main_frame, text="Chemical Compound Lookup", 
                               font=('Arial', 16, 'bold'))
        title_label.grid(row=0, column=0, columnspan=2, pady=(0, 20))
        
        # Search frame
        search_frame = ttk.Frame(main_frame)
        search_frame.grid(row=1, column=0, columnspan=2, sticky="ew", pady=(0, 20))
        
        # Search label and entry
        ttk.Label(search_frame, text="Compound Name:").grid(row=0, column=0, sticky="w", padx=(0, 10))
        self.search_var = tk.StringVar()
        self.search_entry = ttk.Entry(search_frame, textvariable=self.search_var, width=30)
        self.search_entry.grid(row=0, column=1, sticky="ew", padx=(0, 10))
        self.search_entry.bind('<Return>', lambda e: self.search_compound())
        
        # Search button
        self.search_button = ttk.Button(search_frame, text="Search", command=self.search_compound)
        self.search_button.grid(row=0, column=2)
        
        # Progress bar
        self.progress = ttk.Progressbar(main_frame, mode='indeterminate')
        self.progress.grid(row=2, column=0, columnspan=2, sticky="ew", pady=(0, 20))
        
        # Results frame
        results_frame = ttk.LabelFrame(main_frame, text="Results", padding="10")
        results_frame.grid(row=3, column=0, columnspan=2, sticky="nsew", pady=(0, 20))
        
        # Compound info
        self.info_text = tk.Text(results_frame, height=8, width=60, wrap=tk.WORD)
        self.info_text.grid(row=0, column=0, columnspan=2, sticky="ew", pady=(0, 10))
        
        # Image label
        self.image_label = ttk.Label(results_frame, text="No image available")
        self.image_label.grid(row=1, column=0, columnspan=2)
        
        # Configure grid weights
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(3, weight=1)
        search_frame.columnconfigure(1, weight=1)
        results_frame.columnconfigure(0, weight=1)
        results_frame.rowconfigure(0, weight=1)
    
    def setup_smiles_tab(self, frame):
        title_label = ttk.Label(frame, text="SMILES to IUPAC Name", font=('Arial', 16, 'bold'))
        title_label.grid(row=0, column=0, columnspan=2, pady=(0, 20))
        ttk.Label(frame, text="Paste SMILES string:").grid(row=1, column=0, sticky="w", padx=(0, 10))
        self.smiles_var = tk.StringVar()
        self.smiles_entry = ttk.Entry(frame, textvariable=self.smiles_var, width=40)
        self.smiles_entry.grid(row=1, column=1, sticky="ew", padx=(0, 10))
        self.smiles_entry.bind('<Return>', lambda e: self.smiles_to_iupac())
        self.smiles_button = ttk.Button(frame, text="Get IUPAC Name", command=self.smiles_to_iupac)
        self.smiles_button.grid(row=1, column=2)
        self.drawer_button = ttk.Button(frame, text="Open Structure Drawer (Browser)", command=self.open_ketcher)
        self.drawer_button.grid(row=2, column=1, pady=(10, 5))
        self.iupac_result = tk.Text(frame, height=4, width=60, wrap=tk.WORD)
        self.iupac_result.grid(row=3, column=0, columnspan=3, pady=(20, 0))
        frame.columnconfigure(1, weight=1)

    def search_compound(self):
        """Search for a compound by name"""
        compound_name = self.search_var.get().strip()
        if not compound_name:
            messagebox.showwarning("Warning", "Please enter a compound name")
            return
        
        # Start search in a separate thread to avoid blocking GUI
        self.search_button.config(state='disabled')
        self.progress.start()
        self.info_text.delete(1.0, tk.END)
        self.image_label.config(text="Searching...", image='')
        
        thread = threading.Thread(target=self._perform_search, args=(compound_name,))
        thread.daemon = True
        thread.start()
    
    def _perform_search(self, compound_name):
        """Perform the actual search in a separate thread"""
        try:
            # First, search for the compound by name
            search_url = f"{self.base_url}/compound/name/{compound_name}/cids/JSON"
            response = requests.get(search_url, timeout=10)
            
            if response.status_code != 200:
                # Try fuzzy search
                search_url = f"{self.base_url}/compound/name/{compound_name}/cids/JSON?name_type=word"
                response = requests.get(search_url, timeout=10)
                
                if response.status_code != 200:
                    self.root.after(0, lambda: self._show_error("Compound not found. Please try a different name."))
                    return
            
            data = response.json()
            
            # Try both old and new API response formats
            cids = []
            if 'InformationList' in data:
                cids = data.get('InformationList', {}).get('Information', [{}])[0].get('CID', [])
            elif 'IdentifierList' in data:
                cids = data.get('IdentifierList', {}).get('CID', [])
            
            if not cids:
                self.root.after(0, lambda: self._show_error("No compounds found with that name."))
                return
            
            # Get the first (best match) CID
            cid = cids[0]
            
            # Get compound information
            info_url = f"{self.base_url}/compound/cid/{cid}/property/MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES/JSON"
            info_response = requests.get(info_url, timeout=10)
            
            if info_response.status_code != 200:
                self.root.after(0, lambda: self._show_error("Could not retrieve compound information."))
                return
            
            info_data = info_response.json()
            print('DEBUG: info_data =', info_data)  # Debug print
            # Handle both old and new formats
            properties = {}
            if 'PropertyTable' in info_data and 'Properties' in info_data['PropertyTable']:
                props_list = info_data['PropertyTable']['Properties']
                if props_list and isinstance(props_list, list):
                    properties = props_list[0]
            properties.setdefault('MolecularFormula', 'N/A')
            properties.setdefault('MolecularWeight', 'N/A')
            properties.setdefault('IUPACName', 'N/A')
            properties.setdefault('CanonicalSMILES', 'N/A')

            # Fetch synonyms to check for a real match
            synonyms_url = f"{self.base_url}/compound/cid/{cid}/synonyms/JSON"
            synonyms_response = requests.get(synonyms_url, timeout=10)
            synonyms_data = synonyms_response.json() if synonyms_response.status_code == 200 else {}
            synonyms = []
            if 'InformationList' in synonyms_data:
                info = synonyms_data['InformationList'].get('Information', [{}])[0]
                synonyms = info.get('Synonym', [])
            # Prepare all possible names to check
            all_names = [str(properties.get('IUPACName', '')).strip().lower()]
            all_names += [s.strip().lower() for s in synonyms]
            user_input = compound_name.strip().lower()
            # Only accept if the user input matches exactly any synonym or IUPAC name
            if user_input not in all_names:
                self.root.after(0, lambda: self._show_error("No relevant compound found for your input. Please try a different name."))
                return

            # Get compound image
            image_url = f"{self.base_url}/compound/cid/{cid}/PNG"
            image_response = requests.get(image_url, timeout=10)
            
            # Update GUI with results
            self.root.after(0, lambda: self._display_results(compound_name, cid, properties, image_response))
            
        except requests.RequestException as e:
            self.root.after(0, lambda: self._show_error(f"Network error: {str(e)}"))
        except Exception as e:
            self.root.after(0, lambda: self._show_error(f"Error: {str(e)}"))
        finally:
            self.root.after(0, self._finish_search)
    
    def _display_results(self, compound_name, cid, properties, image_response):
        """Display the search results in the GUI"""
        # Display compound information
        info_text = f"Compound: {compound_name}\n"
        info_text += f"CID: {cid}\n"
        info_text += f"Molecular Formula: {properties.get('MolecularFormula', 'N/A')}\n"
        info_text += f"Molecular Weight: {properties.get('MolecularWeight', 'N/A')} g/mol\n"
        info_text += f"IUPAC Name: {properties.get('IUPACName', 'N/A')}\n"
        info_text += f"SMILES: {properties.get('CanonicalSMILES', 'N/A')}\n"
        
        self.info_text.delete(1.0, tk.END)
        self.info_text.insert(1.0, info_text)
        
        # Display image if available
        if image_response.status_code == 200:
            try:
                image_data = image_response.content
                image = Image.open(io.BytesIO(image_data))
                
                # Resize image to fit in GUI (max 300x300)
                max_size = (300, 300)
                image.thumbnail(max_size, Image.Resampling.LANCZOS)
                
                # Convert to PhotoImage for tkinter
                photo = ImageTk.PhotoImage(image)
                self.image_label.config(image=photo, text="")
                # Store reference to prevent garbage collection
                self.current_image = photo
            except Exception as e:
                self.image_label.config(text=f"Could not load image: {str(e)}", image='')
        else:
            self.image_label.config(text="Image not available", image='')
    
    def _show_error(self, message):
        """Show error message"""
        messagebox.showerror("Error", message)
        self.info_text.delete(1.0, tk.END)
        self.image_label.config(text="No image available", image='')
    
    def _finish_search(self):
        """Finish the search process"""
        self.progress.stop()
        self.search_button.config(state='normal')
    
    def smiles_to_iupac(self):
        smiles = self.smiles_var.get().strip()
        self.iupac_result.delete(1.0, tk.END)
        if not smiles:
            self.iupac_result.insert(tk.END, "Please enter a SMILES string.")
            return
        # Try RDKit first
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError("Invalid SMILES string.")
            iupac = rdMolDescriptors.CalcMolFormula(mol)
            # Try to get IUPAC name from PubChem as RDKit does not provide it directly
            iupac_name = self.get_iupac_from_pubchem(smiles)
            if iupac_name:
                self.iupac_result.insert(tk.END, f"IUPAC Name: {iupac_name}\n")
            self.iupac_result.insert(tk.END, f"Molecular Formula: {iupac}\n")
        except Exception as e:
            self.iupac_result.insert(tk.END, f"Error: {str(e)}\n")

    def get_iupac_from_pubchem(self, smiles):
        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/property/IUPACName/JSON"
            r = requests.get(url, timeout=10)
            if r.status_code == 200:
                data = r.json()
                props = data.get('PropertyTable', {}).get('Properties', [{}])[0]
                return props.get('IUPACName', None)
        except Exception:
            pass
        return None

    def open_ketcher(self):
        # Open the Ketcher Single Application Demo in the default browser
        webbrowser.open_new_tab('https://lifescience.opensource.epam.com/KetcherDemoSA/index.html')

    def run(self):
        """Start the GUI application"""
        self.root.mainloop()

def main():
    """Main function to run the chemical compound lookup tool"""
    app = ChemicalCompoundLookup()
    app.run()

if __name__ == "__main__":
    main()