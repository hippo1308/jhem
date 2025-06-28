# Chemical Compound Lookup Tool

This tool allows you to search for chemical compounds by name and retrieve their images and information from PubChem.

## Features

- **Name-based search**: Search compounds by their common names, IUPAC names, or synonyms
- **Fuzzy matching**: If an exact match isn't found, the tool will try to find the best match
- **Rich information**: Get molecular formula, molecular weight, IUPAC name, and SMILES notation
- **Chemical structure images**: View and save 2D chemical structure images
- **Two interfaces**: GUI (tkinter) and command-line versions available

## Installation

1. Make sure you have Python 3.6+ installed
2. Install the required dependencies:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

### GUI Version (Recommended)

Run the graphical interface:
```bash
python xtra.py
```

This will open a window where you can:
1. Enter a compound name in the search box
2. Click "Search" or press Enter
3. View the compound information and image
4. The image will be displayed in the application window

### Command-Line Version

#### Interactive Mode
```bash
python chemical_lookup_cli.py
```

This will start an interactive session where you can enter compound names one by one.

#### Single Search
```bash
python chemical_lookup_cli.py "aspirin"
python chemical_lookup_cli.py "ethanol"
python chemical_lookup_cli.py "caffeine"
```

## Examples

Try searching for these common compounds:
- **Aspirin** (acetylsalicylic acid)
- **Caffeine** (1,3,7-trimethylxanthine)
- **Ethanol** (ethyl alcohol)
- **Glucose** (dextrose)
- **Methane** (CH4)
- **Water** (H2O)

## How It Works

1. **Search**: The tool sends your compound name to PubChem's REST API
2. **Match**: It finds the best matching compound (exact or fuzzy match)
3. **Retrieve**: Gets detailed information including molecular properties
4. **Display**: Shows the chemical structure image and information

## API Information

This tool uses PubChem's PUG REST API:
- **Base URL**: https://pubchem.ncbi.nlm.nih.gov/rest/pug
- **No API key required**: PubChem provides free access to their database
- **Rate limiting**: Please be respectful and don't make too many requests too quickly

## Output Information

For each compound, you'll get:
- **CID**: PubChem Compound Identifier
- **Molecular Formula**: Chemical formula (e.g., C9H8O4 for aspirin)
- **Molecular Weight**: Mass in g/mol
- **IUPAC Name**: Official chemical name
- **SMILES**: Simplified molecular input line entry system notation
- **2D Structure**: Chemical structure image

## Error Handling

The tool handles various error cases:
- **Compound not found**: Suggests trying a different name
- **Network errors**: Provides clear error messages
- **Image loading issues**: Gracefully handles missing images

## Tips for Better Searches

1. **Use common names**: "aspirin" works better than "acetylsalicylic acid"
2. **Try synonyms**: "ethanol" or "ethyl alcohol" both work
3. **Be specific**: "caffeine" is better than "coffee"
4. **Check spelling**: Chemical names can be complex

## Files

- `xtra.py` - GUI version of the tool
- `chemical_lookup_cli.py` - Command-line version
- `requirements.txt` - Python dependencies
- `README_chemical.md` - This documentation

## Dependencies

- `requests` - For HTTP requests to PubChem API
- `Pillow` - For image processing
- `tkinter` - For GUI (included with Python)

## License

This tool is provided as-is for educational and research purposes. Please respect PubChem's terms of service when using their API. 