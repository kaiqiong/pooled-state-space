CGMacros v1.0.0

Source: https://physionet.org/content/cgmacros/1.0.0/
DOI: https://doi.org/10.13026/3z8q-x658
The dataset is distributed by PhysioNet under CC BY-NC-SA 4.0; see the source page for its terms and citation. This code repository does not include participant records or meal photographs.

Run from the repository root:
  python3 scripts/download_cgmacros.py

Or use an existing copy of the official archive:
  python3 scripts/download_cgmacros.py --archive /path/to/CGMacros_dateshifted365.zip

The script checks the ZIP against the official v1.0.0 SHA-256 manifest and extracts only bio.csv and the 45 participant CSVs. About 627 MB must be downloaded because the source bundles recordings with photographs; photographs are not extracted.
Official checksum source:
  https://physionet.org/files/cgmacros/1.0.0/SHA256SUMS.txt
SHA-256 of CGMacros_dateshifted365.zip:
  05c8b0e6f1a2757050aced55ce4bf6ab2ac9b30f2fd8ca193056812d9c621d4d

Expected layout:
  data/raw/CGMacros/bio.csv
  data/raw/CGMacros/CGMacros-NNN/CGMacros-NNN.csv

Alternatively point the application runner at an existing extracted folder using --data DIR. DIR must contain bio.csv and the participant subdirectories. Preserve source files unchanged; all generated windows, fitted objects, and bootstrap draws are written under results/, which is ignored by Git.

Dataset citation:
Gutierrez-Osuna R, Kerr D, Mortazavi B, Das A (2025). CGMacros: a scientific dataset for personalized nutrition and diet monitoring (version 1.0.0). PhysioNet. doi:10.13026/3z8q-x658.
