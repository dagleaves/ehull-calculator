# Energy Above the Convex Hull Calculator

Interfaces with the Materials Project database to calculate the energy above the convex hull (eV/atom) for each material in a given csv
containing material formulas and their corresponding total energies per atom (must be per atom, not total overall).

## Dependencies

There are three dependencies:

```
mp-api==0.30.8
mpcontribs-client==5.1.1
pymatgen==2023.1.30
```

A requirements.txt file is also provided.

```
pip install -r requirements.txt
```

## Input file format

The input file should be formatted as follows:

```
formula,total_energy_per_atom
La9ErO15,-1.323633793
```

## Running the program

```
python ehull.py -i input.csv -o output.csv
```

## Expected output
```
formula,ev_atom
La9ErO15,7.497252317
```
