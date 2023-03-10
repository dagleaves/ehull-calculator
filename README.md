# Energy Above the Convex Hull Calculator

Interfaces with the Materials Project database to calculate the e_hull for each material in a given csv
containing material formulas and their corresponding total energies.

## Dependencies

There are two dependencies:

```
mp_api==0.30.8
pymatgen==2023.1.30
```

A requirements.txt file is also provided.

```
pip install -r requirements.txt
```

## Input file format

The input file should be formatted as follows:

```
formula,total_energy
Eu2HgO6,-6.666567894
```

## Running the program

```
python ehull.py -i input.csv -o output.csv
```


