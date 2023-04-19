import os
import re
import csv
import sys
import time
import argparse
from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.composition import Composition
from mp_api.client import MPRester

API_KEY = "MP-API-Key"


def formulaExpander(formula):
    while len(re.findall('\(\w*\)', formula)) > 0:
        parenthetical = re.findall('\(\w*\)[0-9]+', formula)
        for i in parenthetical:
            p = re.findall('[0-9]+', str(re.findall('\)[0-9]+', i)))
            j = re.findall('[A-Z][a-z]*[0-9]*', i)
            for n in range(0, len(j)):
                numero = re.findall('[0-9]+', j[n])
                if len(numero) != 0:
                    for k in numero:
                        nu = re.sub(k, str(int(int(k) * int(p[0]))), j[n])
                else:
                    nu = re.sub(j[n], j[n] + p[0], j[n])
                j[n] = nu
        newphrase = ""
        for m in j:
            newphrase += str(m)
        formula = formula.replace(i, newphrase)
        if (len((re.findall('\(\w*\)[0-9]+', formula))) == 0) and (len(re.findall('\(\w*\)', formula)) != 0):
            formula = formula.replace('(', '')
            formula = formula.replace(')', '')
    return formula


def get_num_atoms(formula):
    s = re.findall('([A-Z][a-z]?)([0-9]?\.?[0-9]*)', formula)
    num_atoms = 0
    for sa, sb in s:
        if sb == '':
            num_atoms += 1
        else:
            num_atoms += int(sb)
    return num_atoms


def powerset(s):
    ret = []
    x = len(s)
    for i in range(1 << x):
        ret.append([s[j] for j in range(x) if (i & (1 << j))])
    return ret


def get_E_above_hull(formula, E_atom):
    formula = formulaExpander(formula)

    num_atoms = get_num_atoms(formula)
    E_total = E_atom * num_atoms

    comp = Composition(formula)
    entry = ComputedEntry(composition=comp, energy=E_total)

    elements = comp.elements
    elements = [elem.symbol for elem in elements]

    ps = ['-'.join(e) for e in powerset(elements)[1:]]
    mpr = MPRester(API_KEY)

    competitive_compositions = mpr.summary.search(
        chemsys=ps, fields=['composition', 'energy_per_atom'])

    competitive_entries = []
    for indx1, comp in enumerate(competitive_compositions):
        competitive_comp = comp.composition
        if str(competitive_comp).replace(' ', '') == formula:
            continue

        num_atoms = get_num_atoms(str(competitive_comp).replace(' ', ''))
        competitive_energy = comp.energy_per_atom * num_atoms

        competitive_entry = ComputedEntry(

            composition=competitive_comp, energy=competitive_energy)
        competitive_entries.append(competitive_entry)
    all_entries = [entry] + competitive_entries

    phase_diagram = PhaseDiagram(all_entries)
    decomp, energy_above_hull = phase_diagram.get_decomp_and_e_above_hull(
        entry)

    return (decomp, energy_above_hull)


def run_debug():
    print('Running in DEBUG mode')
    mpr = MPRester(API_KEY)
    print('Collecting materials from MP database')
    mats = mpr.summary.search(
        # material_ids=['mp-149', 'mp-505569', 'mp-936295', 'mp-862690'], fields=['composition', 'energy_per_atom', 'energy_above_hull', 'material_id'])
        material_ids=['mp-785008'], fields=['composition', 'energy_per_atom', 'energy_above_hull', 'material_id'])

    print('ID', 'Real E_hull', 'Calc E_hull')
    for mat in mats:
        print(mat.energy_per_atom)
        print(mat.material_id, mat.energy_above_hull,
              get_E_above_hull(str(mat.composition).replace(' ', ''), mat.energy_per_atom)[1])


def main(args):
    save = []
    with open(args.input, 'r') as f:
        reader = csv.reader(f)
        next(reader)

        for row in reader:
            if '-' in row:
                formula = row[0].split('-')[1]
            else:
                formula = row[0]
            print(formula)
            E_atom = float(row[1])

            decomp, energy_above_hull = get_E_above_hull(formula, E_atom)
            save.append([row[0], energy_above_hull])
            time.sleep(1)

    with open(args.output, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerows(save)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Calculates energy above the convex hull for a csv of formulas and total energies (with header)')
    parser.add_argument('--debug', action='store_true',
                        help='Enable debug mode (test imports)')
    parser.add_argument('-i', '--input', type=str,
                        required=True, help='Input csv filename')
    parser.add_argument('-o', '--output', type=str,
                        required=True, help='Output csv filename')
    args = parser.parse_args(sys.argv[1:])

    assert os.path.exists(args.input), 'Input file does not exist'
    if args.debug:
        run_debug()
    else:
        main(args)
