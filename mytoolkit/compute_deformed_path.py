#!/usr/bin/env python3

import sys
import itertools
import collections
import json

import pandas as pd
import numpy as np

from pymatgen.core.composition import Composition
from pymatgen.analysis.phase_diagram import PDEntry
from pymatgen.analysis.reaction_calculator import Reaction, ComputedReaction, ReactionError


def get_entry(products, reactants):
    products_entries = []
    for prod in products:
        cmp = Composition(prod.split("-")[-1])
        eng = 0
        entry = PDEntry(cmp, eng)
        entry.entry_id = prod.split("-")[0]
        entry.entry_name = prod
        # print( entry.entry_name ) ; input()
        products_entries.append(entry)
    reactants_entries = []
    for react in reactants:
        cmp = Composition(react.split("-")[-1])
        eng = 0
        entry = PDEntry(cmp, eng)
        entry.entry_id = react.split("-")[0]
        entry.entry_name = react
        # print( entry.entry_name ) ; input()
        reactants_entries.append(entry)
    return products_entries, reactants_entries

def compute_form_energy(react_compounds_energy_coeffs:list, prod_compounds_energy_coeffs:list):
    # product_idx 是生成物在product_reacants_entries中的索引号组成的列表
    form_energy = 0
    left_total_numatm = 0
    right_total_numatm = 0
    for comp, energy, coeff in react_compounds_energy_coeffs:
        natm = comp.num_atoms
        left_total_numatm += natm*np.abs(coeff)
        form_energy  += natm*np.abs(coeff)*energy
 
    for comp, energy, coeff in prod_compounds_energy_coeffs:
        natm = comp.num_atoms
        right_total_numatm += natm*np.abs(coeff)
        form_energy  += - natm*np.abs(coeff)*energy

    form_energy_peratom = None
    if   np.isclose(left_total_numatm, right_total_numatm, atol=1e-3) and left_total_numatm != 0:
        form_energy_peratom = 1000*form_energy/right_total_numatm
    elif np.isclose(left_total_numatm, right_total_numatm, atol=1e-3) and left_total_numatm == 0:
        form_energy_peratom = None
    else:
        print(f"The left part numberofatom:{left_total_numatm} != the left part numberofatom:{right_total_numatm}")
        sys.exit
    return form_energy_peratom


def check_repeatability(composition_entries:list, numberofreactions):
    compsition_set = set(composition_entries)
    elements = []
    for comp in compsition_set:
        ele = comp.composition.elements
        elements.extend(ele)
    if len(set(elements)) == numberofreactions:
        return True
    else:
        return False

def check_product_number(all_coefs):
    count_positive = 0

    for num in all_coefs:
        if num > 0:
            count_positive += 1
            if count_positive > 1:
                return False
    return True

def check_product_coeffs(products_entries):
    for prod in products_entries:
        product_coef = reaction.get_coeff(products_entries[0].composition)
        if np.isclose(product_coef, 0.0, atol=1e-5):
            return False
        else:
            return True


if __name__ == '__main__':

    enthalpy_csvfile =     sys.argv[1]
    numberofproducts = int(sys.argv[2])
    numberofreactions= int(sys.argv[3])
    enthalpy_datas   = pd.read_csv(enthalpy_csvfile, index_col=0, header=1, na_values=['--'])
    products   = enthalpy_datas.iloc[:,-numberofproducts:].columns
    reactants  = enthalpy_datas.iloc[:,:-numberofproducts].columns
    products_entries, reactants_entries = get_entry(products, reactants)
    reactions = []
    validnum = 0
    totalnum = 0
    for numb in range(2, numberofreactions+1):
        for reactant_entry in list(itertools.combinations(reactants_entries, numb)):
            reactant_entry = list(reactant_entry)
            totalnum += 1
            if check_repeatability(reactant_entry, numberofreactions):
                try:
                    reaction = ComputedReaction(reactant_entry, products_entries)
                    if check_product_coeffs(products_entries) and check_product_number(reaction.coeffs):
                        print(reaction.__str__())
                        reactions.append(reaction)
                        validnum += 1
                except ReactionError:
                    # print(f"Reaction cannot be balanced. {reactant_entry} {products_entries}")
                    pass
            else:
                # print(f"Reapted!!{reactant_entry}")
                pass
    print(f"{totalnum} attempts in total computing reaction")
    print(f"{validnum} attempts in valid computing reaction")
    # input("partially successful")
    result_dt = collections.defaultdict(list)
    presses = []
    for press, comp_series in enthalpy_datas.iterrows():
        presses.append(press)
        for reaction in reactions:
            react_compounds_energy_coeffs = []
            left_ids = []
            for comp, entry in zip(reaction.reactants, reaction._reactant_entries):
                name = entry.entry_name
                enid = entry.entry_id
                eg_peratoms = comp_series.loc[entry.entry_name]
                coefficient = reaction.get_coeff(comp)
                react_compounds_energy_coeffs.append([comp, eg_peratoms, coefficient])
                left_ids.append(enid)
            prod_compounds_energy_coeffs = []
            right_ids = []
            for comp, entry in zip(reaction.products, reaction._product_entries):
                name = entry.entry_name
                enid = entry.entry_id
                eg_peratoms = comp_series.loc[entry.entry_name]
                coefficient = reaction.get_coeff(comp)
                prod_compounds_energy_coeffs.append([comp, eg_peratoms, coefficient])
                right_ids.append(enid)
            
            form_energy_peratom = compute_form_energy(
                react_compounds_energy_coeffs, 
                prod_compounds_energy_coeffs,
                )
            
            left_string  = '+'.join(left_ids)
            right_string = '+'.join(right_ids)
            total_string = '[' + left_string + '->' + right_string + ']' + reaction.__str__()
            if form_energy_peratom is not None:
                # result_dt[reaction.__str__()].append(form_energy_peratom)
                result_dt[total_string].append(form_energy_peratom)

    result_json = json.dumps(result_dt)
    result_file = open("formed-enthalpy.json", 'w')
    result_file.write(result_json)
    result_file.close()

    total_result_pd = pd.DataFrame(data=result_dt, index=presses)
    total_result_pd.to_csv("formed-enthalpy.csv")
    