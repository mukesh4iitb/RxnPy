import re
import os
import sys
import json
import numpy as np
import matplotlib.pyplot as plt
import pyexcel_ods3 as ods
from pyexcel_ods3 import get_data, save_data
from collections import OrderedDict



def en_component(component, energy_dict):
    """
    Note in reaction.in, it can handle the number factor but multiplication should be denoted by ** and division by /.
    Also, it can handle parenthesis around numerical factor.
    """
    if "**" in component:
        factor_species = component.split("**")
        factor_str = factor_species[0].strip()

        if "/" in factor_str:
            factor_str = factor_str.strip("()")
            num, denom = map(float, factor_str.split('/'))
            factor = num / denom
        else:
            factor = float(factor_str)
        return factor*float(energy_dict[factor_species[1].strip()])
    else:
        factor = 1
        return factor*float(energy_dict[component.strip()])



#def en_component_requirement(component, energy_dict):
#    energy_dict_name = [name for name, value in globals().items() if value is energy_dict][0]
#    if "**" in component:
#        factor_species = component.split("**")
#        species = factor_species[1].strip()
#        if species not in energy_dict.keys():
#            return "{}[{}]".format(energy_dict_name, species)
#    else:
#        species = component.strip()
#        if species not in energy_dict.keys():
#            return "{}[{}]".format(energy_dict_name, species)


def Rxn_energy(reaction, energy_dict):
    
    half_reactions = reaction.split('=')
    reactants = half_reactions[0].split('+')
    products = half_reactions[1].split('+')

    en_reactants = 0
    for component in reactants:
        en = en_component(component, energy_dict) #float(energy_dict[component.strip()])
        en_reactants += en
        
    en_products = 0
    for component in products:
        en = en_component(component, energy_dict)
        #en = float(energy_dict[component.strip()])
        en_products += en
    return en_products - en_reactants


#def Rxn_energy_requirement(reaction, energy_dict):
#
#    half_reactions = reaction.split('=')
#    reactants = half_reactions[0].split('+')
#    products = half_reactions[1].split('+')
#
#    requirement = []
#    for component in reactants:
#        req = en_component_requirement(component, energy_dict) #float(energy_dict[component.strip()])
#        if req:
#            requirement.append(req)
#
#    for component in products:
#        req = en_component_requirement(component, energy_dict)
#        if req:
#            requirement.append(req)
#
#    return requirement


def E0_reader(filename):
    with open(filename) as reader:
        energy=reader.read()
        energy=re.sub(r"\[.*?\]", "", energy).strip()  # remove anything inside [], which maybe useful while keeping the full name of the file in E0 etc.
        energy=energy.splitlines()
    
    energy_dict={}
    for en in energy:
        
        en = en.strip()
    
        if not en or en.startswith("#"):  # Skip empty or commented lines
            continue  
        
        part=en.split()
        energy_dict[part[-1]]=float(part[-2])
    
    return energy_dict

def E0_new_reader(filename):
    with open(filename) as reader:
        energy=reader.read()
        energy=re.sub(r"\[.*?\]", "", energy).strip()  # remove anything inside [], which maybe useful while keeping the full name of the file in E0 etc.
        energy=energy.splitlines()
    #print(energy)

    energy_dict={}
    for en in energy:

        en = en.strip()

        if not en or en.startswith("#"):  # Skip empty or commented lines
            continue
        #print("en:", en)
        #print(en.split("Energy:")[1].split()[0])
        #print(en.split("Name:")[1].split()[0])
        energy_str = en.split("Energy:")[1].split()[0]
        energy = float(energy_str)

        name = en.split("Name:")[1].split()[0]
        energy_dict[name]=float(energy)

    with open("results.jsonl", "a") as f:
        json.dump({"E0": energy_dict}, f)
        f.write("\n")


    return energy_dict


def thermal_correction_reader(filename):
    with open(filename) as reader:
        en_corr = reader.read()
        en_corr = re.sub(r"\[.*?\]", "", en_corr).strip()   # remove anything inside [], which maybe useful while keeping the full name of the file in E0 etc.
        en_corr = en_corr.splitlines()

    thermal_corrections_dict = {}
    for corr in en_corr:
        
        corr = corr.strip()
        if not corr or corr.startswith("#"):  # Skip empty or commented lines
            continue  

        parts = corr.split()
        sys = parts[-1]
        zpe, del_u, ts, g_corr = map(float, parts[1:-1]) 
        
        thermal_corrections_dict[sys] = {
            "ZPE": zpe,
            "Del_U": del_u,
            "TS": ts,
            "G_corr": g_corr}
    return thermal_corrections_dict

def thermal_correction_new_reader(filename):
    with open(filename) as reader:
        en_corr = reader.read()
        en_corr = re.sub(r"\[.*?\]", "", en_corr).strip()   # remove anything inside [], which maybe useful while keeping the full name of the file in E0 etc.
        en_corr = en_corr.splitlines()

    thermal_corrections_dict = {}
    for corr in en_corr:

        corr = corr.strip()
        if not corr or corr.startswith("#"):  # Skip empty or commented lines
            continue

        name = corr.split("Name:")[1].split()[0]
        zpe_str = corr.split("ZPE:")[1].split()[0]
        #print("zpe_str:", zpe_str)
        zpe = float(zpe_str)
        delu_str = corr.split("DelU:")[1].split()[0]
        delu = float(delu_str)
        print(corr.split("TdelS:"))
        tdels_str = corr.split("TdelS:")[1].split()[0]
        tdels = float(tdels_str)
        delg_str = corr.split("DelG:")[1].split()[0]
        delg = float(delg_str)

        thermal_corrections_dict[name] = {
            "ZPE": zpe,
            "Del_U": delu,
            "TS": tdels,
            "G_corr": delg}
    return thermal_corrections_dict

def check_Gibbs_energy():
    required_files = ["Reaction1.in", "E0.dat", "Thermal_correction.dat"]

    # Checking presence of files
    missing_files = [f for f in required_files if not os.path.exists(f)]

    gibbs_checklist = []

    if missing_files:
        for f in missing_files:
            print(f"Error: '{f}' is not present.")
        sys.exit(1)  # Exit the program with error code 1
    print("Gibbs energy calc required files:")
    for n, file in enumerate(required_files):
        print("{}-{}".format(n+1, file))
    print("are present.")
    print("-"*80)


    Gibbs = Gibbs_energy("Reaction1.in", "E0.dat", "Thermal_correction.dat")
    print(Gibbs)
    print("Gibbs energy of species:\n", json.dumps(Gibbs[0], indent=4, sort_keys=True))

    if len(Gibbs[1])==0:
        print("You can calculate all Reaction Gibbs free energy!")
    else:
        print("You can't calculate all Reaction Gibbs free energy because of absence of:\n")
        print("-"*80)
        print("Missing Gibbs energy of species because of absence of:\n", json.dumps(Gibbs[1], indent=4, sort_keys=True))

        for comp in Gibbs[1]:
            gibbs_checklist0="{}".format(comp)
            #print(gibbs_checklist0)
            gibbs_checklist.append(gibbs_checklist0)
        write_or_append_checklist("calculation_checklist.ods", "Gibbs", checklist=gibbs_checklist)


def Gibbs_energy(reaction_filename, E0_filename, thermal_corrections_filename):
    reaction_species = reaction_components(reaction_filename)
    E0 = E0_new_reader(E0_filename)
    thermal_corrections = thermal_correction_new_reader(thermal_corrections_filename)
    print(thermal_corrections)

    Gibbs = {}
    sorted_Gibbs = {}
    Gibbs_missing = []
    #sorted_Gibbs_missing = {}
    for rsp in reaction_species:
        missing = False  
        
        if rsp not in thermal_corrections:
            Gibbs_missing.append(f"Thermal_corrections['{rsp}']")
            missing = True  
            
        if rsp not in E0:
            Gibbs_missing.append(f"E0['{rsp}']")
            missing = True 
        
        if missing:
            continue  # Skip only if any key is missing

        # Only execute if both values are available
        G_corr = thermal_corrections[rsp]['ZPE'] -  thermal_corrections[rsp]['TS']
        # Import Note, you want to include delU correction, we can directly read from delG as below:
        # G_corr = thermal_corrections[rsp]['G_corr']
        E0_value = E0[rsp]
        Gibbs[rsp] = G_corr + E0_value

        sorted_Gibbs = {k: Gibbs[k] for k in sorted(Gibbs)}
        #sorted_Gibbs_missing = sorted(set(Gibbs_missing))
    
    return sorted_Gibbs, sorted(set(Gibbs_missing))

def check_absorption_energy():
    required_files = ["Reaction.in", "E0.dat"]

    # Checking presence of files
    missing_files = [f for f in required_files if not os.path.exists(f)]

    abs_checklist = []

    if missing_files:
        for f in missing_files:
            print(f"Error: '{f}' is not present.")
        sys.exit(1)  # Exit the program with error code 1
    print("Absorption energy calc required files:")
    for n, file in enumerate(required_files):
        print("{}-{}".format(n+1, file))
    print("are present.")
    print("-"*80)
    Abs_en = absorption_energy("Reaction.in", "E0.dat")
    for i in range(len(Abs_en[0])):
        print(Abs_en[0][i], Abs_en[1][i])
    print("-"*80)
    print("Missing E0 of the component(s):\n")
    print(Abs_en[2])
    print("-"*80)
    # for writing the abs_checklist:
    for comp in sorted(set(Abs_en[2])):
        abs_checklist0="E0[{}]".format(comp)
        print(abs_checklist0)
        abs_checklist.append(abs_checklist0)
    write_or_append_checklist("calculation_checklist.ods", "Absorption", checklist=abs_checklist)




def absorption_energy(reaction_filename, E0_filename):
    E0 = E0_new_reader(E0_filename)
    calc_reactions = []
    absorption_energy = []
    E0_missing = []
    reactions = reaction_reader(reaction_filename)
    E0 = E0_new_reader(E0_filename)

    for reaction in reactions:
        reaction_species = reaction_components0(reaction)
        calc_absorption_energy = True
        for rsp in reaction_species:
            if rsp not in E0:
                E0_missing.append(rsp)
                calc_absorption_energy = False
        
        if calc_absorption_energy:
            calc_reactions.append(reaction)
            absorption_energy.append(Rxn_energy(reaction, E0))
    return calc_reactions, absorption_energy, set(E0_missing)

def Rxn_Gibbs_energy(reaction_filename, E0_filename, thermal_corrections_filename):
    
    delG_reactions = []
    reactions = reaction_reader(reaction_filename)
    gibbs_dict, _ = Gibbs_energy(reaction_filename, E0_filename, thermal_corrections_filename)
    
    for reaction in reactions:
        delG_reactions.append(Rxn_energy(reaction, gibbs_dict))

    return reactions, delG_reactions



def reaction_reader(filename):
    reactions_list = []
    with open(filename) as reactions:
        reactions = reactions.readlines()


    for reaction in reactions:
        reaction = reaction.strip()  # Remove leading/trailing spaces

        if not reaction or reaction.startswith("#"):  # Skip empty or commented lines
            continue  # Move to next iteration
        reactions_list.append(reaction)
    return reactions_list

def reaction_components0(reaction):
    reaction_species0 = []
    half_reactions = reaction.split('=')
    for half_reaction in half_reactions:
        comps = half_reaction.strip().split("+")
        for comp in comps:
            #print(comp)
            if "**" in comp:
                reaction_species0.append(comp.split("**")[1].strip())
            else:
                reaction_species0.append(comp.strip())
    return reaction_species0


def reaction_components(filename):
    reaction_species = []
    reactions = reaction_reader(filename)

    for reaction in reactions:
        half_reactions = reaction.split('=')
        for half_reaction in half_reactions:
            comps = half_reaction.strip().split("+")
            for comp in comps:
                #print(comp)
                if "**" in comp:
                    reaction_species.append(comp.split("**")[1].strip())
                else:
                    reaction_species.append(comp.strip())
    return reaction_species

def Calc_absorption_en(reaction_energy_list=False):
    print()
    print("\033[1;32mAbsorption energy:\033[0m")
    if os.path.exists("Reaction.in"):
        if os.path.exists("E0.dat"):
            Abs_en = absorption_energy("Reaction.in", "E0.dat")
            if Abs_en[2]:
                check_absorption_energy()
            else:
                for i in range(len(Abs_en[0])):
                    print("--"*40)
                    react = "" 
                    re_parts = Abs_en[0][i].split("=")
                    re_parts = Abs_en[0][i].split("=")
                    left_side = "  +  ".join(part.strip() for part in re_parts[0].split("+"))
                    right_side = "  +  ".join(part.strip() for part in re_parts[1].split("+"))
                    react = "{:^30} = {:^30}".format(left_side, right_side)
                    print("{}. {}    {:>.6}".format(i+1, react, Abs_en[1][i]))
                print("--"*40)
                if reaction_energy_list==True:
                    print(Abs_en[0])
                    print(Abs_en[1])
        else:
            print("Please create an empty E0.dat file, to create checklist file.")
    
    else:
        print("To calculate the absorption energy required files are:\n1-Reaction.in\n2-E0.dat")


def Calc_absorption_en1(reaction_energy_list=False):
    with open("energy_results.dat", 'w') as energy_results:
        print()
        print("\033[1;32mAbsorption energy:\033[0m")
        energy_results.write("Absorption energy:\n")
        if os.path.exists("Reaction.in"):
            if os.path.exists("E0.dat"):
                Abs_en = absorption_energy("Reaction.in", "E0.dat")
                if Abs_en[2]:
                    check_absorption_energy()
                else:
                    abs_energy_str = ""
                    for i in range(len(Abs_en[0])):
                        print("--"*40)
                        energy_results.write("--"*40 + "\n")
                        react = "" 
                        re_parts = Abs_en[0][i].split("=")
                        re_parts = Abs_en[0][i].split("=")
                        left_side = "  +  ".join(part.strip() for part in re_parts[0].split("+"))
                        right_side = "  +  ".join(part.strip() for part in re_parts[1].split("+"))
                        react = "{:^30} = {:^30}".format(left_side, right_side)
                        print("{}. {}    {:>.6}".format(i+1, react, Abs_en[1][i]))
                        abs_energy_str += str(Abs_en[1][i]) + " "
                        energy_results.write("{}. {}    {:>.6} \n".format(i+1, react, Abs_en[1][i]))
                    print("--"*40)
                    energy_results.write("--"*40 + "\n")
                    energy_results.write("Absorption energy list: {} \n".format(abs_energy_str))
                    if reaction_energy_list==True:
                        print(Abs_en[0])
                        print(Abs_en[1])
            else:
                print("Please create an empty E0.dat file, to create checklist file.")
        
        else:
            print("To calculate the absorption energy required files are:\n1-Reaction.in\n2-E0.dat")

def Calc_Gibbs_energy(reaction_energy_list=False):
    print()
    print("\033[1;32mGibbs free energy:\033[0m")
    if os.path.exists("Reaction1.in"):
        if os.path.exists("E0.dat") and os.path.exists("Thermal_correction.dat"):
            try:
                rxn = Rxn_Gibbs_energy("Reaction1.in", "E0.dat", "Thermal_correction.dat")
                #print(rxn)
                for i in range(len(rxn[0])):
                    print("--"*40)
                    react = "" 
                    re_parts = rxn[0][i].split("=")
                    re_parts = rxn[0][i].split("=")
                    left_side = "  +  ".join(part.strip() for part in re_parts[0].split("+"))
                    right_side = "  +  ".join(part.strip() for part in re_parts[1].split("+"))
                    react = "{:^30} = {:^30}".format(left_side, right_side)
                    print("{}. {}    {:>.6}".format(i+1, react, rxn[1][i]))
                print("--"*40)
                if reaction_energy_list==True:
                    print(rxn[0])
                    print(rxn[1])
                #print(Rxn_Gibbs_energy("Reaction1.in", "E0.dat", "Thermal_correction.dat"))
            except KeyError:
                print("KeyError")
                check_Gibbs_energy()
        else:
            print("Please create an empty E0.dat and Thermal_correction.dat files, to create checklist file.")

    else:
        print("To calculate the Gibbs energy required files are:\n1-Reaction1.in\n2-E0.dat\n3-Thermal_correction.dat")

def Calc_Gibbs_energy1(reaction_energy_list=False):
    print()
    with open("energy_results.dat", 'a') as energy_results:
        print("\033[1;32mGibbs free energy:\033[0m")
        energy_results.write("Gibbs free energy:\n")
        if os.path.exists("Reaction1.in"):
            if os.path.exists("E0.dat") and os.path.exists("Thermal_correction.dat"):
                try:
                    rxn = Rxn_Gibbs_energy("Reaction1.in", "E0.dat", "Thermal_correction.dat")
                    #print(rxn)
                    gibbs_free_str = ""
                    for i in range(len(rxn[0])):
                        print("--"*40)
                        energy_results.write("--"*40 + "\n")
                        react = "" 
                        re_parts = rxn[0][i].split("=")
                        re_parts = rxn[0][i].split("=")
                        left_side = "  +  ".join(part.strip() for part in re_parts[0].split("+"))
                        right_side = "  +  ".join(part.strip() for part in re_parts[1].split("+"))
                        react = "{:^30} = {:^30}".format(left_side, right_side)
                        print("{}. {}    {:>.6}".format(i+1, react, rxn[1][i]))
                        energy_results.write("{}. {}    {:>.6}".format(i+1, react, rxn[1][i]) + "\n")
                        gibbs_free_str += str(rxn[1][i]) + " "
                    print("--"*40)
                    energy_results.write("--"*40 + "\n")
                    energy_results.write("gibbs free energy list:{}\n".format(gibbs_free_str))
                    if reaction_energy_list==True:
                        print(rxn[0])
                        print(rxn[1])
                    #print(Rxn_Gibbs_energy("Reaction1.in", "E0.dat", "Thermal_correction.dat"))
                except KeyError:
                    print("KeyError")
                    check_Gibbs_energy()
            else:
                print("Please create an empty E0.dat and Thermal_correction.dat files, to create checklist file.")

        else:
            print("To calculate the Gibbs energy required files are:\n1-Reaction1.in\n2-E0.dat\n3-Thermal_correction.dat")

def writing_checklist(abs_checklist="", gibbs_checklist=""):
    header = [
    "S.N. (Abs)", "Absorption Energy Components", "Checklist (Abs)",
    "S.N. (Gibbs)", "Gibbs Components", "Checklist (Gibbs)"
    ]

    # Build the full data list
    data = []

    # Add 5 blank rows
    for _ in range(5):
        data.append([""] * len(header))

    # Add header row
    data.append(header)

    # Add component rows
    for i in range(max(len(abs_checklist), len(gibbs_checklist))):
        row = [
            i+1 if i < len(abs_checklist) else "",
            abs_checklist[i] if i < len(abs_checklist) else "",
            "",
            i+1 if i < len(gibbs_checklist) else "",
            gibbs_checklist[i] if i < len(gibbs_checklist) else "",
            ""
        ]
        data.append(row)

    # Write to ODS
    ods.save_data("calculation_checklist.ods", {"Sheet1": data})


def write_or_append_checklist(file_name, label, checklist):
    try:
        data = get_data(file_name)
    except FileNotFoundError:
        data = OrderedDict()
        data["Sheet1"] = []

    sheet = data.get("Sheet1", [])

    # Add 3 blank rows before new section
    sheet.extend([[""] * 3 for _ in range(3)])

    # Add header
    header = [f"S.N. ({label})", f"{label} Components", f"Checklist ({label})"]
    sheet.append(header)

    # Add checklist rows
    for i, item in enumerate(checklist, 1):
        sheet.append([i, item, ""])

    # Save updated content
    data["Sheet1"] = sheet
    save_data(file_name, data)



if __name__ == "__main__":

    Calc_absorption_en1()

    Calc_Gibbs_energy1(reaction_energy_list=True)

    #writing_checklist(abs_checklist=abs_checklist, gibbs_checklist=gibbs_checklist)
