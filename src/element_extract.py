#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to extract data from materials project database"""
import os
import glob
import warnings
import pandas as pd
from mp_api.client import MPRester
from htepc import MpConnect
from check_json import config
warnings.filterwarnings('ignore')
# Make sure that you have the Materials API key.
# MPRester if needed, e.g, MPRester("API_KEY")
#get API_KEY from material projects, go to dashboard and generate the key.

def create_folder(parent_folder):
    """
    Function to create a 'download' folder
    parameters
    --------------
    parent_folder : path to current working directory
    """
    if not os.path.isdir(parent_folder+"/download"):
        os.mkdir(parent_folder+"/download")

def download(elm,num_el,exclude_el,properties):
    """
    Extracts various properties for compounds that satisfy certain criteria from the Materials Project database.

    Parameters:
    -----------
    elm : str or list of str
        Element(s) always to include in the compounds. For example, for hydrogen, elm = 'H'.
        If multiple elements are desired, provide a list with up to size 2. For example, elm = ['B', 'C'] for boron and carbon.

    num_el : int
        Number of elements in the compound.

    exclude_el : list of str
        List of elements to exclude from the compound.

    properties : list of str
        List of properties to extract.

    Returns:
    --------
    data : pandas DataFrame
        DataFrame containing the extracted data.
    Example:
    --------
    >>> download('H', 2, ['O', 'F'], ['material_id', 'formation_energy_per_atom'])
    """
    parent_folder=os.getcwd()
    create_folder(parent_folder)
    # Check for the presence of the API key file
    if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
        key = input_data["mpi_key"]["API_KEY"]
    else:
        print("config.json file not found. Please provide with your materials project api key\n")
    # Initialize MPRester with API key
    mpr = MPRester(key["key"])
    # Search for materials matching specified criteria
    mpr_search = mpr.materials.summary.search(elements=[elm],
                                     exclude_elements=exclude_el,
                                     fields=properties,
                                     num_elements=num_el)
    # Write extracted data to CSV file
    with open(parent_folder+"/download/"+ "data-"+ elm +".csv", "w") as data_elm:
        for i,propti in enumerate(properties):
            if propti == "structure":
                propty = "spacegroup"
                data_elm.write(propty + ",")
            else:
                propty = propti
                data_elm.write(propty + ",")
            #if i < len(properties) - 1:
            #else:
            #    data_elm.write(propty)
        data_elm.write("composition\n")
        for search in mpr_search:
            property_list = []
            for propty in properties:
                if propty == "structure":
                    property_list.append(search.structure.get_space_group_info()[0])
                else:
                    property_list.append(search.dict()[propty])
            for j,prop in enumerate(property_list):
                #data_elm.write(str(prop) + ",")
                if j < len(property_list) - 1:
                    data_elm.write(str(prop) + ",")
                else:
                    data_elm.write(str(prop) + ",")
            data_elm.write(str(search.structure.composition.formula.replace(" ", "")))
            data_elm.write("\n")
    # Read the CSV file into a DataFrame
    data=pd.read_csv(parent_folder+"/download/"+"data-" + elm + ".csv")
    print(data['material_id'])
    return data

def stable(data):
    """
    Filters the compounds for those having negative formation energy.

    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing information about compounds, including formation energy per atom.

    Returns:
    --------
    data : pandas DataFrame
        DataFrame containing compounds with negative formation energy per atom.
    """
    stable_filter = data["formation_energy_per_atom"] < 0
    data = data[stable_filter]
    data = data.reset_index(drop=True)
    return data

def convexhull(data):
    """
    Filters the compounds close to the convex hull.

    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing information about compounds, including energy above hull.

    Returns:
    --------
    data : pandas DataFrame
        DataFrame containing compounds close to the convex hull.
    """
    stable_fil = data["energy_above_hull"] < 0.001
    data = data[stable_fil].reset_index(drop=True)
    return data

def metal_filter(data):
    """
    Filters metallic compounds from the input DataFrame.

    Parameters:
    -----------
    data : pandas DataFrame
        DataFrame containing information about compounds, including band_gap.

    Returns:
    --------
    data : pandas DataFrame
        DataFrame containing metallic compounds (band gap <= 0.00001).
    """
    zero_band_gap = data['band_gap'] <= 0.00001
    data = data[zero_band_gap]
    data = data.reset_index(drop=True)
    return data
def data_combine(data1,data2):
    """
    Combines two pandas DataFrames into a single DataFrame.

    Parameters:
    -----------
    data1 : pandas DataFrame
        The first DataFrame to be combined.
    data2 : pandas DataFrame
        The second DataFrame to be combined.

    Returns:
    --------
    data : pandas DataFrame
        Combined DataFrame containing data from both data1 and data2.
    """
    data = pd.merge(data1, data2, how='outer')
    data = data.reset_index(drop=True)
    return data
#def data_2_prefix(data):
#    prefix = []
#    for j in range(data.shape[0]):
#        s = ""
#        elm = list(data['composition'][j].keys())
#        count = list(data['composition'][j].values())
#        for i in range(len(elm)):
#            s += elm[i]+str(int(count[i]))
#        prefix.append(s)
#    data.drop(columns=['composition'],axis=1)
#    data['composition'] = prefix
#    return data
def remove(data,element_list):
    """
    Removes compounds containing specified elements from the DataFrame.

    Parameters:
    -----------
    data : pandas DataFrame
        The DataFrame containing compounds to be filtered.
    element_list : str, optional
        File with elements to exclude. Default is 'remove.list'.
        The file should contain elements separated by commas.
        For example, to remove oxygen and nitrogen, write 'O,N' in 'remove.list'.

    Returns:
    --------
    data : pandas DataFrame
        Processed DataFrame with compounds containing specified elements removed.
    """
    if not os.path.isfile(element_list):
        os.system("""echo "NA" > remove.list""")
    with open(element_list, "r") as read_remove:
        lines = read_remove.readlines()
    remove_elements = lines[0].replace("\n", "").split(',')
    pattern_remove = '|'.join(remove_elements)
    print(pattern_remove)
    filter_remove = data["formula_pretty"].str.contains(pattern_remove)
    filter_temp = []
    for rem in filter_remove:
        filter_temp.append(not rem)
    data = data[filter_temp].reset_index(drop=True)
    #data.to_csv(filename)
    return data

def data_one_element_compound(elm,ntype,exclude_el,properties):
    """
    Extracts information for compounds containing only one element.

    Parameters:
    -----------
    elm : str
        Element to search for in compounds. For example, 'B' for boron.
    ntype : int or tuple
        Number of unique elements. Can be a single integer or a tuple (e.g., (1, 2) for 2 different types).
    exclude_el : list
        List of elements to exclude. For example, ['O', 'N'].
    properties : list
        List of properties to extract.

    Returns:
    --------
    data : pandas DataFrame
        DataFrame containing information for compounds with only one element.
    """
    data = download(elm,ntype,exclude_el,properties)
    #data = download(el,ntype,properties)
    #data = stable(data)
    #data = metal_filter(data)
    #data = remove(data,'remove.list')
    return data

def data_two_element_compound(el1,el2,ntype,exclude_el,properties):
    """
    Extracts information for compounds containing two elements.

    Parameters:
    -----------
    el1 : str
        First element to search for in compounds (e.g., 'B' for boron).
    el2 : str
        Second element to search for in compounds (e.g., 'C' for carbon).
    ntype : int or tuple
        Number of unique elements in compounds. Can be a single integer or a tuple (e.g., (1, 2) for 2 different types).
    exclude_el : list
        List of elements to exclude. For example, ['O', 'N'].
    properties : list
        List of properties to extract.

    Returns:
    --------
    data : pandas DataFrame
        DataFrame containing information for compounds with two elements.
    """
    data1=download(el1,ntype,exclude_el,properties)
    data2=download(el2,ntype,exclude_el,properties)
    data = data_combine(data1,data2)
    #data = stable(data)
    #data = metal_filter(data)
    #data = remove(data,'remove.list')
    return data

def create_input():
    """
    Reads 'download.csv' file inside 'download' folder and creates 'input.in' and 'mpid-list.in' files
    for further downloading and calculations.
    """
    data_file = pd.DataFrame(pd.read_csv('download/download.csv'))
    nrow = data_file.shape[0]
    with open("mpid-list.in", "a") as mpfile_append:
        for i in range(nrow):
            mpid = data_file['material_id'][i]
            #comp = data_file['formula_pretty'][i]
            comp = data_file['composition'][i]
            mpfile_append.write("v{} {} {}".format(i+1,mpid,comp) + "\n")
    with open('input.in', 'w') as input_write:
        if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
            #import download as d
            d = input_data['download']
            input_write.write(str(d['inp']['start']) + "\n")
            input_write.write(str(d['inp']['end']) + "\n")
            input_write.write("{} 0".format(d['inp']['nkpt']) + "\n")
            input_write.write("mpid-list.in\n")
            input_write.write("{}".format(d['inp']['plot']) + "\n")
            input_write.write("DFT = {}".format(d['inp']['calc']) + "\n")
        else:
            input_write.write(str(1) + "\n")
            input_write.write(str(nrow) + "\n")
            input_write.write("200 0\n")
            input_write.write("mpid-list.in\n")
            input_write.write("phband\n")
            input_write.write("DFT = QE\n")

def extract(ntype,properties,elm,exclude_el,nelm=1,metal=False,
            neg_fe=False,
            thermo_stable=False,
            ordering='NM',
            nsites=10,
            spacegroup=None,
            out='download/download.csv'):
    """
    Function to extract the data and apply filters, then write 'download.csv' file inside 'download' folder.

    Parameters:
    -----------
    ntype : tuple
        Number of unique elements in the compound. For example: (1, 3) for 3 different unique elements in compounds.
    properties : list
        List of properties to extract.
    elm : list
        List of elements used in search.
    exclude_el : list
        List of elements to exclude.
    nelm : int, optional
        Length of list elm. Default is 1.
    metal : bool, optional
        True to download zero bandgap compounds. Default is False.
    neg_fe : bool, optional
        True to download compounds with negative formation energy. Default is False.
    thermo_stable : bool, optional
        True to download compounds at the convex hull. Default is False.
    ordering : str, optional
        Magnetic ordering of the compound. Default is 'NM'.
    nsites : int, optional
        Maximum number of sites in the compound. Default is 10.
    spacegroup : int or str, optional
        Spacegroup number or name. Default is None.
    out : str, optional
        Output file to write. Default is 'download/download.csv'.

    Returns:
    --------
    data : pandas DataFrame
        Extracted data after applying filters.
    """
    if os.path.isdir("download_old"):
        os.system("rm -r download_old")
    if os.path.isdir("download"):
        print("A download folder is found, renaming download_old\n")
        os.system("mv download download_old")
    if nelm == 1:
        data = data_one_element_compound(elm[0],ntype,exclude_el,properties)
    elif nelm == 2:
        data = data_two_element_compound(elm[0],elm[1],ntype,exclude_el,properties)
    else:
        print("Upto 2 elements are allowed\n")
    if metal:
        data = metal_filter(data)
    if neg_fe:
        data = stable(data)
    if thermo_stable:
        data = convexhull(data)
    data = data[data['ordering'] == ordering].reset_index(drop=True)
    data = data[data['nsites'] <= nsites].reset_index(drop=True)
    if spacegroup:
        data = data[data['spacegroup'] == spacegroup].reset_index(drop=True)
    os.system("rm download/data*")
    #data = remove(data, 'remove.list',elm)
    data.to_csv(out)
    return data

def download_by_entry(entries,must_include,size_constraint=20,ntype_constraint=5,FE=False,thermo_stable=True,metal=False,magnetic=False,spacegroup=None,properties=None):
    """
    Function to extract and create input files using "mp_api.client.MPRester.get_entries_in_chemsys" Function of the materials project API package (pip install mp_api).
    This mode is turned on when using 'mode':'chemsys' in 'download.py' file.

    Parameters:
    -----------
    entries : list
        List of elements ==> elements and compounds (combination of elements) to search.
    size_constraint : int, optional
        Size of the compounds (total number of ions). Upper bound not included. Default is 20.
    ntype_constraint : int, optional
        Number of different types of ions. Upper bound not included. Default is 5.
    must_include : list
        Elements that must be included in the compounds.
    FE : bool, optional
        True if the formation energy is negative. Default is False.
    metal : bool, optional
        True if the compound is a metal. Default is False.
    magnetic : bool, optional
        True if the compound has a non-zero magnetic moment. Default is False.
    spacegroup : int or str, optional
        Spacegroup number or name. Default is None.
    properties : list, optional
        List of properties to extract.

    Returns:
    --------
    None
    """
    if os.path.isdir("download_old"):
        os.system("rm -r download_old")
    if os.path.isdir("download"):
        print("A download folder is found, renaming download_old\n")
        os.system("mv download download_old")
    parent_folder=os.getcwd()
    create_folder(parent_folder)
    obj = MpConnect()
    # Create condition for must_include elements
    must_in = ""
    for i,elm in enumerate(must_include):
        if i == len(must_include) - 1:
            must_in += "'{}'".format(elm) + " in elm_list"
        else:
            must_in += "'{}'".format(elm) + " in elm_list or "
    must_in = "nelm < ntype_constraint and " + "({})".format(must_in)
    # Get entries in chemical system
    entries = obj.mpr.get_entries_in_chemsys(entries)
    entry = 1
    # Write header to CSV file
    with open(parent_folder+"/download/"+ "download.csv", "w") as data_elm:
        for i,propti in enumerate(properties):
            if propti == "structure":
                propty = "spacegroup"
            else:
                propty = propti
            if i < len(properties) - 1:
                data_elm.write(propty + ",")
            else:
                data_elm.write(propty)
        data_elm.write("\n")
    # Iterate over entries
    for i,_ in enumerate(entries):
        # Extract data for each entries
        mpid = entries[i].data['material_id']
        obj.setting(mpid)
        band_gap = obj.data['band_gap']
        form_energy = obj.data['formation_energy_per_atom']
        ordering = obj.data['ordering']
        nelm = len(entries[i].composition.elements)
        comp = entries[i].composition.formula.replace(' ','')
        count = int(entries[i].composition.num_atoms)
        elm_list = list(entries[i].composition.as_dict().keys())
        energy_above_hull = obj.data['energy_above_hull']
        # Define logic for filtering based on optional parameters
        if isinstance(thermo_stable, bool):
            if thermo_stable:
                thermo_logic = energy_above_hull < 0.0001
            else:
                thermo_logic = True
        elif isinstance(thermo_stable, (int, float)):
            thermo_logic = energy_above_hull < thermo_stable
        else:
            thermo_logic = True
        if metal:
            gap_logic = band_gap < 0.0001
        else:
            gap_logic = True
        if FE:
            fe_logic = form_energy < 0.0
        else:
            fe_logic = True
        if not magnetic:
            mag_logic = ordering == 'NM'
        else:
            mag_logic = True
        if not spacegroup:
            sg_logic = True
        else:
            sg_logic = obj.data['symmetry']['symbol'] == spacegroup
        print("Extracting {}".format(comp) + "\n")
        # Write data to input file and CSV file if conditions are met
        with open("mpid-list.in", "a") as mplist_write:
            if count < size_constraint and gap_logic and fe_logic and mag_logic and sg_logic and thermo_logic:
                if eval(must_in):
                    mplist_write.write("v{} {} {}".format(entry,mpid,obj.prefix) + "\n")
                    entry += 1
        property_list = []
        with open(parent_folder+"/download/"+ "download.csv", "a") as data_elm:
            for propty in properties:
                if propty == "structure":
                    property_list.append(obj.data['symmetry']['symbol'])
                else:
                    property_list.append(obj.data[propty])
            for j,prop in enumerate(property_list):
                if j < len(property_list) - 1:
                    data_elm.write(str(prop) + ",")
                else:
                    data_elm.write(str(prop))
            data_elm.write("\n")
def main():
    """
    Main function to orchestrate the data extraction and input file creation process.

    If 'mpid-list.in' file does not exist, the function reads settings from 'config.json' to
    create mpid-list.in file.

    Parameters:
    -----------------
    None

    Returns:
    -----------------
    None
    """
    # Check if 'mpid-list.in' file exists
    if not os.path.isfile("mpid-list.in"):
        # Check if 'config.json' exists
        if os.path.isfile("config.json") or os.path.isfile("../../config.json"):
            # Read settings from 'config.json' to create 'mpid-list.in' file
            d = input_data["download"]
            mode = d['mode']
            ntype = d['element']['ntype']
            exclude_el = d['element']['exclude']
            elm_list = d['element']['elm']
            nelm = len(elm_list)
            properties = d['element']['prop']
            metal=d['element']['metal']
            neg_fe=d['element']['FE']
            thermo_stable=d['element']['thermo_stable']
            ordering=d['element']['ordering']
            nsites=d['element']['nsites']
            spacegroup=d['element']['spacegroup']
        else:
            # Provide default settings if 'config.json' doesn't exist
            print("input file config.json not found\n")
            print("Create one with following format\n")
            msg="""element={'metal':True, 'FE':True, 'exclude':["O", "N", "F", "Cl", "Br", "I"],'ntype':(1,2), 'elm':['B'], 'prop':["material_id", "formula_pretty", "structure", "formation_energy_per_atom", "band_gap", "energy_above_hull","nsites","ordering","nsites"],'ordering':'NM','nsites':10,'spacegroup':None}
inp=    {'start':1, 'end':50, 'nkpt':200, 'evenkpt': False, 'plot':'phband', 'calc':'QE'}
chemsys={'entries':['B'],'size_constraint':20,'ntype_constraint':5,'must_include':['Mg'],'FE':False,'metal':False,'magnetic':False,'spacegroup':None}"""
            print(msg + "\n")
            print("Utilizing default settings\n")
            ntype = (1,2) #Number of different types of element in the compound.
            exclude_el = ["Lu"]
            nsites = 10
            #exclude_el = ["O", "N", "F", "Cl", "Br", "I"]
            elm = 'B'
            nelm = 1
            elm_list = [elm]
            metal = False
            neg_fe = False
            thermo_stable = False
            ordering = 'FM'
            spacegroup = None
            properties=["material_id", "formula_pretty", "structure", "formation_energy_per_atom", "band_gap", "energy_above_hull","total_magnetization","ordering",'total_magnetization_normalized_formula_units', 'num_magnetic_sites','theoretical','nsites']
            default1={'metal':metal,'FE':neg_fe, 'thermo_stable':thermo_stable, 'exclude':exclude_el,'ntype':(1,2),'elm':[elm],'prop':properties,'ordering':ordering,'nsites':nsites,'spacegroup':spacegroup}
            default2={'start':1, 'end':2, 'nkpt':200, 'evenkpt': False, 'plot':'phband','calc':'QE'}
            chemsys={'entries':['B','Mg'],'size_constraint':20,'ntype_constraint':5,'must_include':['Mg','B'],'FE':False,'metal':False, 'magnetic':False,'spacegroup':spacegroup}
            d = {
                 'element':default1,
                 'inp':default2,
                 'chemsys':chemsys
                 }
            # Default mode is 'element'
            mode = 'element'
        # Perform actions based on mode
        if mode == 'element':
            # Extract data and create input files
            data = extract(ntype,properties,elm_list,exclude_el,
                           nelm=nelm,
                           metal=metal,
                           neg_fe=neg_fe,
                           thermo_stable=thermo_stable,
                           ordering=ordering,
                           nsites=nsites,
                           spacegroup=spacegroup)
            create_input()
        elif mode == 'chemsys':
            # Download data for compounds based on chemical system
            download_by_entry(d['chemsys']['entries'],d['chemsys']['must_include'],d['chemsys']['size_constraint'],d['chemsys']['ntype_constraint'],
                              d['chemsys']['FE'],d['chemsys']['thermo_stable'],d['chemsys']['metal'],d['chemsys']['magnetic'],d['chemsys']['spacegroup'],properties)
        elif mode == 'fromcif':
            # List CIF files
            list_cif = glob.glob("*.cif",recursive=True)
            if len(list_cif) > 0:
                print("These cif files are found\n")
                for cif in list_cif:
                    print(cif + "\n")
        elif mode == 'fromvasp':
            # List VASP files
            list_vasp = glob.glob("*.vasp",recursive=True)
            if len(list_vasp) > 0:
                print("These .vasp files are found\n")
                for vasp in list_vasp:
                    print(vasp + "\n")
        else:
            print("mode = element, chemsys, fromcif, or fromvasp available\n")
if __name__ == "__main__":
    input_data = config()
    main()
