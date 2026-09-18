#!/usr/bin/env python
#"""Writen by Niraj K. Nepal, Ph.D."""
"""Module to extract data from materials project database"""
import os
import glob
from enum import Enum
import shutil
import warnings
import pandas as pd
from htesp.htepc import MpConnect, mprester
from htesp.check_json import config, require_api_key
warnings.filterwarnings('ignore')
# Make sure that you have the Materials API key.
# MPRester if needed, e.g, MPRester("API_KEY")
#get API_KEY from material projects, go to dashboard and generate the key.

def legacy_mpid(value):
    """A Materials Project id in the legacy ``mp-<integer>`` spelling.

    FIX(29): Materials Project is migrating identifiers from ``MPID``
    (``mp-763``) to ``AlphaID`` (``mp-bdj``).  The two spellings address the
    same material and the API accepts either, but only the legacy one matches
    the ``R<mpid>-<compound>/`` directories, ``scf_dir/scf-<mpid>.in`` and the
    tracking files of a campaign started before the migration -- and it is
    what a person recognises.

    The 1.x code wrote ``search['material_id'].string``.  The rewrite replaced
    that with ``search.dict()['material_id']``, which drops two things:

    * ``.dict()`` serialises the id to a **plain str already in the alphabetic
      spelling** -- ``'mp-bdj'``, with no ``.string`` attribute left to
      recover the legacy form.  Attribute access (``search.material_id``)
      keeps the ``MPID``/``AlphaID`` object, so callers should prefer it;
    * ``.string``, which is the conversion itself.

    This helper therefore tries, in order: the object's own ``.string``;
    decoding the text through :class:`~emmet.core.mpid.AlphaID`, which accepts
    either spelling and preserves the prefix (``mvc-`` stays ``mvc-``); and
    finally the text unchanged, so a non-MP identifier (OQMD, AFLOW) passes
    through untouched.
    """
    string = getattr(value, "string", None)
    if string:
        return str(string)
    text = str(value)
    if "-" not in text:
        return text
    try:
        from emmet.core.mpid import AlphaID
        return str(AlphaID(text).string)
    except Exception:            # not an AlphaID-shaped id, or no emmet-core
        return text


def plain_value(value):
    """An ``Enum`` member -> its value; anything else unchanged.

    FIX(28): ``emmet.core`` declares ``Ordering`` as a *plain* ``Enum``, so
    ``Ordering.NM`` is not equal to ``"NM"`` and ``str(Ordering.NM)`` is
    ``"Ordering.NM"``.  Both bite the magnetic-ordering filter:

    * ``download()`` wrote ``str(prop)`` into ``download/data-<elm>.csv``, so the
      ``ordering`` column held ``Ordering.NM`` and the ``extract()`` filter
      ``data['ordering'] == 'NM'`` matched **no row at all** -- a search that
      returned 465 compounds wrote 0 of them to mpid-list.in;
    * ``download_by_entry()`` compares the raw field, so ``mag_logic`` was
      False for every entry whenever ``chemsys.magnetic`` was false.

    Normalising here keeps the CSV and the comparisons in the documented
    ``"NM"`` / ``"FM"`` / ``"AFM"`` / ``"FiM"`` vocabulary of config.json.
    """
    return value.value if isinstance(value, Enum) else value


def filter_ordering(data, ordering):
    """Keep the rows whose magnetic ordering matches ``ordering``.

    FIX(30): this was an unconditional ``data['ordering'] == ordering``, the
    only filter in :func:`extract` with no way to switch it off -- ``metal``,
    ``FE`` and ``thermo_stable`` are each guarded by an ``if``.  That was
    survivable while Materials Project reported a definite ordering for
    everything, but it now reports ``"Unknown"`` for any material with no
    magnetism calculation: in a B-binary search, 202 of 465 hits, and 158 of
    the 222 that survive the metal and formation-energy filters.  With
    ``ordering: "NM"`` every one of those is discarded silently, which is how
    a search returning 465 compounds wrote 2 rows to mpid-list.in.

    Accepted values of ``ordering``:

    ``null`` / ``None``
        apply no filter at all -- keep every ordering;
    a string (``"NM"``)
        keep exactly that ordering, the historical behaviour;
    a list (``["NM", "Unknown"]``)
        keep any row whose ordering is one of them.  This is usually what is
        wanted: "Unknown" mostly means nobody ran the magnetism calculation,
        not that the material is magnetic.

    An empty list is treated as ``null`` rather than as "match nothing", since
    a filter that rejects everything is never what a configuration means.
    """
    if ordering is None:
        return data.reset_index(drop=True)
    if isinstance(ordering, (list, tuple, set, frozenset)):
        wanted = [str(plain_value(item)) for item in ordering]
        if not wanted:
            return data.reset_index(drop=True)
        return data[data['ordering'].isin(wanted)].reset_index(drop=True)
    return data[data['ordering'] == str(plain_value(ordering))].reset_index(drop=True)


def create_folder(parent_folder):
    """
    Function to create a 'download' folder
    parameters
    --------------
    parent_folder : path to current working directory
    """
    # FIX(all): os.makedirs(..., exist_ok=True) instead of isdir-then-mkdir
    os.makedirs(os.path.join(parent_folder, "download"), exist_ok=True)

def download(elm,num_el,exclude_el,properties,input_data=None):
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

    input_data : dict, optional
        Parsed ``config.json``.  Loaded with :func:`htesp.check_json.config`
        when omitted.

    Returns:
    --------
    data : pandas DataFrame
        DataFrame containing the extracted data.
    Example:
    --------
    >>> download('H', 2, ['O', 'F'], ['material_id', 'formation_energy_per_atom'])
    """
    # FIX(1): ``input_data`` used to be a module-level name bound only inside
    # ``if __name__ == "__main__":``, so importing this module and calling
    # download() raised NameError.  The configuration is loaded here instead.
    if input_data is None:
        input_data = config()
    parent_folder=os.getcwd()
    create_folder(parent_folder)
    # FIX(1): the key now comes from $MP_API_KEY / ~/.config/htesp/credentials
    # / config.json via require_api_key(), which raises a message that says
    # how to set it instead of leaving ``key`` unbound.
    mpr = mprester(require_api_key(input_data))
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
                elif propty == "material_id":
                    # FIX(3): an MPDataDoc is not subscriptable, so 1.x's
                    # ``search['material_id']`` no longer works.
                    # FIX(29): read the *attribute*, not ``.dict()[...]`` --
                    # the dict form is already a plain 'mp-bdj' string with no
                    # ``.string`` left on it, while the attribute is still an
                    # MPID whose legacy spelling legacy_mpid() can take.
                    property_list.append(
                        legacy_mpid(getattr(search, propty, search.dict()[propty])))
                else:
                    # FIX(28): enum -> its value, so 'ordering' reaches the CSV
                    # as 'NM' rather than 'Ordering.NM'
                    property_list.append(plain_value(search.dict()[propty]))
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
        # FIX(all): os.system replaced by a plain file write
        with open(element_list, "w") as write_remove:
            write_remove.write("NA\n")
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

def data_one_element_compound(elm,ntype,exclude_el,properties,input_data=None):
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
    data = download(elm,ntype,exclude_el,properties,input_data)
    #data = download(el,ntype,properties)
    #data = stable(data)
    #data = metal_filter(data)
    #data = remove(data,'remove.list')
    return data

def data_two_element_compound(el1,el2,ntype,exclude_el,properties,input_data=None):
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
    data1=download(el1,ntype,exclude_el,properties,input_data)
    data2=download(el2,ntype,exclude_el,properties,input_data)
    data = data_combine(data1,data2)
    #data = stable(data)
    #data = metal_filter(data)
    #data = remove(data,'remove.list')
    return data

def create_input(input_data=None):
    """
    Reads 'download.csv' file inside 'download' folder and creates 'input.in' and 'mpid-list.in' files
    for further downloading and calculations.
    """
    # FIX(1): ``input_data`` was a name bound only under ``__main__``.
    if input_data is None:
        input_data = config()
    data_file = pd.DataFrame(pd.read_csv('download/download.csv'))
    nrow = data_file.shape[0]
    with open("mpid-list.in", "a") as mpfile_append:
        for i in range(nrow):
            mpid = data_file['material_id'][i]
            #comp = data_file['formula_pretty'][i]
            comp = data_file['composition'][i]
            mpfile_append.write("v{} {} {}".format(i+1,mpid,comp) + "\n")
    # FIX(1): check_json.config() now always returns a fully populated dict,
    # so the 'config.json exists?' branch is gone; a missing key still falls
    # back to the documented defaults.
    try:
        inp = input_data['download']['inp']
    except (KeyError, TypeError):
        inp = {}
    with open('input.in', 'w') as input_write:
        input_write.write(str(inp.get('start', 1)) + "\n")
        input_write.write(str(inp.get('end', nrow)) + "\n")
        input_write.write("{} 0".format(inp.get('nkpt', 200)) + "\n")
        input_write.write("mpid-list.in\n")
        input_write.write("{}".format(inp.get('plot', 'phband')) + "\n")
        input_write.write("DFT = {}".format(inp.get('calc', 'QE')) + "\n")

def extract(ntype,properties,elm,exclude_el,nelm=1,metal=False,
            neg_fe=False,
            thermo_stable=False,
            ordering='NM',
            nsites=10,
            spacegroup=None,
            out='download/download.csv',
            input_data=None):
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
    ordering : str, list or None, optional
        Magnetic ordering to keep. Default is 'NM'. ``None`` applies no
        filter; a list keeps any of its values, e.g. ``["NM", "Unknown"]``
        (Materials Project reports ``"Unknown"`` whenever no magnetism
        calculation was run). See :func:`filter_ordering`.
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
    # FIX(1): configuration is loaded here rather than read from a module
    # global that only ``__main__`` used to set.
    if input_data is None:
        input_data = config()
    # FIX(all): os.system() replaced by pure-Python file operations, so a
    # folder name can never be interpreted by a shell.
    if os.path.isdir("download_old"):
        shutil.rmtree("download_old")
    if os.path.isdir("download"):
        print("A download folder is found, renaming download_old\n")
        shutil.move("download", "download_old")
    if nelm == 1:
        data = data_one_element_compound(elm[0],ntype,exclude_el,properties,input_data)
    elif nelm == 2:
        data = data_two_element_compound(elm[0],elm[1],ntype,exclude_el,properties,input_data)
    else:
        raise ValueError("Upto 2 elements are allowed, got {}".format(nelm))
    if metal:
        data = metal_filter(data)
    if neg_fe:
        data = stable(data)
    if thermo_stable:
        data = convexhull(data)
    # FIX(30): null skips the filter, a list keeps any of its values;
    # see filter_ordering() for why this had to become optional.
    data = filter_ordering(data, ordering)
    data = data[data['nsites'] <= nsites].reset_index(drop=True)
    if spacegroup:
        data = data[data['spacegroup'] == spacegroup].reset_index(drop=True)
    # FIX(all): os.system("rm download/data*") -> glob + os.remove
    for stale in glob.glob(os.path.join("download", "data*")):
        os.remove(stale)
    #data = remove(data, 'remove.list',elm)
    data.to_csv(out)
    return data

def download_by_entry(entries,must_include,size_constraint=20,ntype_constraint=5,FE=False,thermo_stable=True,metal=False,magnetic=False,spacegroup=None,properties=None,input_data=None):
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
    # FIX(all): os.system() replaced by pure-Python file operations, so a
    # folder name can never be interpreted by a shell.
    if os.path.isdir("download_old"):
        shutil.rmtree("download_old")
    if os.path.isdir("download"):
        print("A download folder is found, renaming download_old\n")
        shutil.move("download", "download_old")
    parent_folder=os.getcwd()
    create_folder(parent_folder)
    if input_data is None:
        input_data = config()
    obj = MpConnect()
    # FIX(2): the filter used to be assembled as a Python source string
    # ("nelm < ntype_constraint and ('Mg' in elm_list or 'B' in elm_list)")
    # and handed to eval(), i.e. arbitrary config.json content was executed.
    # The predicate it emulated is written out explicitly below; note that an
    # empty must_include now selects nothing, where eval() raised SyntaxError.
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
        # FIX(29): legacy mp-<integer> spelling (see legacy_mpid)
        mpid = legacy_mpid(entries[i].data['material_id'])
        obj.setting(mpid)
        band_gap = obj.data['band_gap']
        form_energy = obj.data['formation_energy_per_atom']
        # FIX(28): Ordering.NM != 'NM'; compare the value
        ordering = plain_value(obj.data['ordering'])
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
                # FIX(2): explicit predicate, replacing eval(must_in)
                if nelm < ntype_constraint and any(e in elm_list for e in must_include):
                    mplist_write.write("v{} {} {}".format(entry,mpid,obj.prefix) + "\n")
                    entry += 1
        property_list = []
        with open(parent_folder+"/download/"+ "download.csv", "a") as data_elm:
            for propty in properties:
                if propty == "structure":
                    property_list.append(obj.data['symmetry']['symbol'])
                elif propty == "material_id":
                    property_list.append(mpid)
                else:
                    # FIX(28): enum -> its value (see plain_value)
                    property_list.append(plain_value(obj.data[propty]))
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
    # FIX(1): the configuration is loaded once, here, instead of being read
    # from a module global bound only under ``if __name__ == "__main__":``.
    input_data = config()
    # Check if 'mpid-list.in' file exists
    if not os.path.isfile("mpid-list.in"):
        # check_json.config() always returns a populated dict now; only a
        # structurally broken file falls through to the defaults below.
        try:
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
        except (KeyError, TypeError):
            # Provide default settings if 'config.json' is unusable
            print("usable download settings not found in config.json\n")
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
                           spacegroup=spacegroup,
                           input_data=input_data)
            create_input(input_data)
        elif mode == 'chemsys':
            # Download data for compounds based on chemical system
            download_by_entry(d['chemsys']['entries'],d['chemsys']['must_include'],d['chemsys']['size_constraint'],d['chemsys']['ntype_constraint'],
                              d['chemsys']['FE'],d['chemsys']['thermo_stable'],d['chemsys']['metal'],d['chemsys']['magnetic'],d['chemsys']['spacegroup'],properties,input_data)
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
    # FIX(1): main() loads the configuration itself, so no module global.
    main()
