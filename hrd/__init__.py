# -*- coding: utf-8 -*-
"""
Created on  : Mon Oct 24 17:32:59 2016
Author      : Guus Rongen
Project     : PR0000.00
Description :

"""

import sqlite3
import os
import pandas as pd
import shutil
import numpy as np
import datetime
import re
import linecache

def create_empty_copy(srcdb, dstdb):
    """
    Create an empty copy of a database, that can be used to fill a new
    database with the same configuration.

    Parameters
    ----------
    srcdb : str
        path to database from which the empty copy is made
    dstdb : str
        path where to save the new empty database
    """

    # Check if the source database exists
    if not os.path.exists(srcdb):
        raise NameError('The system cannot find the path specified: \'{}\''.format(srcdb))

    # Check if the destination database exists
    if os.path.exists(dstdb):
        raise NameError('The path: \'{}\' already exists.'.format(dstdb))

    # Connect to existing database
    conn = sqlite3.connect(srcdb)
    # Connect to new database
    connnew = sqlite3.connect(dstdb)

    # Haal de tabelstructuur op
    sql_tbl_struct = "SELECT sql FROM sqlite_master WHERE type='table'"
    sql_createtbl  = conn.execute(sql_tbl_struct).fetchall()

    # Maak de nieuwe tabel
    for itbl in sql_createtbl:
        connnew.execute(itbl[0])

    # Sluit de databases af
    connnew.close()
    conn.close()

def export_locations(db, add_db_name = False):
    """
    Function to export the HRDLocation tabel from the HRD database

    Parameters
    ----------
    db : str
        Path to source database
    add_db_name : boolean
        Add the database name as an extra column to the table

    Returns
    -------
    pandas dataframe
    """

    # Check if the source database exists
    if not os.path.exists(db):
        raise NameError('The system cannot find the path specified: \'{}\''.format(db))

    # Connect to existing database
    conn = sqlite3.connect(db)
    query = r'SELECT * FROM HRDLocations;'
    # Get location table
    tablocaties = pd.read_sql(query, conn)
    # Close database
    conn.close()
    # Add column with database name
    if add_db_name:
        tablocaties['Database'] = db.split('\\')[-1]

    return tablocaties

def find_HRDLocationId(name, conn):
    """
    Find HRDLocationId that matches with a name

    Parameters
    ----------
    name : string
        Name as present in table "HRDLocations", column "Name"
    conn : SQLite-connection / path
        Default None. If None, a new connection is made with the given path.

    Returns
    -------
    HRDLocationID
    """

    # If a new connection should be made
    if isinstance(conn, str):
        conn = sqlite3.connect(conn)
        locid = conn.execute('SELECT HRDLocationId FROM HRDLocations WHERE Name=?;', (name,)).fetchall()[0][0]
        conn.close()
    # If a connection is given
    elif isinstance(conn, sqlite3.Connection):
        locid = conn.execute('SELECT HRDLocationId FROM HRDLocations WHERE Name=?;', (name,)).fetchall()[0][0]

    # Else
    else:
        raise TypeError("Type of input argument \'conn\' should be a path (string) or sqlite3.Connection")

    return locid


def slice_from_db(srcdatabases, dstdatabase, locationids, overwrite = False, emptydatabase=None, drop_duplicates=False):
    """
    Function to take a slice from a one or more existing HRD databases.

    Parameters
    ---------
    srcdatabases : (list of) paths
        database(s) from which the locations are selected
    dstdatabase : path
        location on which to create the new database
    locationids: list of integers
        list with HRDLocationId's that are usedto select the data
    overwrite : bool
        Whether to overwrite the existing database (default = False)
    emptydatabase : path
        Path with an empty database with the same schema. If this is given no
        copy has to be made, which can save processingtime (default = None)
    drop_duplicates : boolean
        Whether to drop duplicate HydroDynamicInputData per location, closure
        situation and wind direction (default = False)
    """

    # Delete the exisiting table if overwrite is True
    if os.path.exists(dstdatabase) and overwrite:
        os.remove(dstdatabase)
    elif os.path.exists(dstdatabase) and not overwrite:
        print('"{}" already exists. Set overwrite to True if you want to overwrite the existing database.'.format(dstdatabase))
        return None

    # Make srcdatabases a list if only one path is given
    if type(srcdatabases) is str:
        srcdatabases = [srcdatabases]

    # Maak een lege database om te vullen
    if not emptydatabase:
        create_empty_copy(srcdatabases[0], dstdatabase)
    # Copy the empty database
    else:
        shutil.copy2(emptydatabase, dstdatabase)

    # Connect to the new database
    conn = sqlite3.connect(dstdatabase)
    # Get the tablenames from the new table
    table_names_query = "SELECT tbl_name FROM sqlite_master WHERE type='table'"
    table_names = conn.execute(table_names_query).fetchall()
    table_names = [naam[0] for naam in table_names]
    conn.close()

    # Create query based on location id's
    if not isinstance(locationids[0], list):
        locatiequery = r'SELECT * FROM {tn} WHERE LocationId IN ('+', '.join(['{}'.format(i) for i in locationids]) + ');'
        hrdlocatiequery = r'SELECT * FROM {tn} WHERE HRDLocationId IN ('+', '.join(['{}'.format(i) for i in locationids]) + ');'

    # Attach the database
    for dbnum, srcdatabase in enumerate(srcdatabases):

        # If de locationids are a list of lists:
        if isinstance(locationids[0], list):
            locatiequery = r'SELECT * FROM {tn} WHERE LocationId IN ('+', '.join(['{}'.format(i) for i in locationids[dbnum]]) + ');'
            hrdlocatiequery = r'SELECT * FROM {tn} WHERE HRDLocationId IN ('+', '.join(['{}'.format(i) for i in locationids[dbnum]]) + ');'

        # Check if the source database exists
        if not os.path.exists(srcdatabase):
            print('"{}" not found.'.format(srcdatabase))
            continue

        # Connect to database
        conn = sqlite3.connect(srcdatabase)
        # Attach the database
        conn.execute("ATTACH DATABASE ? AS db2", (dstdatabase,))

        # Get HydroDynamicData
        HydroDynamicData = pd.read_sql(hrdlocatiequery.format(tn = 'HydroDynamicData'), conn)

        if not drop_duplicates:
            # Select the HydroDynamicDataIds from the table
            hydrodynamicdataids = HydroDynamicData.HydroDynamicDataId.values.tolist()

        else:
            hydrodynamicdataids = []
            # Select the hydro dynamic data ids where the entires are unique
            for index, group in HydroDynamicData.groupby(['HRDLocationId', 'ClosingSituationId', 'HRDWindDirectionId']):
                # Select de HydroDynamicDataIds for each combination of the location and discrete stochasts
                group_hddids = group.HydroDynamicDataId
                # Select HydroDynamicInputData
                HydroDynamicInputData = pd.read_sql('SELECT * FROM HydroDynamicInputData WHERE HydroDynamicDataId IN ('+', '.join(['{}'.format(i) for i in group_hddids]) + ');', conn)
                # Use 'HydroDynamicDataId', 'HRDInputColumnId' as index, and unstack.
                # This gives a table with every InputData combination as entry
                HydroDynamicInputData = HydroDynamicInputData.set_index(['HydroDynamicDataId', 'HRDInputColumnId']).unstack()
                # Drop the duplicates and select the index, add this directly to the list
                hydrodynamicdataids += HydroDynamicInputData.drop_duplicates().index.tolist()




        # Create the Hydrodataquery
        hydrodataquery = r'SELECT * FROM {tn} WHERE HydroDynamicDataId IN (' + ','.join(['{}'.format(hdid) for hdid in hydrodynamicdataids]) + ');'


        # Loop through tables and copy a selection
        for tabel in table_names:
            columns = conn.execute('PRAGMA table_info({});'.format(tabel)).fetchall()
            columns = [col[1] for col in columns]
            # Determine what needs to be collected from the table
            # If the table has a column HRDLocationId, use the locatiequery
            if 'HRDLocationId' in columns:
                query = hrdlocatiequery.format(tn = tabel)
            # Europoort contains a seiches table with LocationId instead of HRDLocationId
            elif 'LocationId' in columns:
                query = locatiequery.format(tn = tabel)
            # Else if the table has a column HRDDynamicDataId, use the hydrodataquery
            elif 'HydroDynamicDataId' in columns:
                query = hydrodataquery.format(tn = tabel)
            # Else, get all
            else:
                # Similar tables only have to be copied once, the first time.
                if dbnum > 0:
                    continue
                # It is the first round, get all.
                query = r'SELECT * FROM {tn}'.format(tn = tabel)

            # Insert query into the new database
            conn.execute('INSERT INTO db2.{} '.format(tabel) + query.format(tn = tabel))

        # Disconnect
        conn.commit()
        conn.close()

def reset_hrdids(database, start_loc_id = 1, start_data_id = 1):
    """
    Function to reset the HRDLocationId's and HydroDynamicDataId's in a hrd.

    Parameters
    ----------
    database : str
        Path to database in which the correction is made
    start_loc_id : integer
        Start id for the HRDLocationId
    start_data_id : integer
        Start id for the HydroDynamicDataId

    Returns
    -------
    last_loc_id : integer
    last_data_id : integer
    """

    # Connect to database
    conn = sqlite3.connect(database)

    # Find table names
    table_names_query = "SELECT tbl_name FROM sqlite_master WHERE type='table'"
    table_names = conn.execute(table_names_query).fetchall()
    table_names = [naam[0] for naam in table_names]

    # Find all location ids
    loc_id_query = "SELECT HRDLocationId FROM HRDLocations"
    loc_ids_old = conn.execute(loc_id_query).fetchall()
    loc_ids_old = [loc_id[0] for loc_id in loc_ids_old]
    loc_ids_tmp = [-1 * loc_id for loc_id in loc_ids_old]
    loc_ids_new = range(start_loc_id, start_loc_id+len(loc_ids_old))

    # Find all data ids
    data_id_query = "SELECT HydroDynamicDataId FROM HydroDynamicData"
    data_ids_old = conn.execute(data_id_query).fetchall()
    data_ids_old = [data_id[0] for data_id in data_ids_old]
    data_ids_tmp = [-1 * data_id for data_id in data_ids_old]
    data_ids_new = range(start_data_id, start_data_id+len(data_ids_old))

    # Loop through tables and replace the ids
    for tabel in table_names:

        # Find all column namens
        columns = conn.execute('PRAGMA table_info({});'.format(tabel)).fetchall()
        columns = [col[1] for col in columns]

        # If the table has a column HRDLocationId, use the locatiequery
        if 'HRDLocationId' in columns:
            for loc_id_old, loc_id_tmp in zip(loc_ids_old, loc_ids_tmp):
                conn.execute('UPDATE {tn} SET HRDLocationId = ? WHERE HRDLocationId = ?'.format(tn=tabel), (loc_id_tmp, loc_id_old))
            for loc_id_tmp, loc_id_new in zip(loc_ids_tmp, loc_ids_new):
                conn.execute('UPDATE {tn} SET HRDLocationId = ? WHERE HRDLocationId = ?'.format(tn=tabel), (loc_id_new, loc_id_tmp))

        # Europoort contains a seiches table with LocationId instead of HRDLocationId
        elif 'LocationId' in columns:
            for loc_id_old, loc_id_tmp in zip(loc_ids_old, loc_ids_tmp):
                conn.execute('UPDATE {tn} SET LocationId = ? WHERE LocationId = ?'.format(tn=tabel), (loc_id_tmp, loc_id_old))
            for loc_id_tmp, loc_id_new in zip(loc_ids_tmp, loc_ids_new):
                conn.execute('UPDATE {tn} SET LocationId = ? WHERE LocationId = ?'.format(tn=tabel), (loc_id_new, loc_id_tmp))

        # Else if the table has a column HRDDynamicDataId, use the hydrodataquery
        if 'HydroDynamicDataId' in columns:
            for data_id_old, data_id_tmp in zip(data_ids_old, data_ids_tmp):
                conn.execute('UPDATE {tn} SET HydroDynamicDataId = ? WHERE HydroDynamicDataId = ?'.format(tn=tabel), (data_id_tmp, data_id_old))
            for data_id_tmp, data_id_new in zip(data_ids_tmp, data_ids_new):
                conn.execute('UPDATE {tn} SET HydroDynamicDataId = ? WHERE HydroDynamicDataId = ?'.format(tn=tabel), (data_id_new, data_id_tmp))

    # Commit changes
    conn.commit()

    # Close connection
    conn.close()

    # Return the last indices
    return loc_ids_new[-1], data_ids_new[-1]




def join_hrds(srcdatabase, adddatabase):
    """
    Function two add one database to another database with the same schema.

    Parameters
    ----------
    srcdatabase : str
        path to database to which a database is added
    adddatabase : str
        path to database which is added
    """

    # Connect to database
    conn = sqlite3.connect(adddatabase)
    # Query table names
    table_names_query = "SELECT tbl_name FROM sqlite_master WHERE type='table'"
    table_names = conn.execute(table_names_query).fetchall()
    table_names = [naam[0] for naam in table_names]

    # Attach the database
    conn.execute("ATTACH DATABASE ? AS srcdb", (srcdatabase,))

    # Loop through tables and copy a selection
    for tabel in table_names:
        columns = conn.execute('PRAGMA table_info({});'.format(tabel)).fetchall()
        columns = [col[1] for col in columns]
        checkcols = ['HRDLocationId', 'LocationId', 'HydroDynamicDataId']
        # If the table has a column HRDLocationId, LocationId or HydroDynamicDataId, add table
        if bool(set(checkcols) & set(columns)):
            # Insert query into the new database
            conn.execute('INSERT INTO srcdb.{} SELECT * FROM {}'.format(tabel, tabel))
    # Disconnect
    conn.commit()
    conn.execute("DETACH srcdb")
    conn.close()

def loc_to_file(conn, naam):
    """
    Function to export the loadcombinations of a location to one excel or csv.

    Parameters
    ----------
    conn : sqlite3.connection
        Connection to database from which the empty copy is made
    naam : str
        Locationname
    """

    if isinstance(naam, (str, int)):
        naam = [naam]

    # Zoek location id
    if isinstance(naam[0], str):
        locid = [find_HRDLocationId(n, conn) for n in naam]
    else:
        locid = naam


    query = ('SELECT * FROM HydroDynamicData WHERE HRDLocationId IN ({});'.format(','.join([str(i) for i in locid])))
    HydroDynamicData = pd.read_sql(query, conn, index_col='HydroDynamicDataId')
    
    # Inputdata
    hydrodynamicdataids = ','.join(HydroDynamicData.index.values.astype(str).tolist())
    query = """
    SELECT
    ID.HydroDynamicDataId,
    IV.ColumnName,
    ID.Value
    FROM HydroDynamicInputData ID
    INNER JOIN HRDInputVariables IV
    ON ID.HRDInputColumnId = IV.HRDInputColumnId
    WHERE HydroDynamicDataId IN ({});
    """.format(hydrodynamicdataids)

    HydroDynamicInputData = pd.read_sql(query, conn, index_col=['HydroDynamicDataId', 'ColumnName']).unstack()
    HydroDynamicInputData.columns = HydroDynamicInputData.columns.get_level_values(1)
    
    # Uitvoerdata
    query = """
    SELECT
    RD.HydroDynamicDataId,
    RV.ColumnName,
    RD.Value
    FROM HydroDynamicResultData RD
    INNER JOIN HRDResultVariables RV
    ON RD.HRDResultColumnId = RV.HRDResultColumnId
    WHERE HydroDynamicDataId IN ({});
    """.format(hydrodynamicdataids)

    HydroDynamicResultData = pd.read_sql(query, conn, index_col=['HydroDynamicDataId', 'ColumnName']).unstack()
    HydroDynamicResultData.columns = HydroDynamicResultData.columns.get_level_values(1)
    
    # Voeg samen
    resultaat = HydroDynamicData.join(HydroDynamicInputData).join(HydroDynamicResultData)
    
    return resultaat

def loc_to_file_JWS(srcdb, dstpath, naam, overwrite = False):
    """
    Function to export the loadcombinations of a location to one excel or csv.

    Parameters
    ----------
    srcdb : str
        path to database from which the empty copy is made
    dstpath : str
        path where to save the new empty database
    naam : str
        Locationname
    overwrite : boolean
        Whether to overwrite an existing file. Default is False.
    """
    # Definieer dictionaries met alle mogelijke stochasten aan de in- en uitvoerkant, gebaseerd op de HLCD
    raw_data_invoer = {'Stochastnr':      list(range(1,16)),
                       'Stochastnaam_NL': ["Afvoer Lobith","Afvoer Lith","Afvoer Borgharen","Afvoer Olst",
                                           "Afvoer Dalfsen","Waterstand Maasmond","Waterstand IJsselmeer",
                                           "Waterstand Markermeer", "Windsnelheid", "Waterstand", "Golfperiode",
                                           "Zeewaterstand", "Golfhoogte", "Zeewaterstand (u)", "Onzekerheid zeewaterstand (u)"],
                       'Stochastnaam_EN': ["Discharge Lobith","Discharge Lith","Discharge Borgharen","Discharge Olst",
                                           "Discharge Dalfsen","Water level Maasmond","Water level IJssel lake",
                                           "Water level Marker lake", "Wind speed", "Water level", "Wave period",
                                           "Sea water level", "Wave height", "Sea water level (u)", "Uncertainty water level (u)"],
                       'Eenheid':         ["m3/s","m3/s","m3/s","m3/s","m3/s","m3/s","m3/s", "m3/s", "m/s",
                                           "m+NAP", "s","m+NAP", "m", "m+NAP", "m"]}

    raw_data_uitvoer = {'Stochastnr':      list(range(1,8)),
                        'Stochastnaam_NL': ["Waterstand","Significante golfhoogte","Significante golfperiode (Ts)",
                                            "Piekgolfperiode (Tp)", "Gemiddelde piekgolfperiode (Tpm)",
                                            "Spectrale golfperiode (Tm-1,0)", "Golfrichting"],
                        'Stochastnaam_EN': ["Water level","Significant wave height","Significant wave period (Ts)",
                                            "Peak wave period (Tp)", "Average peak wave period (Tpm)",
                                            "Spectral wave period (Tm-1,0)", "Wave direction"],
                        'Eenheid':         ["m+NAP", "m", "s", "s", "s", "s", "degrees"]}

    stochasten_invoer  = pd.DataFrame(raw_data_invoer , columns = ['Stochastnr', 'Stochastnaam_NL', 'Stochastnaam_EN', 'Eenheid'])
    stochasten_uitvoer = pd.DataFrame(raw_data_uitvoer, columns = ['Stochastnr', 'Stochastnaam_NL', 'Stochastnaam_EN', 'Eenheid'])


    if not os.path.exists(srcdb):
        raise NameError('The system cannot find the path specified: \'{}\''.format(srcdb))

    # Delete the exisiting table if overwrite is True
    if dstpath:
        if os.path.exists(dstpath) and overwrite:
            os.remove(dstpath)
        elif os.path.exists(dstpath) and not overwrite:
            print('"{}" already exists. Set overwrite to True if you want to overwrite the existing database.'.format(dstpath))
            return None

    # Maak verbinding
    conn = sqlite3.connect(srcdb)

    # Zoek location id
    locid = find_HRDLocationId(naam, conn)

    query = ('SELECT * FROM HydroDynamicData WHERE HRDLocationId IN ({});'.format(int(locid)))
    HDD = pd.read_sql(query, conn)
    HDD.columns = ['data_id', 'locatie_id', 'sluitsituatie', 'windrichting']
    HDD.set_index(['data_id'], inplace = True)

    # Invoerdata
    query = 'SELECT * FROM {tn} WHERE HydroDynamicDataId IN (' + ','.join(['{}'.format(hdid) for hdid in HDD.index]) + ');'
    HDID = pd.read_sql(query.format(tn = 'HydroDynamicInputData'), conn)
    HDID.columns = ['data_id', 'invoer_id', 'invoer']
    HDID.set_index(['data_id', 'invoer_id'], inplace = True)
    HDID = HDID.unstack()
    HDID = HDD[['sluitsituatie', 'windrichting']].join(HDID)
    kolomnamen_invoer             = pd.read_sql('SELECT InputVariableId FROM HRDInputVariables', conn)
    Selectie_invoerStochasten     = kolomnamen_invoer.values.squeeze().tolist()
#    print(kolomnamen_invoer)
#    print(Selectie_invoerStochasten)
#    Gebruikte_invoerStochastnamen = []
#    for stochast in Selectie_invoerStochasten:
#        hulp  = stochasten_invoer.loc[stochasten_invoer['Stochastnr'] == stochast]
#        Gebruikte_invoerStochastnamen.append(hulp["Stochastnaam_NL"].values.squeeze().tolist())
#
#    HDID.columns = pd.MultiIndex.from_product((['invoer'], ['Sluitsituatie', 'Windrichting'] + Gebruikte_invoerStochastnamen))
    HDID.columns = pd.MultiIndex.from_product((['invoer'], ['Sluitsituatie', 'Windrichting'] + Selectie_invoerStochasten))

    # Uitvoerdata
    HDOD = pd.read_sql(query.format(tn = 'HydroDynamicResultData'), conn)
    HDOD.columns = ['data_id', 'uitvoer_id', 'uitvoer']
    HDOD.set_index(['data_id', 'uitvoer_id'], inplace = True)
    HDOD = HDOD.unstack()
    kolomnamen_uitvoer             = pd.read_sql('SELECT ResultVariableId FROM HRDResultVariables', conn)

    Selectie_uitvoerStochasten     = kolomnamen_uitvoer.values.squeeze().tolist()
#    Gebruikte_uitvoerStochastnamen = []
#    for stochast in Selectie_uitvoerStochasten:
#        hulp  = stochasten_uitvoer.loc[stochasten_uitvoer['Stochastnr'] == stochast]
#        Gebruikte_uitvoerStochastnamen.append(hulp["Stochastnaam_NL"].values.squeeze().tolist())
#
#    HDOD.columns = pd.MultiIndex.from_product((['uitvoer'], Gebruikte_uitvoerStochastnamen))
    try:
        HDOD.columns = pd.MultiIndex.from_product((['uitvoer'], Selectie_uitvoerStochasten))
    except:
        print('Could not replace output columns')

    # Voeg samen
    HDD = HDID.join(HDOD)

    conn.close()
    if dstpath:
        # Schrijf naar excel
        if dstpath.endswith('xlsx'):
            HDD.to_excel(dstpath)
        elif dstpath.endswith('csv'):
            HDD.to_csv(dstpath)
        else:
            raise TypeError('Exportformat of destination path \'{}\' not understood'.format(dstpath))
    else:
        return HDD

def add_to_table(tablename, dataframe, conn, ifexists='append'):

    # First check if table exists
    n = conn.execute("SELECT name FROM sqlite_master WHERE type='table' AND name=?;", (tablename,)).fetchall()
    if not np.size(n):
        raise ValueError('Table "{}" does not exist'.format(tablename))

    # Secondly check if all columns are present in the table
    columns = [r[1] for r in conn.execute("PRAGMA TABLE_INFO('{}');".format(tablename)).fetchall()]
    try:
        dataframe = dataframe[columns]
    except ValueError as err:
        print('The input dataframe does not have the right columns:', err.args)

    # Voeg toe aan tabel
    dataframe.to_sql(tablename, conn, if_exists=ifexists, chunksize=100000, index=False)

def bretschneider (d, fe, u):
    """
    Bereken golfcondities met Bretschneider.
    Gebaseerd op "subroutine Bretschneider" in Hydra-NL, geprogrammeerd door
    Matthijs Duits

    Parameters
    ----------
    d : float
        Waterdiepte
    fe : float
        Strijklengte
    u : float
        Windsnelheid

    Returns
    -------
    hs : float
        Significante golfhoogte
    tp : float
        Piekperiode
    """

    # Corrigeer als input is float
    if isinstance(d, float):
        if ((d <= 0.0) | (fe <= 0.0) | (u <= 0.0)):
            hs = 0.0
            tp = 0.0

    # Corrigeer als input is array
    else:
        hs_arr = np.empty_like(d)
        tp_arr = np.empty_like(d)

        indices = ((d <= 0.0) | (fe <= 0.0) | (u <= 0.0))
        hs_arr[indices] = 0.0
        tp_arr[indices] = 0.0

        u = u[~indices]
        fe = fe[~indices]
        d = d[~indices]

    # Initialiseer constanten
    g = 9.81

    # Bereken Hs en Tp
    dtilde = (d * g) / (u * u)
    v1 = np.tanh(0.53 * (dtilde ** 0.75))
    v2 = np.tanh(0.833 * (dtilde ** 0.375))

    ftilde = (fe * g) / (u * u)

    hhulp = (0.0125 / v1) * (ftilde ** 0.42)
    htilde = 0.283 * v1 * np.tanh(hhulp)
    hs = (u * u * htilde) / g

    thulp = (0.077 / v2) * (ftilde ** 0.25)
    ttilde = 2.4 * np.pi * v2 * np.tanh(thulp)
    tp = (1.08 * u * ttilde) / g

    if isinstance(d, float):
        return hs, tp

    else:
        hs_arr[~indices] = hs
        tp_arr[~indices] = tp
        return hs_arr, tp_arr

def vul_tabel_aan(tabel, controlekolom, aanvulkolom, controlewaarde, aanvulwaarde):
    """
    Functie om tabel aan te vullen met missende waarden. Ontworpen om
    ontbrekende windsnelheden op te vullen.

    Parameters
    ----------
    tabel : pandas.DataFrame
        Tabel die opgevuld moet worden.
    controlekolom : string
        Kolomnaam van de variabele waarvan de volledigheid wordt gecontroleerd.
    aanvulkolom : string
        Kolomnaam met variabele waarover de ontbrekende waarde gedupliceerd
        worden.
    controlewaarde : float
        Waarde waaraan de correcte hoeveelheid entries wordt gecontroleerd.
    aanvulwaarde : float
        Waarde waarbij de missende entries worden aangevuld.

    Returns
    -------
    tabel : pandas.DataFrame
        aangevulde tabel
    """

    # Bepaal eerst hoeveel gegevens er per aanwezige windrichting zijn, bij de controlewindsnelheid
    nentries = tabel.where(tabel[controlekolom] == controlewaarde).dropna(how='all').groupby(aanvulkolom).count()[controlekolom].to_dict()

    if len(tabel.where(tabel[controlekolom] == aanvulwaarde).dropna(how='all')[aanvulkolom].unique()) != 1:
        raise ValueError('Geen unieke waarde van aanvulkolom en aanvulwaarde (Bijvoorbeeld, meerdere windrichtingen aanwezig voor windsnelheid 0)')

    # Voeg vervolgens de tabellen toe aan op basis van het aantal entries bij 10.0 m/s
    # Richting = windrichgint, group = table waar ws = 0 bij wr
    for richting, group in tabel.where(tabel[controlekolom] == aanvulwaarde).dropna(how='all').groupby(aanvulkolom):
        
        # Voeg group toe aan elke windrichting met evenveel entries:
        for windrichting, aantal in nentries.items():
        
            # Als de richtingen overeenkomen, sla over want al aanwezig.
            if richting == windrichting:
                continue
            # Als de aantallen overeenkomen, verander windrichting en voeg toe.
            if aantal == len(group):
                group.loc[:, aanvulkolom] = windrichting
                tabel = pd.concat([tabel, group])

    return tabel

def vul_sqlite(conn, HRDLocations, HRDInputVariables, HRDResultVariables, resultaattabel, UncertaintyModelFactor, vacuum=True, verbose=True):
    """
    Vul sqlite-database met nieuwe gegevens van locaties, resultaten en
    onzekerheden. Deze tabellen worden uit de meegegeven database verwijderd,
    en vervolgens opnieuw gevuld met de meegegeven tabellen..

    Parameters
    ----------
    conn : sqlite3.connection
        Verbinding met de database die ingevuld moet worden. De opzet (PRAGMA)
        van de te vullen tabellen moet overeenkomen met de aangeleverde
        tabellen.
    HRDLocations : pandas.DataFrame
        Tabel met de locaties.
    HRDInputVariables : pandas.DataFrame
        Tabel met invoervariabelen.
    HRDResultVariables : pandas.DataFrame
        Tabel met uitvoervariabelen.
    resultaattabel : pandas.DataFrame
        Tabel met de resultaten. Deze tabel bevat de gegevens voor de drie
        tabellen HydroDynamicData, HydroDynamicInputData en HydroDynamicResultData.
        De kolommen moeten overeenkomen met de ColumnName's in de HRDInputVariables
        en HRDResultVariables (en dus ook van de verbonden database).
    UncertaintyModelFactor: pandas.DataFrame
        Tabel met de modelonzekerheden.
    vacuum : boolean
        Of de database achteraf gecompactiseerd moet worden (Default: True).
    verbose : boolean
        Of extra tekstuitvoer gegeven moet worden (default: True)
    """
    # Remove tables that need to be filled again
    for table in [
        'HydroDynamicData',
        'HydroDynamicInputData',
        'HydroDynamicResultData',
        'HRDLocations',
        'UncertaintyModelFactor',
        'HRDResultVariables',
        'HRDInputVariables'
    ]:
        conn.execute('DELETE FROM {table};'.format(table=table))

    # Add locations to database
    nlocations = len(HRDLocations.loc[:, 'HRDLocationId'].unique())
    if verbose:
        print('Adding HRDLocations ({} locations)...'.format(nlocations))
    add_to_table('HRDLocations', HRDLocations, conn)

    # If the result table does not contain a HRDLocationId column, the locations will be joined on X, Y
    # Round the x and y coordinates to the nearest integer, for better join.
    if 'HRDLocationId' not in resultaattabel.columns:
        for tabel in [HRDLocations, resultaattabel]:
            tabel.loc[:, ['XCoordinate', 'YCoordinate']] = tabel.loc[:, ['XCoordinate', 'YCoordinate']].round().astype(int)
        
        # Merge with locations on x and y and convert to HRDLocationId
        resultaattabel = resultaattabel.merge(
            HRDLocations[['HRDLocationId', 'XCoordinate', 'YCoordinate']],
            on=['XCoordinate', 'YCoordinate']).dropna().drop(['XCoordinate', 'YCoordinate'], axis=1)

    # HRDInputVariables
    if verbose:
        print('Adding HRDInputVariables')
    add_to_table('HRDInputVariables', HRDInputVariables, conn)

    # HRDResultVariables
    if verbose:
        print('Adding HRDResultVariables')
    add_to_table('HRDResultVariables', HRDResultVariables, conn)

    # HydroDynamicData
    HydroDynamicData = resultaattabel[['HydroDynamicDataId', 'HRDLocationId', 'ClosingSituationId', 'HRDWindDirectionId']]
    if nlocations != len(HydroDynamicData.loc[:, 'HRDLocationId'].unique()):
        raise ValueError("Not all locations are filled in database")
    if verbose:
        print('Adding HydroDynamicData ({} locations)...'.format(nlocations))
    add_to_table('HydroDynamicData', HydroDynamicData, conn)

    # HydroDynamicInputData
    HydroDynamicInputData = resultaattabel[['HydroDynamicDataId'] + HRDInputVariables['ColumnName'].tolist()]
    HydroDynamicInputData.columns = ['HydroDynamicDataId'] + HRDInputVariables['HRDInputColumnId'].tolist()
    HydroDynamicInputData.set_index('HydroDynamicDataId', inplace=True)
    HydroDynamicInputData = pd.DataFrame(HydroDynamicInputData.stack()).reset_index()
    HydroDynamicInputData.columns = ['HydroDynamicDataId', 'HRDInputColumnId', 'Value']

    if verbose:
        print('Adding HydroDynamicInputData...')
    add_to_table('HydroDynamicInputData', HydroDynamicInputData, conn)

    # HydroDynamicResultData
    HydroDynamicResultData = resultaattabel[['HydroDynamicDataId'] + HRDResultVariables['ColumnName'].tolist()]
    HydroDynamicResultData.columns = ['HydroDynamicDataId'] + HRDResultVariables['HRDResultColumnId'].tolist()
    HydroDynamicResultData.set_index('HydroDynamicDataId', inplace=True)
    HydroDynamicResultData = pd.DataFrame(HydroDynamicResultData.stack()).reset_index()
    HydroDynamicResultData.columns = ['HydroDynamicDataId', 'HRDResultColumnId', 'Value']

    add_to_table('HydroDynamicResultData', HydroDynamicResultData, conn)

    # UncertaintyModelFactor
    if verbose:
        print('Adding UncertaintyModelFactor...')
    add_to_table('UncertaintyModelFactor', UncertaintyModelFactor, conn)
    if nlocations != len(UncertaintyModelFactor.loc[:, 'HRDLocationId'].unique()):
        raise ValueError("Not all (or too much) locations are present in the uncertaintymodelfactor table")

    # Compact
    if vacuum:
        if verbose:
            print('Data added, vacuuming...')
        conn.execute('VACUUM;')

    conn.commit()

def pas_datum_aan(conn):
    # Pas datum aan
    now = datetime.datetime.now()
    conn.execute('UPDATE General SET CreationDate=?;', ('{:04d}-{:02d}-{:02d} {:02d}:{:02d}'.format(now.year, now.month, now.day, now.hour, now.minute), ))


def lees_swan(filelocation, U=None, R=None, M=None):
    """
    Leest SWAN-bestand in het formaat gebruikt voor WBI2023

    Parameters
    ----------
    filelocation : str
        Locatie waar het SWAN-bestand staat.
    U : int or float
        Windsnelheid, als dit argument niet aanwezig is wordt het afgeleid uit
        de bestandsnaam ( re.findall('U(.+)D(.+)L(.+)NZ', file) )
    R : int or float
        Windrichting, als dit argument niet aanwezig is wordt het afgeleid uit
        de bestandsnaam ( re.findall('U(.+)D(.+)L(.+)NZ', file) )
    M : int or float
        Zeewaterstand, als dit argument niet aanwezig is wordt het afgeleid uit
        de bestandsnaam ( re.findall('U(.+)D(.+)L(.+)NZ', file) )

    Returns
    -------
    resultaat : pd.DataFrame
        resultaattabel
    """

    # lees headerregel
    select_kolommen = ['Xp', 'Yp', 'Hsig', 'TPsmoo', 'Tm_10', 'Dir']
    with open(filelocation, 'r') as f:
        for line in f.readlines():
            if select_kolommen[0] in line:
                header = line[1:].split()
                break

    select_indices = [header.index(col) for col in select_kolommen]

    # Converteer de kolomnamen zo dat ze overeen komen met de algemene namen in de HRD's
    naam_conv = {
        'Xp' : 'XCoordinate',
        'Yp' : 'YCoordinate',
        'Hsig' : 'Hs',
        'TPsmoo' : 'Tp',
        'Tm_10' : 'Tm-1,0',
        'Dir' : 'Wave Direction'
    }
    kolomnamen = list(map(naam_conv.get, select_kolommen))

    # lees file
    resultaat = pd.read_table(filelocation, delim_whitespace=True, comment='%', header=None, usecols=select_indices, names=kolomnamen)

    # Parse bestandsnaam, als de belasting niet meegegeven is
    if (not U) and (not R) and (not M):
        file = os.path.split(filelocation)[1]
        [(ws, wr, zws)] = re.findall('U(.+)D(.+)L(.+)OO', file)
#    if not U:
#        U = float(ws)
#    if not R:
#        # De windrichting is afgerond, draai dit terug.
#        R = int(wr) // 22.5 * 22.5
#    if not M:
#        # De zeewaterstand is opgegeven in centimeters, maak hier meter van
#        M = float(M.replace('m', '-').replace('p', '')) / 100.

    # Voeg de belasting toe aan de resultaattabel
    resultaat.insert(2, 'HRDWindDirection', R)
    resultaat.insert(3, 'Wind speed', U)
    resultaat.insert(4, 'Water level OS11', M)

    return resultaat


def export_tables(hrdlocationids, hrd_database, as_dict=False):
    """
    Export HRDLocations, result table, UncertaintyModelFactor, HRDResultVariables, HRDInputVariables, HRDWindDirections from database

    Parameters
    ----------
    hrdlocationids : int, tuple, list
        HRDLocationIds for locations to export. If tuple, the nearest coordinate is searched for
    hrd_database : str
        Path to HRD
    """

    if isinstance(hrdlocationids, (int, np.integer)):
        hrdlocatieids = [hrdlocationids]

    if not isinstance(hrdlocationids, (tuple, list, np.array)):
        raise ValueError('Datatype "{}" for HRDLocationIds not understood.'.format(type(hrdlocationids)))

    conn = sqlite3.connect(hrd_database)

    if isinstance(hrdlocationids, tuple):
        # Find nearest location
        HRDLocations = pd.read_sql(
            'SELECT HRDLocationId, XCoordinate as X, YCoordinate AS Y FROM HRDLocations;',
            con=conn,
            index_col='HRDLocationId'
        )

        nearest = np.argmin(np.hypot(
            HRDLocations.Y.values - hrdlocationids[0], HRDLocations.Y.values - hrdlocationids[1]))
        
        hrdlocationids = [HRDLocations.index.values[nearest]]

    # HRDLocations
    hrdlocidstr = ','.join([str(i) for i in hrdlocationids])
    HRDLocations = pd.read_sql('SELECT * FROM HRDLocations WHERE HRDLocationId IN ({});'.format(hrdlocidstr), con=conn)
    
    # result table
    resultaattabel = loc_to_file(conn, hrdlocationids)

    # UncertaintyModelFactor
    UncertaintyModelFactor = pd.read_sql('SELECT * FROM UncertaintyModelFactor WHERE HRDLocationId IN ({});'.format(hrdlocidstr), con=conn)

    # HRDResultVariables
    HRDResultVariables = pd.read_sql('SELECT * FROM HRDResultVariables;', con=conn)

    # HRDInputVariables
    HRDInputVariables = pd.read_sql('SELECT * FROM HRDInputVariables;', con=conn)

    # HRDWindDirections
    winddirdict = pd.read_sql('SELECT Direction, HRDWindDirectionId FROM HRDWindDirections;', con=conn, index_col='Direction')['HRDWindDirectionId'].to_dict()

    if as_dict:
        return {
            'HRDLocations': HRDLocations,
            'resultaattabel': resultaattabel,
            'UncertaintyModelFactor': UncertaintyModelFactor,
            'HRDResultVariables': HRDResultVariables,
            'HRDInputVariables': HRDInputVariables,
            'winddirdict': winddirdict
        }
    else:
        return HRDLocations, resultaattabel, UncertaintyModelFactor,  HRDResultVariables, HRDInputVariables, winddirdict